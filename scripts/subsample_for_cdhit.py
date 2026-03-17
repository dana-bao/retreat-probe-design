"""
Subsample eprobe-filtered FASTA files based on a ratio threshold,
outputting to a deduplicate_input directory for downstream CD-HIT.

Usage:
    python3 subsample_for_cdhit.py \
        --plan         data/subsample_plan.csv \
        --filtered_dir /rds/.../eprobe_filtered \
        --out_dir      /rds/.../deduplicate_input \
        [--threshold   5.0] \
        [--seed        SEED]

For each species:
    if ratio > threshold -> subsample to 625 * threshold
    ratio <= threshold -> copy all sequences unchanged

Target is fixed at 625 for all species.
"""

import argparse
import csv
import glob
import os
import random
import sys


def read_fasta(path):
    records = []
    with open(path) as f:
        header = None
        seq = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header = line
                seq = []
            else:
                seq.append(line)
        if header is not None:
            records.append((header, ''.join(seq)))
    return records


def write_fasta(records, path):
    with open(path, 'w') as f:
        for header, seq in records:
            f.write(header + '\n')
            f.write(seq + '\n')


def find_filtered_fa(filtered_dir, species_id):
    matches = glob.glob(os.path.join(filtered_dir, f'{species_id}*.filtered.fa'))
    if len(matches) > 1:
        print(f'[WARN] multiple filtered files for {species_id}: {matches}', file=sys.stderr)
    return matches[0] if matches else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plan',         required=True, help='subsample_plan.csv')
    ap.add_argument('--filtered_dir', required=True, help='eprobe_filtered directory')
    ap.add_argument('--out_dir',      required=True, help='Output directory (deduplicate_input)')
    ap.add_argument('--threshold',    type=float, default=5.0,
                    help='Ratio threshold: species with ratio > threshold are subsampled '
                         'to round(target * threshold) (default: 5.0)')
    ap.add_argument('--seed',         type=int, default=None, help='Random seed for reproducibility (default: None, truly random)')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    if args.seed is not None:
        random.seed(args.seed)

    with open(args.plan) as f:
        plan = list(csv.DictReader(f))

    total_written = 0
    skipped = []

    TARGET = 625

    for row in plan:
        sid   = row['species_id']
        ratio = float(row['ratio'])

        fasta_path = find_filtered_fa(args.filtered_dir, sid)
        if fasta_path is None:
            print(f'[WARN] {sid}: no filtered FASTA found in {args.filtered_dir} — skipping',
                  file=sys.stderr)
            skipped.append(sid)
            continue

        records = read_fasta(fasta_path)
        n_avail = len(records)

        if ratio > args.threshold:
            subsample_to = round(TARGET * args.threshold)
            if n_avail <= subsample_to:
                print(f'[WARN] {sid}: only {n_avail} records available, '
                      f'wanted {subsample_to} — using all', file=sys.stderr)
                sampled = records
            else:
                sampled = random.sample(records, subsample_to)
            action = f'reduce->{subsample_to}'
        else:
            sampled = records
            action = 'keep_all'

        out_path = os.path.join(args.out_dir, f'{sid}.fasta')
        write_fasta(sampled, out_path)
        total_written += len(sampled)
        print(f'{sid} ({action}): {n_avail} -> {len(sampled)}  [{os.path.basename(fasta_path)}]')

    print(f'\nDone. Total sequences written: {total_written:,}')
    if skipped:
        print(f'Skipped ({len(skipped)}): {skipped}')


if __name__ == '__main__':
    main()
