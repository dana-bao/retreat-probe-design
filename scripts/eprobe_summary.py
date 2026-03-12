"""
Summarize eProbe input sequence counts and biophysical filter results per species.

Usage:
    python3 scripts/eprobe_summary.py \
        --eprobe_input_dir    /rds/.../eprobe_input \
        --eprobe_filtered_dir /rds/.../eprobe_filtered \
        --out                 data/eprobe_summary.csv
"""

import argparse
import csv
import os

# count sequences in fasta
def count_fasta_seqs(path):
    if not path or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count

# count lines in lca
def count_lca_lines(path):
    if not path or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if line.strip():
                count += 1
    return count

# handle different input suffixes
def parse_sid(sid):
    for suffix, itype in [
        ('_bamdam_sub300k', 'bamdam_sub300k'),
        ('_bamdam',         'bamdam'),
        ('_sub300k',        'sub300k'),
    ]:
        if sid.endswith(suffix):
            return sid[:-len(suffix)], itype
    return sid, 'plain'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--eprobe_input_dir',    required=True, help='eprobe_input directory')
    parser.add_argument('--eprobe_filtered_dir', required=True, help='eprobe_filtered directory')
    parser.add_argument('--out',                 default='eprobe_summary.csv')
    args = parser.parse_args()

    # collect species from eprobe filtered output
    sids = sorted(
        f.replace('.filtered.fa', '')
        for f in os.listdir(args.eprobe_filtered_dir)
        if f.endswith('.filtered.fa')
    )

    rows = []
    for sid in sids:
        numeric_id, input_type = parse_sid(sid)
        input_path    = os.path.join(args.eprobe_input_dir,    f'{sid}.fasta')
        filtered_path = os.path.join(args.eprobe_filtered_dir, f'{sid}.filtered.fa')
        rejected_path = os.path.join(args.eprobe_filtered_dir, f'{sid}.rejected.fa')

        eprobe_input    = count_fasta_seqs(input_path)
        eprobe_passed   = count_fasta_seqs(filtered_path)
        eprobe_rejected = count_fasta_seqs(rejected_path)
        pct_passed = f'{100 * eprobe_passed / eprobe_input:.1f}' if eprobe_input else 'N/A'

        rows.append({
            'species_id':      numeric_id,
            'input_type':      input_type,
            'eprobe_input':    eprobe_input,
            'eprobe_passed':   eprobe_passed,
            'eprobe_rejected': eprobe_rejected,
            'pct_passed':      pct_passed,
        })
        print(f'{numeric_id} ({input_type}): input={eprobe_input:,}  passed={eprobe_passed:,}  '
              f'rejected={eprobe_rejected:,}  ({pct_passed}%)')

    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['species_id', 'input_type', 'eprobe_input',
                                               'eprobe_passed', 'eprobe_rejected', 'pct_passed'])
        writer.writeheader()
        writer.writerows(rows)

    print(f'\nSaved: {args.out}')


if __name__ == '__main__':
    main()
