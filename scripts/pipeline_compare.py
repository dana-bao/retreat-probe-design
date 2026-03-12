'''
Compare read counts across K-mer-based approach and global alignment-based approach pipeline stages
for species found in both databases, for both simulated-with-deamination and simulated-without-deamination reads.

For each species x mode (aeDNA/modern), counts reads at each stage:
  kraken_total           - all reads in K-mer-based approach .out (classified + unclassified)
  kraken_classified      - reads assigned by K-mer-based approach to any taxon (C lines)
  kraken_genus           - reads in kraken_filter genus_level.out
  kraken_correct         - reads in kraken_filter correct_genus.out
  competitive_classified - reads assigned by ngsLCA to any taxon (taxid != 0 in .lca)
  competitive_genus      - reads in bamdam .shrunk.lca (genus-filtered by bamdam)
  competitive_correct    - reads in bamdam_filter correct_genus .lca

Output: CSV with one row per (species_id, modern) combination.

Usage:
  python pipeline_compare.py \
    --species          data/species_in_both_dbs.txt \
    --kraken_dir       /path/to/kraken2/results_uniq \
    --kraken_filt_dir  /path/to/kraken2/filtered \
    --ngslca_dir       /path/to/ngslca_out \
    --ngslca_mod_dir   /path/to/modern_ngslca_out \
    --bamdam_dir       /path/to/bamdam_out \
    --bamdam_mod_dir   /path/to/modern_bamdam_out \
    --bamdam_filt_dir      /path/to/bamdam_filtered \
    --bamdam_mod_filt_dir  /path/to/modern_bamdam_filtered \
    --out                  pipeline_comparison.csv

All directory arguments are optional; missing files report as empty cells.
'''

import argparse
import csv
import glob
import os
import sys


# counting functions
def count_lines(path):
    '''Count non-empty lines in a file. Returns None if file missing.'''
    if path is None or not os.path.isfile(path):
        return None
    count = 0
    with open(path) as f:
        for line in f:
            if line.strip():
                count += 1
    return count


def count_kraken_out(path):
    '''Count total and classified (C) reads in a Kraken2 .out file.'''
    if path is None or not os.path.isfile(path):
        return None, None
    total = 0
    classified = 0
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            total += 1
            if line[0] == 'C':
                classified += 1
    return total, classified


def count_ngslca_assigned(path):
    '''Count reads assigned by ngsLCA (column 2 taxid != 0).'''
    if path is None or not os.path.isfile(path):
        return None
    count = 0
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 2 and parts[1].strip() != '0':
                count += 1
    return count


# file discovery functions
def find_kraken_out(kraken_dir, species_id, modern):
    '''Find {species_id}[_modern]_task*.k2.0.2.core_nt.out'''
    if kraken_dir is None:
        return None
    suffix = '_modern' if modern else ''
    pattern = os.path.join(kraken_dir, f'{species_id}{suffix}_task*.k2.0.2.core_nt.out')
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def find_kraken_filtered(kraken_filt_dir, species_id, modern, ftype):
    '''Find {species_id}[_modern]_task*.k2.0.2.core_nt.{ftype}.out'''
    if kraken_filt_dir is None:
        return None
    suffix = '_modern' if modern else ''
    pattern = os.path.join(kraken_filt_dir, f'{species_id}{suffix}_task*.k2.0.2.core_nt.{ftype}.out')
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def find_ngslca_lca(ngslca_dir, species_id, modern):
    '''Find {species_id}[_modern].lca directly in ngslca_dir.'''
    if ngslca_dir is None:
        return None
    suffix = '_modern' if modern else ''
    path = os.path.join(ngslca_dir, f'{species_id}{suffix}.lca')
    return path if os.path.isfile(path) else None


def find_bamdam_shrunk(bamdam_dir, species_id, modern):
    '''Find {species_id}[_modern].shrunk.lca under bamdam_dir (one subdir deep).'''
    if bamdam_dir is None:
        return None
    suffix = '_modern' if modern else ''
    fname = f'{species_id}{suffix}.shrunk.lca'
    for candidate in [
        os.path.join(bamdam_dir, fname),
        os.path.join(bamdam_dir, f'{species_id}{suffix}', fname),
    ]:
        if os.path.isfile(candidate):
            return candidate
    return None


def find_bamdam_filtered(bamdam_filt_dir, species_id, modern):
    '''Find {species_id}[_modern].shrunk.correct_genus.lca in bamdam_filt_dir.'''
    if bamdam_filt_dir is None:
        return None
    suffix = '_modern' if modern else ''
    path = os.path.join(bamdam_filt_dir, f'{species_id}{suffix}.shrunk.correct_genus.lca')
    return path if os.path.isfile(path) else None


# main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species',         required=True,  help='File with one species taxid per line')
    parser.add_argument('--kraken_dir',      default=None,   help='Kraken2 results_uniq directory')
    parser.add_argument('--kraken_filt_dir', default=None,   help='kraken_filter output directory')
    parser.add_argument('--ngslca_dir',      default=None,   help='ngsLCA output directory (regular)')
    parser.add_argument('--ngslca_mod_dir',  default=None,   help='ngsLCA output directory (modern)')
    parser.add_argument('--bamdam_dir',      default=None,   help='bamdam output directory (regular)')
    parser.add_argument('--bamdam_mod_dir',  default=None,   help='bamdam output directory (modern)')
    parser.add_argument('--bamdam_filt_dir',     default=None, help='bamdam_filter output directory (regular)')
    parser.add_argument('--bamdam_mod_filt_dir', default=None, help='bamdam_filter output directory (modern)')
    parser.add_argument('--out',             default='pipeline_comparison.csv')
    args = parser.parse_args()

    with open(args.species) as f:
        species_ids = [int(line.strip()) for line in f if line.strip()]

    print(f"Processing {len(species_ids)} species x 2 modes = {len(species_ids)*2} rows", flush=True)

    rows = []
    for species_id in species_ids:
        for modern in (False, True):
            mode_label = 'modern' if modern else 'aeDNA'
            ngslca_dir    = args.ngslca_mod_dir     if modern else args.ngslca_dir
            bamdam_dir    = args.bamdam_mod_dir     if modern else args.bamdam_dir
            bamdam_filt   = args.bamdam_mod_filt_dir if modern else args.bamdam_filt_dir

            print(f"  {species_id} [{mode_label}]", flush=True)

            # k-mer based
            k_out = find_kraken_out(args.kraken_dir, species_id, modern)
            k_total, k_classified = count_kraken_out(k_out)
            k_genus  = count_lines(find_kraken_filtered(args.kraken_filt_dir, species_id, modern, 'genus_level'))
            k_correct = count_lines(find_kraken_filtered(args.kraken_filt_dir, species_id, modern, 'correct_genus'))

            # global alignment based
            n_assigned = count_ngslca_assigned(find_ngslca_lca(ngslca_dir, species_id, modern))
            b_genus    = count_lines(find_bamdam_shrunk(bamdam_dir, species_id, modern))
            b_correct  = count_lines(find_bamdam_filtered(bamdam_filt, species_id, modern))

            rows.append({
                'species_id':             species_id,
                'mode':                   mode_label,
                'kraken_total':           k_total,
                'kraken_classified':      k_classified,
                'kraken_genus':           k_genus,
                'kraken_correct':         k_correct,
                'competitive_classified': n_assigned,
                'competitive_genus':      b_genus,
                'competitive_correct':    b_correct,
            })

    fieldnames = [
        'species_id', 'mode',
        'kraken_total', 'kraken_classified', 'kraken_genus', 'kraken_correct',
        'competitive_classified', 'competitive_genus', 'competitive_correct',
    ]

    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults written to: {args.out}", flush=True)

    # summary
    print(f"\n{'species':>10} {'mode':>8} {'k_total':>10} {'k_correct':>10} {'c_correct':>10}")
    print('-' * 55)
    for r in rows:
        def fmt(v): return str(v) if v is not None else '-'
        print(f"{r['species_id']:>10} {r['mode']:>8} "
              f"{fmt(r['kraken_total']):>10} {fmt(r['kraken_correct']):>10} "
              f"{fmt(r['competitive_correct']):>10}")


if __name__ == '__main__':
    main()
