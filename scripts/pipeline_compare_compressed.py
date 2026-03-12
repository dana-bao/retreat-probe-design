'''
Compare read counts across K-mer-based approach and global alignment-based approach pipeline stages
for species found in both databases, for both simulated-with-deamination and simulated-without-deamination reads.
Uses compressed-pipeline outputs (chunk-split ngsLCA and bamdam files).

For each species x mode (aeDNA/modern), counts reads at each stage:
  kraken_total           - all reads in K-mer-based approach .out (classified + unclassified)
  kraken_classified      - reads assigned by K-mer-based approach to any taxon (C lines)
  kraken_genus           - reads in kraken_filter genus_level.out
  kraken_correct         - reads in kraken_filter correct_genus.out
  competitive_classified - reads assigned by ngsLCA to any taxon (summed across chunks)
  competitive_genus      - reads in bamdam .shrunk.lca (summed across chunks)
  competitive_correct    - reads in bamdam_filter correct_genus .lca

Output: CSV with one row per (species_id, modern) combination.

Usage:
  python pipeline_compare_compressed.py \
    --species          data/species_in_both_dbs.txt \
    --kraken_dir       /path/to/kraken2/results_uniq \
    --kraken_filt_dir  /path/to/kraken2/filtered \
    --ngslca_dir       /path/to/ngslca_out_compressed \
    --ngslca_mod_dir   /path/to/modern_ngslca_out_compressed \
    --bamdam_dir       /path/to/bamdam_out_compressed \
    --bamdam_mod_dir   /path/to/modern_bamdam_out_compressed \
    --bamdam_filt_dir      /path/to/bamdam_filtered_compressed \
    --bamdam_mod_filt_dir  /path/to/modern_bamdam_filtered_compressed \
    --out                  pipeline_comparison_compressed.csv

All directory arguments are optional; missing files report as empty cells.

Variant of pipeline_compare.py, adapted for the compressed pipeline where some outputs are split into chunks and/or merged. See count_stats_compressed.py for annotation.
'''

import argparse
import csv
import glob
import os


def count_lines(path):
    if path is None or not os.path.isfile(path):
        return None
    count = 0
    with open(path) as f:
        for line in f:
            if line.strip():
                count += 1
    return count


def count_kraken_out(path):
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
    if path is None or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 2 and parts[1].strip() != '0':
                count += 1
    return count


def find_kraken_out(kraken_dir, species_id, modern):
    if kraken_dir is None:
        return None
    suffix = '_modern' if modern else ''
    pattern = os.path.join(kraken_dir, f'{species_id}{suffix}_task*.k2.0.2.core_nt.out')
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def find_kraken_filtered(kraken_filt_dir, species_id, modern, ftype):
    if kraken_filt_dir is None:
        return None
    suffix = '_modern' if modern else ''
    pattern = os.path.join(kraken_filt_dir, f'{species_id}{suffix}_task*.k2.0.2.core_nt.{ftype}.out')
    matches = glob.glob(pattern)
    return matches[0] if matches else None


def find_ngslca_lca(ngslca_dir, species_id, modern):
    '''Return list of chunk .lca files for this species.'''
    if ngslca_dir is None:
        return []
    suffix = '_modern' if modern else ''
    return sorted(glob.glob(os.path.join(ngslca_dir, f'{species_id}{suffix}_chunk*.lca')))


def find_bamdam_shrunk(bamdam_dir, species_id, modern):
    '''Return list of chunk .shrunk.lca files for this species.'''
    if bamdam_dir is None:
        return []
    suffix = '_modern' if modern else ''
    pattern = os.path.join(bamdam_dir, f'{species_id}{suffix}', f'{species_id}{suffix}_chunk*.shrunk.lca')
    return sorted(glob.glob(pattern))


def find_bamdam_filtered(bamdam_filt_dir, species_id, modern):
    '''Return list of correct_genus filtered .lca files (single merged or chunks).'''
    if bamdam_filt_dir is None:
        return []
    suffix = '_modern' if modern else ''
    path = os.path.join(bamdam_filt_dir, f'{species_id}{suffix}.shrunk.correct_genus.lca')
    if os.path.isfile(path):
        return [path]
    return sorted(glob.glob(os.path.join(bamdam_filt_dir, f'{species_id}{suffix}_chunk*.shrunk.correct_genus.lca')))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species',             required=True)
    parser.add_argument('--kraken_dir',          default=None)
    parser.add_argument('--kraken_filt_dir',     default=None)
    parser.add_argument('--ngslca_dir',          default=None)
    parser.add_argument('--ngslca_mod_dir',      default=None)
    parser.add_argument('--bamdam_dir',          default=None)
    parser.add_argument('--bamdam_mod_dir',      default=None)
    parser.add_argument('--bamdam_filt_dir',     default=None)
    parser.add_argument('--bamdam_mod_filt_dir', default=None)
    parser.add_argument('--out',                 default='pipeline_comparison_compressed.csv')
    args = parser.parse_args()

    with open(args.species) as f:
        species_ids = [int(line.strip()) for line in f if line.strip()]

    print(f"Processing {len(species_ids)} species x 2 modes = {len(species_ids)*2} rows", flush=True)

    rows = []
    for species_id in species_ids:
        for modern in (False, True):
            mode_label = 'modern' if modern else 'aeDNA'
            ngslca_dir  = args.ngslca_mod_dir      if modern else args.ngslca_dir
            bamdam_dir  = args.bamdam_mod_dir      if modern else args.bamdam_dir
            bamdam_filt = args.bamdam_mod_filt_dir if modern else args.bamdam_filt_dir

            print(f"  {species_id} [{mode_label}]", flush=True)

            k_out = find_kraken_out(args.kraken_dir, species_id, modern)
            k_total, k_classified = count_kraken_out(k_out)
            k_genus   = count_lines(find_kraken_filtered(args.kraken_filt_dir, species_id, modern, 'genus_level'))
            k_correct = count_lines(find_kraken_filtered(args.kraken_filt_dir, species_id, modern, 'correct_genus'))

            n_assigned = sum(count_ngslca_assigned(p) for p in find_ngslca_lca(ngslca_dir, species_id, modern))
            b_genus    = sum(count_lines(p) or 0 for p in find_bamdam_shrunk(bamdam_dir, species_id, modern))
            b_correct  = sum(count_lines(p) or 0 for p in find_bamdam_filtered(bamdam_filt, species_id, modern))

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

    print(f"\n{'species':>10} {'mode':>8} {'k_total':>10} {'k_correct':>10} {'c_correct':>10}")
    print('-' * 55)
    for r in rows:
        def fmt(v): return str(v) if v is not None else '-'
        print(f"{r['species_id']:>10} {r['mode']:>8} "
              f"{fmt(r['kraken_total']):>10} {fmt(r['kraken_correct']):>10} "
              f"{fmt(r['competitive_correct']):>10}")


if __name__ == '__main__':
    main()
