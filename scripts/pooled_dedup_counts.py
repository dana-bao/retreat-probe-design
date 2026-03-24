"""
Count probes per species in pooled_dedup FASTA and report alongside
species name and phylum, sorted by probe count descending.

Usage:
    python3 scripts/pooled_dedup_counts.py \
        --pooled_dedup  /rds/.../cdhit_output/pooled_dedup \
        --names         /rds/.../kraken2_taxonomy/names.dmp \
        --cache         data/taxonomy_cache.tsv \
        --out           data/pooled_dedup_counts.tsv
"""

import argparse
import csv
import os
from collections import defaultdict


def parse_species_counts(fasta_path):
    counts = defaultdict(int)
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                sid = line[1:].strip().split('|')[0]
                counts[sid] += 1
    return counts


def load_names(names_dmp):
    names = {}
    with open(names_dmp) as f:
        for line in f:
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 4 and parts[3] == 'scientific name':
                names[parts[0]] = parts[1]
    return names


def load_cache(path):
    cache = {}
    if os.path.isfile(path):
        with open(path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                cache[row['taxid']] = row['group']
    return cache


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pooled_dedup', required=True,
                    help='Path to pooled_dedup FASTA from cd-hit-est')
    ap.add_argument('--names',        required=True,
                    help='NCBI names.dmp for species name lookup')
    ap.add_argument('--cache',        default='data/taxonomy_cache.tsv',
                    help='TSV cache mapping taxid -> phylum group')
    ap.add_argument('--out',          default='data/pooled_dedup_counts.tsv')
    args = ap.parse_args()

    print(f'[1/3] Counting probes per species from {args.pooled_dedup}...')
    counts = parse_species_counts(args.pooled_dedup)
    print(f'      {sum(counts.values()):,} probes across {len(counts)} species')

    print(f'[2/3] Loading names and taxonomy cache...')
    names = load_names(args.names)
    cache = load_cache(args.cache)

    print(f'[3/3] Writing {args.out}...')
    rows = sorted(
        [{'species_id': sid, 'species_name': names.get(sid, sid),
          'phylum': cache.get(sid, 'Unclassified'), 'probe_count': count}
         for sid, count in counts.items()],
        key=lambda r: -r['probe_count']
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['species_id', 'species_name', 'phylum', 'probe_count'],
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    # also print to stdout
    print(f"\n{'species_id':<12} {'phylum':<18} {'probe_count':>12}  species_name")
    print('-' * 80)
    for r in rows:
        print(f"{r['species_id']:<12} {r['phylum']:<18} {r['probe_count']:>12,}  {r['species_name']}")
    print('-' * 80)
    print(f"{'TOTAL':<12} {'':<18} {sum(r['probe_count'] for r in rows):>12,}")
    print(f'\nSaved: {args.out}')


if __name__ == '__main__':
    main()
