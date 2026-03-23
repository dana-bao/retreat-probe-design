"""
Stacked bar plot: retreat species by taxonomic group, coloured by assembly availability.

Fetches higher-level taxonomy from NCBI Entrez (requires internet).
Results are cached to --cache so subsequent runs are instant.

Assembly availability categories (derived from all_availability.tsv):
  Species assembly    -- _merge == 'both'  (exact species-level assembly in our DB)
  Genus assembly only -- _merge == 'left_only' AND an Assembly Accession is present
                         (genus representative available but no species-level assembly)
  No assembly         -- _merge == 'left_only' AND no Assembly Accession

Usage:
  python scripts/taxonomy_availability_plot.py \\
      --input data/all_availability.tsv \\
      --out   figures/ \\
      --email your@email.com \\
      --rank  phylum
"""

import argparse
import collections
import csv
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio import Entrez

# availability categories
AVAIL_SPECIES = 'Species assembly'
AVAIL_GENUS   = 'Genus assembly only'
AVAIL_NONE    = 'No assembly'

# colors for plotting
COLORS = {
    AVAIL_SPECIES: '#7FC97F',   # mint
    AVAIL_GENUS:   '#BEAED4',   # lavender
    AVAIL_NONE:    '#FDC086',   # peach
}

CATEGORIES = [AVAIL_SPECIES, AVAIL_GENUS, AVAIL_NONE]


# data loading
def load_availability(path):
    """
    Read all_availability.tsv and return list of dicts with keys:
      taxid, name, availability
    Skips 'right_only' rows (strain/subspecies not in the primary retreat list).
    """
    rows = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            merge = row.get('_merge', '').strip()
            if merge == 'right_only':
                continue

            taxid = row['Species Taxonomic ID'].strip()
            name  = row['Species Name'].strip()
            has_accession = bool(row.get('Assembly Accession', '').strip())

            if merge == 'both':
                avail = AVAIL_SPECIES
            elif has_accession:
                avail = AVAIL_GENUS
            else:
                avail = AVAIL_NONE

            rows.append({'taxid': taxid, 'name': name, 'availability': avail})
    return rows


def load_strain_taxids(path):
    """
    Return a list of taxids for 'right_only' rows (strain/subspecies variants).
    These are excluded from the plot but their taxonomy is cached for other scripts.
    """
    taxids = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('_merge', '').strip() == 'right_only':
                taxids.append(row['Species Taxonomic ID'].strip())
    return taxids


# NCBI taxonomy lookup with caching
def fetch_group(taxid, rank, cache):
    """
    Return the scientific name at `rank` for `taxid` from NCBI taxonomy.
    Falls back to kingdom if `rank` is absent in the lineage.
    Results are stored in `cache` (dict taxid -> group name).
    """
    if taxid in cache:
        return cache[taxid]

    try:
        handle  = Entrez.efetch(db='taxonomy', id=taxid, retmode='xml')
        records = Entrez.read(handle)
        handle.close()
        time.sleep(0.35)   # stay within NCBI's 3 req/s limit

        lineage = records[0].get('LineageEx', [])

        for node in lineage:
            if node.get('Rank', '') == rank:
                cache[taxid] = node['ScientificName']
                return cache[taxid]

        # fallback: use kingdom label with "(other)" suffix
        for node in lineage:
            if node.get('Rank', '') == 'kingdom':
                cache[taxid] = node['ScientificName'] + ' (other)'
                return cache[taxid]

        cache[taxid] = 'Unclassified'
        return cache[taxid]

    except Exception as e:
        print(f'  WARNING: taxid {taxid} — {e}')
        cache[taxid] = 'Unclassified'
        return cache[taxid]


def load_cache(path):
    cache = {}
    if os.path.isfile(path):
        with open(path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                cache[row['taxid']] = row['group']
        print(f'Loaded {len(cache)} cached taxonomy entries from {path}')
    return cache


def save_cache(cache, path):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['taxid', 'group'], delimiter='\t')
        writer.writeheader()
        for taxid, group in cache.items():
            writer.writerow({'taxid': taxid, 'group': group})


# plotting
def make_plot(counts, groups, rank, out_path):
    """
    counts : dict  group -> Counter(availability -> n)
    groups : list of group names, sorted by total descending
    """
    x = range(len(groups))
    fig, ax = plt.subplots(figsize=(max(10, len(groups) * 1.1), 6))

    bottoms = [0] * len(groups)
    for cat in CATEGORIES:
        vals = [counts[g][cat] for g in groups]
        ax.bar(x, vals, bottom=bottoms, label=cat,
               color=COLORS[cat], edgecolor='white', linewidth=0.5)
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    for i, g in enumerate(groups):
        total = sum(counts[g].values())
        ax.text(i, bottoms[i] + 0.15, str(total),
                ha='center', va='bottom', fontsize=11)

    ax.set_xticks(list(x))
    ax.set_xticklabels(groups, rotation=40, ha='right', fontsize=12)
    ax.set_ylabel('Number of species', fontsize=13)
    ax.set_title(f'Retreat species by {rank} — assembly availability', fontsize=15)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=13)
    ax.set_ylim(0, max(bottoms) * 1.12)
    ax.yaxis.get_major_locator().set_params(integer=True)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


# main function
def main():
    ap = argparse.ArgumentParser(description='Taxonomic availability bar plot')
    ap.add_argument('--input',  default='data/all_availability.tsv',
                    help='Path to all_availability.tsv')
    ap.add_argument('--out',    default='figures',
                    help='Output directory for the figure')
    ap.add_argument('--email',  required=True,
                    help='Email address for NCBI Entrez (required by NCBI)')
    ap.add_argument('--rank',   default='phylum',
                    help='Taxonomic rank for grouping (default: phylum). '
                         'Try "class" for finer splits within Arthropoda/Tracheophyta.')
    ap.add_argument('--cache',  default='data/taxonomy_cache.tsv',
                    help='TSV cache for NCBI taxonomy lookups')
    ap.add_argument('--out-tsv', default=None,
                    help='Path to write summary TSV '
                         '(default: data/assembly_availability_by_<rank>.tsv)')
    args = ap.parse_args()

    Entrez.email = args.email
    os.makedirs(args.out, exist_ok=True)

    # load data
    species = load_availability(args.input)
    print(f'Loaded {len(species)} retreat species')

    strain_taxids = load_strain_taxids(args.input)
    print(f'Found {len(strain_taxids)} strain/subspecies variants (cached but excluded from plot)')

    # fetch taxonomy (with cache)
    cache = load_cache(args.cache)
    all_taxids = [sp['taxid'] for sp in species] + strain_taxids
    n_to_fetch = sum(1 for t in all_taxids if t not in cache)
    print(f'Fetching taxonomy for {n_to_fetch} taxids from NCBI '
          f'(rank={args.rank})...')

    for i, sp in enumerate(species):
        sp['group'] = fetch_group(sp['taxid'], args.rank, cache)
        if (i + 1) % 20 == 0 or (i + 1) == len(species):
            print(f'  {i + 1}/{len(species)}')

    for taxid in strain_taxids:
        fetch_group(taxid, args.rank, cache)

    save_cache(cache, args.cache)

    # aggregate counts
    counts = collections.defaultdict(collections.Counter)
    for sp in species:
        counts[sp['group']][sp['availability']] += 1

    # sort groups by total species count descending
    groups = sorted(counts, key=lambda g: -sum(counts[g].values()))

    print(f'\nGroups found ({args.rank} level):')
    for g in groups:
        c = counts[g]
        print(f'  {g}: {c[AVAIL_SPECIES]} species / '
              f'{c[AVAIL_GENUS]} genus-only / '
              f'{c[AVAIL_NONE]} none  '
              f'(total {sum(c.values())})')

    out_path = os.path.join(args.out, f'taxonomy_availability_{args.rank}.png')
    make_plot(counts, groups, args.rank, out_path)

    # write summary TSV
    tsv_path = args.out_tsv or f'data/assembly_availability_by_{args.rank}.tsv'
    with open(tsv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([args.rank, 'total',
                         'species_level', 'species_pct',
                         'genus_only',    'genus_pct',
                         'no_assembly',   'no_assembly_pct'])
        for g in groups:
            c     = counts[g]
            total = sum(c.values())
            s     = c[AVAIL_SPECIES]
            gn    = c[AVAIL_GENUS]
            n     = c[AVAIL_NONE]
            writer.writerow([
                g, total,
                s,  round(s  / total * 100, 1),
                gn, round(gn / total * 100, 1),
                n,  round(n  / total * 100, 1),
            ])
    print(f'Saved: {tsv_path}')


if __name__ == '__main__':
    main()
