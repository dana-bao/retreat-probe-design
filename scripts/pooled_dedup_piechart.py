"""
Two pie charts for comparing species composition by taxonomic rank:
  1. Proportion of reads in pooled_dedup per phylum  -> pooled_dedup_piechart.png
  2. Proportion of species count per phylum           -> species_piechart.png

Both use the same color scheme and group ordering so the two charts are directly
comparable side-by-side.

Usage:
    python3 pooled_dedup_piechart.py \
        --pooled_dedup  /rds/.../pooled_dedup \
        --summary       data/eprobe_summary.csv \
        --cache         data/taxonomy_cache.tsv \
        [--email        your@email.com] \
        [--rank         phylum] \
        [--out_reads    figures/pooled_dedup_piechart.png] \
        [--out_species  figures/species_piechart.png]

Sequence headers must be in the format >{species_id}|... (species ID before first |).
If a species taxid is missing from the cache, NCBI is queried (--email required).
"""

import argparse
import csv
import os
import time
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from Bio import Entrez
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

PALETTE = [
    '#addd8e',  # YlGn  light green
    '#78c679',  # YlGn  medium green
    '#238443',  # YlGn  dark green
    '#004529',  # YlGn  very dark green
    '#253494',  # PuBu  dark navy
    '#225ea8',  # PuBu  medium-dark blue
    '#1d91c0',  # YlGnBu sky blue
    '#41b6c4',  # YlGnBu teal-blue
    '#7fcdbb',  # YlGnBu light teal
]


def assign_colors(labels_sorted_by_size):
    """Assign lightest colors to largest slices, darkest to smallest; cycle if >9 groups."""
    n = len(labels_sorted_by_size)
    return [PALETTE[(n - 1 - i) % len(PALETTE)] for i in range(n)]


def parse_species_counts(fasta_path):
    counts = defaultdict(int)
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                sid = line[1:].strip().split('|')[0]
                counts[sid] += 1
    return counts


def load_summary_species(summary_csv):
    """Return list of species IDs from eprobe_summary.csv."""
    sids = []
    with open(summary_csv) as f:
        for row in csv.DictReader(f):
            sids.append(row['species_id'].strip())
    return sids


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


def fetch_group(taxid, rank, cache):
    if taxid in cache:
        return cache[taxid]
    if not HAS_BIOPYTHON:
        cache[taxid] = 'Unclassified'
        return cache[taxid]
    try:
        handle  = Entrez.efetch(db='taxonomy', id=taxid, retmode='xml')
        records = Entrez.read(handle)
        handle.close()
        time.sleep(0.35)
        lineage = records[0].get('LineageEx', [])
        for node in lineage:
            if node.get('Rank', '') == rank:
                cache[taxid] = node['ScientificName']
                return cache[taxid]
        for node in lineage:
            if node.get('Rank', '') == 'kingdom':
                cache[taxid] = node['ScientificName'] + ' (other)'
                return cache[taxid]
        cache[taxid] = 'Unclassified'
    except Exception as e:
        print(f'  WARNING: taxid {taxid} — {e}')
        cache[taxid] = 'Unclassified'
    return cache[taxid]


def make_pie(group_data, title, out_path, rank, colors_by_group):
    """
    group_data : list of (group_label, count) sorted largest first
    colors_by_group : dict mapping group label -> hex color
    """
    labels = [g for g, _ in group_data]
    sizes  = [n for _, n in group_data]
    total  = sum(sizes)
    colors = [colors_by_group[g] for g in labels]

    fig, ax = plt.subplots(figsize=(10, 8))
    wedges, texts = ax.pie(
        sizes,
        labels=None,
        colors=colors,
        startangle=140,
        wedgeprops=dict(linewidth=0.5, edgecolor='white'),
    )

    for i, (wedge, size) in enumerate(zip(wedges, sizes)):
        pct = size / total * 100
        if pct < 1.0:
            continue
        angle = np.radians((wedge.theta1 + wedge.theta2) / 2)
        if pct >= 2.5:
            # text inside the slice
            x = 0.80 * np.cos(angle)
            y = 0.80 * np.sin(angle)
            text_color = 'white' if colors[i] in PALETTE[:5] else '#333333'
            ax.text(x, y, f'{pct:.1f}%', ha='center', va='center',
                    fontsize=9, color=text_color)
        else:
            # arrow pointing from outside into the slice
            x_tip  = 0.88 * np.cos(angle)
            y_tip  = 0.88 * np.sin(angle)
            x_text = 1.30 * np.cos(angle)
            y_text = 1.30 * np.sin(angle)
            ax.annotate(
                f'{pct:.1f}%',
                xy=(x_tip, y_tip), xytext=(x_text, y_text),
                ha='center', va='center', fontsize=8, color='#333333',
                arrowprops=dict(arrowstyle='->', color='#666666', lw=0.7),
            )

    legend_labels = [f'{lab}  ({n:,})' for lab, n in zip(labels, sizes)]
    ax.legend(
        wedges, legend_labels,
        title=f'{rank.capitalize()}',
        loc='center left',
        bbox_to_anchor=(1.0, 0.5),
        fontsize=9,
        title_fontsize=10,
    )

    ax.set_title(f'{title}\n(total = {total:,})', fontsize=13, fontweight='bold', pad=16)

    plt.tight_layout()
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pooled_dedup', required=True,
                    help='Path to pooled_dedup FASTA')
    ap.add_argument('--summary',      required=True,
                    help='data/eprobe_summary.csv (used for species-count pie chart)')
    ap.add_argument('--cache',        default='data/taxonomy_cache.tsv',
                    help='TSV cache from taxonomy_availability_plot.py')
    ap.add_argument('--email',        default=None,
                    help='Email for NCBI Entrez (required if taxids missing from cache)')
    ap.add_argument('--rank',         default='phylum',
                    help='Taxonomic rank for grouping (default: phylum)')
    ap.add_argument('--out_reads',    default='figures/pooled_dedup_piechart.png')
    ap.add_argument('--out_species',  default='figures/species_piechart.png')
    args = ap.parse_args()

    if HAS_BIOPYTHON and args.email:
        Entrez.email = args.email

    # ---- reads per species from pooled_dedup ----
    species_read_counts = parse_species_counts(args.pooled_dedup)
    print(f'Total reads: {sum(species_read_counts.values()):,} across {len(species_read_counts)} species')

    # ---- species list from eprobe_summary ----
    all_species = load_summary_species(args.summary)
    print(f'Total species in pipeline: {len(all_species)}')

    # ---- taxonomy cache ----
    cache = load_cache(args.cache)

    all_sids = set(species_read_counts.keys()) | set(all_species)
    missing  = [sid for sid in all_sids if sid not in cache]
    if missing:
        if not args.email:
            print(f'[WARN] {len(missing)} species not in cache and no --email provided. '
                  f'Grouping as Unclassified: {missing}')
        else:
            print(f'Fetching taxonomy for {len(missing)} species from NCBI...')
        for sid in missing:
            fetch_group(sid, args.rank, cache)
        save_cache(cache, args.cache)

    # ---- aggregate reads by group ----
    reads_by_group = defaultdict(int)
    for sid, n in species_read_counts.items():
        reads_by_group[cache.get(sid, 'Unclassified')] += n

    # ---- aggregate species count by group ----
    species_by_group = defaultdict(int)
    for sid in all_species:
        species_by_group[cache.get(sid, 'Unclassified')] += 1

    # ---- determine shared group order (by species count, high to low) ----
    reads_sorted   = sorted(reads_by_group.items(),   key=lambda x: x[1], reverse=True)
    species_sorted = sorted(species_by_group.items(), key=lambda x: x[1], reverse=True)

    # canonical order = species-count order; reads chart uses same wedge positions
    species_order = [g for g, _ in species_sorted]
    reads_reordered = [(g, reads_by_group.get(g, 0)) for g in species_order]
    # append any groups present only in reads (not in species list)
    species_set = set(species_order)
    for g, n in reads_sorted:
        if g not in species_set:
            reads_reordered.append((g, n))

    # assign colors based on species order so both charts share the same mapping
    # largest slices (by species count) get lighter colors, smallest get darker
    n_groups = len(species_order)
    colors_by_group = {g: PALETTE[(n_groups - 1 - i) % len(PALETTE)] for i, g in enumerate(species_order)}
    extra_i = n_groups
    for g, _ in reads_reordered:
        if g not in colors_by_group:
            colors_by_group[g] = PALETTE[extra_i % len(PALETTE)]
            extra_i += 1

    make_pie(reads_reordered, f'Reads in pooled_dedup by {args.rank}',
             args.out_reads,   args.rank, colors_by_group)
    make_pie(species_sorted,  f'Pipeline species by {args.rank}',
             args.out_species, args.rank, colors_by_group)


if __name__ == '__main__':
    main()
