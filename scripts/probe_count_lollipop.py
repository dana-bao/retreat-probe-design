"""
Manhattan-style lollipop plot: probe count per species in pooled_dedup,
grouped and coloured by phylum.

Usage:
    python3 probe_count_lollipop.py \
        --pooled_dedup  /rds/.../cdhit_output/pooled_dedup \
        --names         /rds/.../taxonomy/names.dmp \
        --cache         data/taxonomy_cache.tsv \
        [--email        your@email.com] \
        [--rank         phylum] \
        [--out          figures/probe_count_lollipop.png]

Sequence headers must be in the format >{species_id}|... (species ID before first |).
"""

import argparse
import csv
import os
import time
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

try:
    from Bio import Entrez
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# Same palette as pooled_dedup_piechart.py
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


def parse_species_counts(fasta_path):
    counts = defaultdict(int)
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                sid = line[1:].strip().split('|')[0]
                counts[sid] += 1
    return counts


def load_names(names_dmp):
    """Return dict taxid -> scientific name from NCBI names.dmp."""
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
        print(f'Loaded {len(cache)} cached taxonomy entries from {path}')
    return cache


def save_cache(cache, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
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


def abbreviate_name(name):
    """
    Shorten a species name for axis labels:
      'Silene acaulis'                            -> 'S. acaulis'
      'Anthyllis vulneraria subsp. alpestris'     -> 'A. vulneraria subsp. alp.'
      'Debaryomyces hansenii var. hansenii MTCC 234' -> 'D. hansenii var. hans. MTCC 234'
      'Debaryomyces hansenii CBS767'              -> 'D. hansenii CBS767'
    Rules:
      1. Genus -> first letter + '.'
      2. Species epithet kept in full
      3. subsp./var./f. marker kept, followed by first 4 chars of infraspecific name + '.'
      4. Strain codes (all-caps tokens or alphanumeric like CBS767) appended as-is
      5. Result capped at 26 chars (with trailing '…' if needed)
    """
    import re
    tokens = name.split()
    if not tokens:
        return name

    result = tokens[0][0] + '.'          # G.
    if len(tokens) > 1:
        result += ' ' + tokens[1]        # G. species

    i = 2
    while i < len(tokens):
        tok = tokens[i]
        if tok.lower() in ('subsp.', 'var.', 'f.', 'subsp', 'var', 'f'):
            marker = tok if tok.endswith('.') else tok + '.'
            if i + 1 < len(tokens):
                infra = tokens[i + 1]
                result += f' {marker} {infra[:4]}.'
                i += 2
            else:
                result += f' {marker}'
                i += 1
        elif re.match(r'^[A-Z]{2,}[0-9]*$', tok) or re.match(r'^[A-Z]+\d+$', tok):
            # strain code like CBS767, MTCC, ATCC
            strain = tok
            if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                strain += ' ' + tokens[i + 1]
                i += 1
            result += ' ' + strain
            i += 1
        else:
            i += 1

    if len(result) > 26:
        result = result[:25] + '…'
    return result


def make_plot(species_data, rank, out_path):
    """
    species_data: list of dicts with keys sid, name, phylum, count
    Sorted by phylum (alphabetically) then by count descending within phylum.
    """
    # sort: phylum alpha, then count descending within phylum
    species_data.sort(key=lambda r: (r['phylum'], -r['count']))

    # count species per phylum to replicate piechart color assignment:
    # largest phylum -> lightest color (same as pooled_dedup_piechart.py)
    from collections import Counter
    phylum_counts = Counter(r['phylum'] for r in species_data)
    phyla = [p for p, _ in phylum_counts.most_common()]  # largest first
    n = len(phyla)
    color_map = {p: PALETTE[(n - 1 - i) % len(PALETTE)] for i, p in enumerate(phyla)}

    fig, ax = plt.subplots(figsize=(max(12, len(species_data) * 0.18), 7))

    xs = np.arange(len(species_data))
    for i, r in enumerate(species_data):
        color = color_map[r['phylum']]
        ax.vlines(i, 1, r['count'], color=color, linewidth=1.2, zorder=2)
        ax.scatter(i, r['count'], color=color, s=25, zorder=3)

    # phylum background bands and labels
    phylum_spans = {}
    for i, r in enumerate(species_data):
        p = r['phylum']
        if p not in phylum_spans:
            phylum_spans[p] = [i, i]
        else:
            phylum_spans[p][1] = i

    for p, (x0, x1) in phylum_spans.items():
        color = color_map[p]
        ax.axvspan(x0 - 0.5, x1 + 0.5, alpha=0.07, color=color, zorder=0)
        pass

    ax.set_yscale('log')
    ax.set_xlim(-0.8, len(species_data) - 0.2)
    ax.set_xticks(xs)
    ax.set_xticklabels(
        [abbreviate_name(r['name']) for r in species_data],
        rotation=90, fontsize=6.5, ha='center'
    )
    ax.set_ylabel('Probe count (pooled_dedup)', fontsize=11)
    ax.set_title(f'Probes per species in pooled_dedup, grouped by {rank}',
                 fontsize=13, fontweight='bold')

    legend_patches = [mpatches.Patch(color=color_map[p], label=p) for p in phyla]
    ax.legend(handles=legend_patches, title=rank.capitalize(),
              bbox_to_anchor=(1.01, 1), loc='upper left',
              fontsize=8, title_fontsize=9, frameon=False)

    ax.grid(axis='y', linestyle='--', alpha=0.3)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    plt.subplots_adjust(bottom=0.25, right=0.82)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pooled_dedup', required=True,
                    help='Path to pooled_dedup FASTA from cd-hit-est')
    ap.add_argument('--names',        required=True,
                    help='NCBI names.dmp for species name lookup')
    ap.add_argument('--cache',        default='data/taxonomy_cache.tsv',
                    help='TSV cache mapping taxid -> phylum group')
    ap.add_argument('--email',        default=None,
                    help='Email for NCBI Entrez (required if taxids missing from cache)')
    ap.add_argument('--rank',         default='phylum',
                    help='Taxonomic rank for grouping (default: phylum)')
    ap.add_argument('--out',          default='figures/probe_count_lollipop.png')
    args = ap.parse_args()

    if HAS_BIOPYTHON and args.email:
        Entrez.email = args.email

    print(f'[1/4] Counting probes per species from {args.pooled_dedup}...')
    species_counts = parse_species_counts(args.pooled_dedup)
    print(f'      {sum(species_counts.values()):,} probes across {len(species_counts)} species')

    print(f'[2/4] Loading species names from {args.names}...')
    names = load_names(args.names)

    print(f'[3/4] Loading taxonomy cache from {args.cache}...')
    cache = load_cache(args.cache)
    missing = [sid for sid in species_counts if sid not in cache]
    if missing:
        if not args.email:
            print(f'[WARN] {len(missing)} species not in cache and no --email provided. '
                  f'Grouping as Unclassified.')
        else:
            print(f'      Fetching {len(missing)} missing entries from NCBI...')
        for sid in missing:
            fetch_group(sid, args.rank, cache)
        save_cache(cache, args.cache)

    print(f'[4/4] Generating lollipop plot...')
    species_data = [
        {
            'sid':    sid,
            'name':   names.get(sid, sid),
            'phylum': cache.get(sid, 'Unclassified'),
            'count':  count,
        }
        for sid, count in species_counts.items()
    ]
    make_plot(species_data, args.rank, args.out)


if __name__ == '__main__':
    main()
