'''
Visualize and compare kraken2 vs competitive mapping pipeline results for the 7 species
found in both databases, for both aeDNA-damaged and modern reads.

Produces three figures:
  1. kraken_summary.png   - read counts at each kraken2 stage per species
  2. competitive_summary.png   - read counts at each competitive mapping stage per species
  3. intersection.png     - bar plots showing read overlap between kraken2 correct_genus and competitive correct_genus per species 

Usage:
  python pipeline_visualize.py \
    --species          data/species_in_both_dbs.txt \
    --kraken_dir       /path/to/kraken2/results_uniq \
    --kraken_filt_dir  /path/to/kraken2/filtered \
    --ngslca_dir       /path/to/ngslca_out \
    --ngslca_mod_dir   /path/to/modern_ngslca_out \
    --bamdam_dir       /path/to/bamdam_out \
    --bamdam_mod_dir   /path/to/modern_bamdam_out \
    --bamdam_filt_dir      /path/to/bamdam_filtered \
    --bamdam_mod_filt_dir  /path/to/modern_bamdam_filtered \
    --out_dir          figures/
'''

import argparse
import csv
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# labels for plot
MODE_LABEL = {
    'aeDNA':  'Simulated with deamination and mutation',
    'modern': 'Simulated without deamination or mutation',
}

def _fmt(v, _):
    v = int(v)
    if v >= 1_000_000: return f'{v/1_000_000:.1f}M'
    if v >= 1_000:     return f'{v/1_000:.0f}k'
    return str(v)



# find relevant files based on species_id and mode (aeDNA or modern)
def find_kraken_out(kraken_dir, species_id, modern):
    if kraken_dir is None:
        return None
    suffix = '_modern' if modern else ''
    matches = glob.glob(os.path.join(kraken_dir, f'{species_id}{suffix}_task*.k2.0.2.core_nt.out'))
    return matches[0] if matches else None


def find_kraken_filtered(kraken_filt_dir, species_id, modern, ftype):
    if kraken_filt_dir is None:
        return None
    suffix = '_modern' if modern else ''
    matches = glob.glob(os.path.join(kraken_filt_dir, f'{species_id}{suffix}_task*.k2.0.2.core_nt.{ftype}.out'))
    return matches[0] if matches else None


def find_bamdam_shrunk(bamdam_dir, species_id, modern):
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
    if bamdam_filt_dir is None:
        return None
    suffix = '_modern' if modern else ''
    path = os.path.join(bamdam_filt_dir, f'{species_id}{suffix}.shrunk.correct_genus.lca')
    return path if os.path.isfile(path) else None


def find_ngslca_lca(ngslca_dir, species_id, modern):
    if ngslca_dir is None:
        return []
    suffix = '_modern' if modern else ''
    path = os.path.join(ngslca_dir, f'{species_id}{suffix}.lca')
    if os.path.isfile(path):
        return [path]
    chunks = sorted(glob.glob(os.path.join(ngslca_dir, f'{species_id}{suffix}_chunk*.lca')))
    return chunks


# acquiring read names within kraken and bamdam outputs for intersection analysis
def read_names_kraken(path):
    names = set()
    if path and os.path.isfile(path):
        with open(path) as f:
            for line in f:
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        names.add(parts[1])
    return names


def read_names_bamdam(path):
    names = set()
    if path and os.path.isfile(path):
        with open(path) as f:
            for line in f:
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 1:
                        names.add(parts[0].split(':')[0])
    return names


# count for number of lines within files
def count_lines(path):
    if path is None or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if line.strip():
                count += 1
    return count

# count for classified reads for kraken2
def count_kraken_classified(path):
    if path is None or not os.path.isfile(path):
        return 0
    classified = 0
    with open(path) as f:
        for line in f:
            if line.strip() and line[0] == 'C':
                classified += 1
    return classified


# count for assigned reads for ngsLCA
def count_ngslca_assigned(path):
    if path is None or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split('\t')
            # skip ngsLCA header lines (version line starts with tab, param line starts with #)
            if not parts[0] or parts[0].startswith('#'):
                continue
            if len(parts) >= 2 and parts[1].strip().split(':')[0] != '0':
                count += 1
    return count


# plotting function for k-mer-based approach, grouped bar chart per species, shared y-axis with log scale
def make_kraken_plot(species_ids, data, out_path, y_top=None):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    modes = ['aeDNA', 'modern']
    colors = {'classified': '#2B8CBE', 'genus_level': '#A6BDDB', 'correct_genus': '#ECE7F2'}
    labels = {'classified': 'Classified', 'genus_level': 'Genus level', 'correct_genus': 'Correct genus'}

    for ax, mode in zip(axes, modes):
        x = range(len(species_ids))
        width = 0.25
        for i, key in enumerate(['classified', 'genus_level', 'correct_genus']):
            vals = [data[sid][mode]['kraken'][key] for sid in species_ids]
            ax.bar([xi + i * width for xi in x], vals, width, label=labels[key], color=colors[key])

        ax.set_title(f'{MODE_LABEL[mode]}', fontsize=13)
        ax.set_xlabel('Species taxid', fontsize=13)
        ax.set_ylabel('Read count (log scale)', fontsize=13)
        ax.set_xticks([xi + width for xi in x])
        ax.set_xticklabels([str(s) for s in species_ids], rotation=0, ha='center', fontsize=12)
        ax.legend(fontsize=13)

    for ax in axes:
        ax.set_yscale('log')
        ax.set_ylim(bottom=1_000_000, top=y_top)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt))

    fig.suptitle('K-mer-based approach read classification per species', fontsize=15, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")

# plotting function for global alignment-based approach, same logic as above, shared y-axis with log scale
def make_competitive_plot(species_ids, data, out_path, y_top=None):
    '''
    Grouped bar chart: simulated with vs without deamination side by side per species.
    Bars: classified (ngsLCA), genus_level (bamdam shrunk), correct_genus (bamdam filtered).
    '''
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    modes = ['aeDNA', 'modern']
    colors = {'classified': '#2B8CBE', 'genus_level': '#A6BDDB', 'correct_genus': '#ECE7F2'}
    labels = {'classified': 'Classified', 'genus_level': 'Genus level', 'correct_genus': 'Correct genus'}

    for ax, mode in zip(axes, modes):
        x = range(len(species_ids))
        width = 0.25
        for i, key in enumerate(['classified', 'genus_level', 'correct_genus']):
            vals = [data[sid][mode]['competitive'][key] for sid in species_ids]
            ax.bar([xi + i * width for xi in x], vals, width, label=labels[key], color=colors[key])

        ax.set_title(f'{MODE_LABEL[mode]}', fontsize=13)
        ax.set_xlabel('Species taxid', fontsize=13)
        ax.set_ylabel('Read count (log scale)', fontsize=13)
        ax.set_xticks([xi + width for xi in x])
        ax.set_xticklabels([str(s) for s in species_ids], rotation=0, ha='center', fontsize=12)
        ax.legend(fontsize=13)

    for ax in axes:
        ax.set_yscale('log')
        ax.set_ylim(bottom=1_000_000, top=y_top)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt))

    fig.suptitle('Local alignment-based approach read classification per species', fontsize=15, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")

# plotting function for read overlap, stacked bar chart per species, shared y-axis with log scale
def make_intersection_plot(species_ids, data, out_path):
    '''
    Stacked bar per species: K-mer-based only | overlap | global alignment only.
    One subplot per mode (simulated with / without deamination).
    '''
    modes = ['aeDNA', 'modern']
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, mode in zip(axes, modes):
        x = list(range(len(species_ids)))
        k_only_vals, both_vals, b_only_vals = [], [], []

        for sid in species_ids:
            k_set = data[sid][mode]['kraken_names']
            b_set = data[sid][mode]['competitive_names']
            k_only_vals.append(len(k_set - b_set))
            both_vals.append(len(k_set & b_set))
            b_only_vals.append(len(b_set - k_set))

        ax.bar(x, k_only_vals, label='K-mer-based only',      color='#2C7FB8')
        ax.bar(x, both_vals,   label='Overlap',                color='#41B6C4',
               bottom=k_only_vals)
        ax.bar(x, b_only_vals, label='Local alignment only',  color='#7FCDBB',
               bottom=[k + b for k, b in zip(k_only_vals, both_vals)])

        ax.set_title(f'{MODE_LABEL[mode]}', fontsize=13)
        ax.set_xlabel('Species taxid', fontsize=13)
        ax.set_ylabel('Read count (log scale)', fontsize=13)
        ax.set_xticks(x)
        ax.set_xticklabels([str(s) for s in species_ids], rotation=0, ha='center', fontsize=12)
        ax.legend(fontsize=13)

    for ax in axes:
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt))

    fig.suptitle('K-mer-based vs local alignment-based: correct-genus read overlap per species (linear)',
                 fontsize=15, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")

# write summary table with all counts and overlaps for reference
def write_summary_table(species_ids, data, out_path):
    fields = [
        'species_id', 'mode',
        'kraken_classified', 'kraken_genus', 'kraken_correct',
        'competitive_classified', 'competitive_genus', 'competitive_correct',
        'kraken_only', 'overlap', 'competitive_only',
    ]
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for sid in species_ids:
            for mode in ['aeDNA', 'modern']:
                d = data[sid][mode]
                k_set = d['kraken_names']
                c_set = d['competitive_names']
                writer.writerow({
                    'species_id':              sid,
                    'mode':                    mode,
                    'kraken_classified':       d['kraken']['classified'],
                    'kraken_genus':            d['kraken']['genus_level'],
                    'kraken_correct':          d['kraken']['correct_genus'],
                    'competitive_classified':  d['competitive']['classified'],
                    'competitive_genus':       d['competitive']['genus_level'],
                    'competitive_correct':     d['competitive']['correct_genus'],
                    'kraken_only':             len(k_set - c_set),
                    'overlap':                 len(k_set & c_set),
                    'competitive_only':        len(c_set - k_set),
                })
    print(f"Saved: {out_path}")


# main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species',         required=True)
    parser.add_argument('--kraken_dir',      default=None)
    parser.add_argument('--kraken_filt_dir', default=None)
    parser.add_argument('--ngslca_dir',          default=None,  help='ngsLCA output dir (aeDNA)')
    parser.add_argument('--ngslca_mod_dir',      default=None,  help='ngsLCA output dir (modern)')
    parser.add_argument('--bamdam_dir',          default=None,  help='bamdam dir for aeDNA')
    parser.add_argument('--bamdam_mod_dir',      default=None,  help='bamdam dir for modern')
    parser.add_argument('--bamdam_filt_dir',     default=None,  help='bamdam_filter output dir (aeDNA)')
    parser.add_argument('--bamdam_mod_filt_dir', default=None,  help='bamdam_filter output dir (modern)')
    parser.add_argument('--out_dir',             default='figures')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.species) as f:
        species_ids = [int(l.strip()) for l in f if l.strip()]

    # collect data
    data = {}
    for sid in species_ids:
        data[sid] = {}
        for modern, mode in [(False, 'aeDNA'), (True, 'modern')]:
            ngslca_dir  = args.ngslca_mod_dir      if modern else args.ngslca_dir
            bamdam_dir  = args.bamdam_mod_dir      if modern else args.bamdam_dir
            bamdam_filt = args.bamdam_mod_filt_dir if modern else args.bamdam_filt_dir

            k_out = find_kraken_out(args.kraken_dir, sid, modern)
            k_classified = count_kraken_classified(k_out)
            k_genus   = count_lines(find_kraken_filtered(args.kraken_filt_dir, sid, modern, 'genus_level'))
            k_correct = count_lines(find_kraken_filtered(args.kraken_filt_dir, sid, modern, 'correct_genus'))

            c_classified = sum(count_ngslca_assigned(p) for p in find_ngslca_lca(ngslca_dir, sid, modern))
            c_genus   = count_lines(find_bamdam_shrunk(bamdam_dir, sid, modern))
            c_correct = count_lines(find_bamdam_filtered(bamdam_filt, sid, modern))

            k_names = read_names_kraken(find_kraken_filtered(args.kraken_filt_dir, sid, modern, 'correct_genus'))
            c_names = read_names_bamdam(find_bamdam_filtered(bamdam_filt, sid, modern))

            data[sid][mode] = {
                'kraken': {
                    'classified':    k_classified,
                    'genus_level':   k_genus,
                    'correct_genus': k_correct,
                },
                'competitive': {
                    'classified':    c_classified,
                    'genus_level':   c_genus,
                    'correct_genus': c_correct,
                },
                'kraken_names':      k_names,
                'competitive_names': c_names,
            }
            print(f"  {sid} [{mode}]: kraken_correct={k_correct} competitive_correct={c_correct} "
                  f"intersection={len(k_names & c_names)}", flush=True)

    # only plot species that have any available data
    plotted_ids = [
        sid for sid in species_ids
        if any(
            data[sid][mode]['kraken']['correct_genus'] > 0 or
            data[sid][mode]['competitive']['correct_genus'] > 0
            for mode in ['aeDNA', 'modern']
        )
    ]
    pending = [sid for sid in species_ids if sid not in plotted_ids]
    if pending:
        print(f"Note: {len(pending)} species have no data yet and are excluded from plots: {pending}")

    all_vals = [
        data[sid][mode][pipeline][key]
        for sid in plotted_ids
        for mode in ['aeDNA', 'modern']
        for pipeline in ['kraken', 'competitive']
        for key in ['classified', 'genus_level', 'correct_genus']
    ]
    y_top = max((v for v in all_vals if v > 0), default=1) * 30

    make_kraken_plot(plotted_ids, data, os.path.join(args.out_dir, 'kraken_summary.png'), y_top=y_top)
    make_competitive_plot(plotted_ids, data, os.path.join(args.out_dir, 'competitive_summary.png'), y_top=y_top)
    make_intersection_plot(plotted_ids, data, os.path.join(args.out_dir, 'intersection.png'))
    write_summary_table(species_ids, data, os.path.join(args.out_dir, 'summary_table.csv'))


if __name__ == '__main__':
    main()
