'''
Visualize and compare kraken2 vs competitive mapping pipeline results for the 7 species
found in both databases, for both aeDNA-damaged and modern reads.
Uses compressed-pipeline outputs (chunk-split ngsLCA and bamdam files).

Produces three figures:
  1. kraken_summary.png        - read counts at each kraken2 stage per species
  2. competitive_summary.png   - read counts at each competitive mapping stage per species
  3. intersection.png          - bar plots showing read overlap between kraken2 and competitive

Usage:
  python pipeline_visualize_compressed.py \
    --species          data/species_in_both_dbs.txt \
    --kraken_dir       /path/to/kraken2/results_uniq \
    --kraken_filt_dir  /path/to/kraken2/filtered \
    --ngslca_dir       /path/to/ngslca_out_compressed \
    --ngslca_mod_dir   /path/to/modern_ngslca_out_compressed \
    --bamdam_dir       /path/to/bamdam_out_compressed \
    --bamdam_mod_dir   /path/to/modern_bamdam_out_compressed \
    --bamdam_filt_dir      /path/to/bamdam_filtered_compressed \
    --bamdam_mod_filt_dir  /path/to/modern_bamdam_filtered_compressed \
    --out_dir          figures_compressed/
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

MODE_LABEL = {
    'aeDNA':  'Simulated with deamination and mutation',
    'modern': 'Simulated without deamination or mutation',
}

def _fmt(v, _):
    v = int(v)
    if v >= 1_000_000: return f'{v/1_000_000:.1f}M'
    if v >= 1_000:     return f'{v/1_000:.0f}k'
    return str(v)


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
    # Try single merged file first
    path = os.path.join(bamdam_filt_dir, f'{species_id}{suffix}.shrunk.correct_genus.lca')
    if os.path.isfile(path):
        return [path]
    # Fall back to chunk files
    return sorted(glob.glob(os.path.join(bamdam_filt_dir, f'{species_id}{suffix}_chunk*.shrunk.correct_genus.lca')))


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


def read_names_bamdam(paths):
    '''Union read names across a list of .lca files.'''
    names = set()
    for path in paths:
        if path and os.path.isfile(path):
            with open(path) as f:
                for line in f:
                    if line.strip():
                        parts = line.split('\t')
                        if len(parts) >= 1:
                            names.add(parts[0].split(':')[0])
    return names


def count_lines(path):
    if path is None or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if line.strip():
                count += 1
    return count


def count_kraken_classified(path):
    if path is None or not os.path.isfile(path):
        return 0
    classified = 0
    with open(path) as f:
        for line in f:
            if line.strip() and line[0] == 'C':
                classified += 1
    return classified


def count_fasta_records(path):
    if path is None or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def count_ngslca_assigned(path):
    if path is None or not os.path.isfile(path):
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split('\t')
            if not parts[0] or parts[0].startswith('#'):
                continue
            if len(parts) >= 2 and parts[1].strip().split(':')[0] != '0':
                count += 1
    return count


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


def make_competitive_plot(species_ids, data, out_path, y_top=None):
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


def make_intersection_plot(species_ids, data, out_path):
    modes = ['aeDNA', 'modern']
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, mode in zip(axes, modes):
        x = list(range(len(species_ids)))
        k_only_vals, both_vals, b_only_vals = [], [], []

        for sid in species_ids:
            k_only_vals.append(data[sid][mode]['k_only'])
            both_vals.append(data[sid][mode]['overlap'])
            b_only_vals.append(data[sid][mode]['c_only'])

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

    fig.suptitle('K-mer-based vs local alignment-based: correct-genus read overlap per species',
                 fontsize=15, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


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
                writer.writerow({
                    'species_id':              sid,
                    'mode':                    mode,
                    'kraken_classified':       d['kraken']['classified'],
                    'kraken_genus':            d['kraken']['genus_level'],
                    'kraken_correct':          d['kraken']['correct_genus'],
                    'competitive_classified':  d['competitive']['classified'],
                    'competitive_genus':       d['competitive']['genus_level'],
                    'competitive_correct':     d['competitive']['correct_genus'],
                    'kraken_only':             d['k_only'],
                    'overlap':                 d['overlap'],
                    'competitive_only':        d['c_only'],
                })
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species',             required=True)
    parser.add_argument('--kraken_dir',          default=None)
    parser.add_argument('--kraken_filt_dir',     default=None)
    parser.add_argument('--ngslca_dir',          default=None, help='ngsLCA compressed output dir (aeDNA)')
    parser.add_argument('--ngslca_mod_dir',      default=None, help='ngsLCA compressed output dir (modern)')
    parser.add_argument('--bamdam_dir',          default=None, help='bamdam_out_compressed dir (aeDNA)')
    parser.add_argument('--bamdam_mod_dir',      default=None, help='bamdam_out_compressed dir (modern)')
    parser.add_argument('--bamdam_filt_dir',     default=None, help='bamdam_filtered_compressed dir (aeDNA)')
    parser.add_argument('--bamdam_mod_filt_dir', default=None, help='bamdam_filtered_compressed dir (modern)')
    parser.add_argument('--eprobe_input_dir',    default=None, help='eprobe_input dir; if set, aeDNA overlap is read from {sid}.fasta record counts')
    parser.add_argument('--modern_overlap_dir',  default=None, help='modern_overlap_counts dir; if set, modern overlap is read from {sid}.txt')
    parser.add_argument('--out_dir',             default='figures')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.species) as f:
        species_ids = [int(l.strip()) for l in f if l.strip()]

    data = {}
    for sid in species_ids:
        data[sid] = {}
        for modern, mode in [(False, 'aeDNA'), (True, 'modern')]:
            ngslca_dir  = args.ngslca_mod_dir      if modern else args.ngslca_dir
            bamdam_dir  = args.bamdam_mod_dir      if modern else args.bamdam_dir
            bamdam_filt = args.bamdam_mod_filt_dir if modern else args.bamdam_filt_dir

            k_out        = find_kraken_out(args.kraken_dir, sid, modern)
            k_classified = count_kraken_classified(k_out)
            k_genus      = count_lines(find_kraken_filtered(args.kraken_filt_dir, sid, modern, 'genus_level'))
            k_correct    = count_lines(find_kraken_filtered(args.kraken_filt_dir, sid, modern, 'correct_genus'))

            c_classified = sum(count_ngslca_assigned(p) for p in find_ngslca_lca(ngslca_dir, sid, modern))
            c_genus      = sum(count_lines(p) for p in find_bamdam_shrunk(bamdam_dir, sid, modern))
            c_correct    = sum(count_lines(p) for p in find_bamdam_filtered(bamdam_filt, sid, modern))

            if not modern and args.eprobe_input_dir:
                eprobe_path = os.path.join(args.eprobe_input_dir, f'{sid}.fasta')
                overlap = count_fasta_records(eprobe_path)
            elif modern and args.modern_overlap_dir:
                overlap_path = os.path.join(args.modern_overlap_dir, f'{sid}.txt')
                if os.path.isfile(overlap_path):
                    with open(overlap_path) as f:
                        overlap = int(f.read().strip())
                else:
                    overlap = 0
            else:
                k_names = read_names_kraken(find_kraken_filtered(args.kraken_filt_dir, sid, modern, 'correct_genus'))
                c_names = read_names_bamdam(find_bamdam_filtered(bamdam_filt, sid, modern))
                overlap = len(k_names & c_names)
                del k_names, c_names

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
                'k_only':  k_correct - overlap,
                'overlap': overlap,
                'c_only':  c_correct - overlap,
            }
            print(f"  {sid} [{mode}]: kraken_correct={k_correct} competitive_correct={c_correct} "
                  f"intersection={overlap}", flush=True)

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

    y_top = 2_000_000_000  # 2,000M ceiling

    make_kraken_plot(plotted_ids, data, os.path.join(args.out_dir, 'kraken_summary_compressed.png'), y_top=y_top)
    make_competitive_plot(plotted_ids, data, os.path.join(args.out_dir, 'competitive_summary_compressed.png'), y_top=y_top)
    make_intersection_plot(plotted_ids, data, os.path.join(args.out_dir, 'intersection_compressed.png'))
    write_summary_table(species_ids, data, os.path.join(args.out_dir, 'summary_table_compressed.csv'))


if __name__ == '__main__':
    main()
