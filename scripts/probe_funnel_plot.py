"""
Sankey-style funnel plot showing read retention across the probe design pipeline:
    eprobe input -> biophysical filter -> subsample -> CD-HIT dedup

Each bar represents 100% of reads entering that stage. The coloured portion
shows what fraction passes; the grey portion shows what is dropped.
Bezier flows between bars contract from 100% on the left to the retained
fraction on the right, making the drop visible in every transition.

Usage:
    python3 probe_funnel_plot.py \
        --summary     data/eprobe_summary.csv \
        --plan        data/subsample_plan.csv \
        --post_dedup  100357 \
        [--threshold  5.0] \
        [--out        figures/probe_funnel.png]
"""

import argparse
import csv
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path


def load_totals(summary_csv, plan_csv, threshold):
    eprobe_input_total  = 0
    eprobe_passed_total = 0
    with open(summary_csv) as f:
        for row in csv.DictReader(f):
            eprobe_input_total  += int(row['eprobe_input'])
            eprobe_passed_total += int(row['eprobe_passed'])

    subsampled_total = 0
    with open(plan_csv) as f:
        for row in csv.DictReader(f):
            passed = int(row['eprobe_passed'])
            ratio  = float(row['ratio'])
            if ratio > threshold:
                subsampled_total += min(passed, round(625 * threshold))
            else:
                subsampled_total += passed

    return eprobe_input_total, eprobe_passed_total, subsampled_total


def bezier_flow(ax, x1, x2, y_bot_l, y_top_l, y_bot_r, y_top_r, color, alpha):
    """Draw a Bezier-connected filled region between two x positions."""
    cx = (x1 + x2) / 2
    verts = [
        (x1, y_bot_l),
        (cx, y_bot_l), (cx, y_bot_r), (x2, y_bot_r),
        (x2, y_top_r),
        (cx, y_top_r), (cx, y_top_l), (x1, y_top_l),
        (x1, y_bot_l),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    ax.add_patch(mpatches.PathPatch(
        Path(verts, codes), facecolor=color, edgecolor='none', alpha=alpha, zorder=2
    ))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--summary',    required=True, help='eprobe_summary.csv')
    ap.add_argument('--plan',       required=True, help='subsample_plan.csv')
    ap.add_argument('--post_dedup', required=True, type=int,
                    help='Read count after CD-HIT deduplication')
    ap.add_argument('--threshold',  type=float, default=5.0,
                    help='Subsample ratio threshold used (default: 5.0)')
    ap.add_argument('--out',        default='figures/probe_funnel.png')
    args = ap.parse_args()

    eprobe_input, eprobe_passed, subsampled = load_totals(
        args.summary, args.plan, args.threshold
    )
    post_dedup = args.post_dedup

    counts = [eprobe_input, eprobe_passed, subsampled, post_dedup]

    # stage-to-stage retention percentages (each relative to previous stage)
    pcts = [
        100.0,
        eprobe_passed / eprobe_input  * 100,
        subsampled    / eprobe_passed * 100,
        post_dedup    / subsampled    * 100,
    ]

    stage_labels = [
        'Eprobe\ninput',
        'Biophysical\nfilter',
        f'Subsample\n({args.threshold:.0f}× target)',
        'CD-HIT\ndedup',
    ]
    drop_labels = [
        None,
        'failed biophysical filter',
        'removed by subsampling',
        'removed by deduplication',
    ]

    # Paul Tol muted palette
    bar_colors = ['#4477AA', '#66CCEE', '#44BB99', '#CCBB44']
    drop_color = '#DDDDDD'

    bar_w = 0.32
    x_pos = [0.0, 1.3, 2.6, 3.9]

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.set_xlim(-0.45, x_pos[-1] + 0.45)
    ax.set_ylim(-22, 120)

    for i in range(len(counts)):
        x     = x_pos[i]
        pct   = pcts[i]
        color = bar_colors[i]

        # coloured retained bar
        ax.bar(x, pct, width=bar_w, color=color, zorder=3, linewidth=0)
        # grey dropped region above
        if pct < 100:
            ax.bar(x, 100 - pct, width=bar_w, bottom=pct,
                   color=drop_color, zorder=3, linewidth=0)

        # bezier flow to next bar:
        # coloured portion of bar i (0→pct) expands to fill bar i+1 (0→100)
        # grey portion of bar i (pct→100) collapses upward to a point at y=100, vanishing at the top
        # the two regions share a clean non-crossing boundary throughout
        if i < len(counts) - 1:
            xn       = x_pos[i + 1]

            # coloured flow: left matches coloured bar i (0→pct), expands to full bar i+1 (0→100)
            bezier_flow(ax, x + bar_w/2, xn - bar_w/2,
                        0, pct, 0, 100,
                        color=color, alpha=0.38)

            # grey flow: left matches grey portion of bar i (pct→100), collapses to point at y=100
            if pct < 100:
                bezier_flow(ax, x + bar_w/2, xn - bar_w/2,
                            pct, 100, 100, 100,
                            color=drop_color, alpha=0.45)

        # percentage label: inside bar if tall enough, just above bar top if small
        label_pct = f'{pct:.1f}%' if i > 0 else '100%\n(baseline)'
        if pct >= 20:
            ax.text(x, pct / 2, label_pct,
                    ha='center', va='center', fontsize=11,
                    fontweight='bold', color='white', zorder=5)
        else:
            # bar too narrow — place label just above the coloured bar
            ax.text(x, pct + 1.5, label_pct,
                    ha='center', va='bottom', fontsize=10,
                    fontweight='bold', color=color, zorder=5)

        # drop percentage label inside grey region
        if i > 0 and pct < 95:
            drop_pct = 100 - pct
            mid_grey = pct + drop_pct / 2
            ax.text(x, mid_grey,
                    f'−{drop_pct:.1f}%\n{drop_labels[i]}',
                    ha='center', va='center', fontsize=8.5,
                    color='#777777', style='italic', zorder=5)

        # absolute count above bar — larger and bold
        ax.text(x, 102, f'{counts[i]:,}',
                ha='center', va='bottom', fontsize=11,
                fontweight='bold', color='#333333', zorder=5)

        # stage label below
        ax.text(x, -4, stage_labels[i],
                ha='center', va='top', fontsize=11,
                fontweight='bold', color='#333333')

    ax.set_ylabel('Percentage retained from previous stage (%)', fontsize=12)
    ax.set_xticks([])
    for spine in ['top', 'right', 'bottom']:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {args.out}')


if __name__ == '__main__':
    main()
