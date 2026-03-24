"""
Correlation plot: bamdam correct rate vs Kraken2 correct rate
for compressed reads (pipeline_stats_compressed.csv).

Usage:
    python3 scripts/correct_rate_correlation.py \
        --input  data/pipeline_stats_compressed.csv \
        --out    figures/correct_rate_correlation_compressed.png
"""

import argparse
import csv
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='data/pipeline_stats_compressed.csv')
    ap.add_argument('--out',   default='figures/correct_rate_correlation_compressed.png')
    args = ap.parse_args()

    data = []
    with open(args.input) as f:
        for row in csv.DictReader(f):
            try:
                total = int(row['total_tiles'])
                bam   = int(row['bamdam_correct_genus'])
                krak  = int(row['kraken_correct_genus'])
                if total > 0 and bam > 0 and krak > 0:
                    data.append((bam / total, krak / total))
            except (ValueError, KeyError):
                pass

    x_bam  = np.array([d[0] for d in data])
    y_krak = np.array([d[1] for d in data])

    rho, pval = stats.spearmanr(x_bam, y_krak)

    # semi-log fit: log10(y) = m*x + b
    ly = np.log10(y_krak)
    m, b = np.polyfit(x_bam, ly, 1)
    x_range = np.linspace(0.9, 1.0, 300)

    n_outside = sum(1 for v in x_bam if v < 0.9)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x_bam, y_krak, s=22, alpha=0.65, color='#225ea8', linewidths=0)
    ax.plot(x_range, 10 ** (m * x_range + b), color='#d62728', linewidth=1.5,
            label=f'Semi-log fit (slope={m:.2f})')

    ax.set_yscale('log')
    ax.set_xlim(0.9, 1.002)
    ax.set_xlabel('bamdam correctly assigned / total tiles', fontsize=11)
    ax.set_ylabel('k-mer correct rate\n(Kraken2 correctly assigned / total tiles)', fontsize=11)
    ax.set_title('Extended Local Alignment Database,\nSimulated with Deamination and Mutation',
                 fontsize=12, fontweight='bold')

    ax.text(0.05, 0.95,
            f"Spearman r = {rho:.3f}\np = {pval:.2e}\nn = {len(x_bam)} ({n_outside} points outside view)",
            transform=ax.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc', alpha=0.8))
    ax.legend(fontsize=9, frameon=False, loc='lower left')
    ax.grid(True, which='both', linestyle='--', alpha=0.3)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {args.out}')


if __name__ == '__main__':
    main()
