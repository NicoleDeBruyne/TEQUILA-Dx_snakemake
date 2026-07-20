#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2024.01.17
# Adapted from Robert Wang (Xing Lab)
# Optimized: 2025

import os, argparse, warnings
import pandas as pd
import numpy as np
from pandas.errors import PerformanceWarning
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore', category=PerformanceWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Filters output from perform_splice_junction_beta_binomial_tests.py to identify '
                    'splice junctions with unusually high or low usage frequencies in the sample of '
                    'interest based on user-supplied padj and deltaPSI thresholds.')
    parser.add_argument('--infile', required=True,
        help='Path to input file (output from perform_splice_junction_beta_binomial_tests.py)')
    parser.add_argument('--outfile', required=True, type=str,
        help='Output file to write outlier junctions.')
    parser.add_argument('--padj-threshold', default=0.05, type=float,
        help='Maximum adjusted p-value to be considered an outlier. (default: 0.05)')
    parser.add_argument('--delta-PSI-threshold', default=0.1, type=float,
        help='Minimum abs(delta PSI) to be considered an outlier. (default: 0.1)')
    parser.add_argument('--plot-volcano', action='store_true',
        help='Plot a volcano plot of delta PSI vs. -log10(p-value) for all junctions.')
    parser.add_argument('--label-top-n-hits', default=0, type=int,
        help='Label the top N hits in the volcano plot.')
    parser.add_argument('--threads', type=int, default=1,
        help='Number of threads to use for parallel processing. Default: 1')
    return parser.parse_args()


def plot_volcano(df, outfile, padj_threshold, delta_PSI_threshold, n=0):
    """Plot a volcano plot of delta PSI vs. -log10(p-value) for all junctions."""

    df = df.copy()
    df['padj'] = df['padj'].replace(0, np.nextafter(0, 1)).astype('float64')
    df['-log10_padj'] = -np.log10(df['padj'])
    df['delta_PSI'] = df['delta_PSI'].astype(float)

    df['rank_delta_PSI'] = df['delta_PSI'].abs().rank(ascending=False)
    df['rank_padj'] = df['-log10_padj'].rank(ascending=False)
    df['combined_score'] = df['rank_delta_PSI'] + df['rank_padj']
    top_hits = df.nsmallest(min(n, len(df)), 'combined_score')

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df['delta_PSI'], df['-log10_padj'], color='black', s=4)
    ax.set_xlabel('Delta PSI')
    ax.set_ylabel('-log10(padj)')
    ax.axhline(-np.log10(padj_threshold), color='red', linestyle='--')
    ax.axvline(delta_PSI_threshold, color='red', linestyle='--')
    ax.axvline(-delta_PSI_threshold, color='red', linestyle='--')
    for _, row in top_hits.iterrows():
        ax.text(row['delta_PSI'], row['-log10_padj'], row['gene'],
                fontsize=6, color='red', ha='left', va='bottom', zorder=100, rotation=45)

    plot_outfile = f'{os.path.splitext(outfile)[0]}_volcano_plot'
    plt.savefig(f'{plot_outfile}_volcano_plot.pdf')


def main():
    """Main script."""

    print("\n\n\n******************************************************************************************")
    print("Identifying splice junction outliers...")
    print("******************************************************************************************\n")

    args = parse_args()

    if not os.path.exists(args.infile):
        raise FileNotFoundError(f"Input file {args.infile} not found.")
    df = pd.read_csv(args.infile, sep='\t')

    outdir = os.path.dirname(args.outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    print('Identifying outliers based on user-defined thresholds...')
    padj_numeric = pd.to_numeric(df['padj'], errors='coerce')
    dpsi_numeric = pd.to_numeric(df['delta_PSI'], errors='coerce')
    outliers = df[(padj_numeric < args.padj_threshold) & (dpsi_numeric.abs() >= args.delta_PSI_threshold)]

    outliers.to_csv(args.outfile, sep='\t', index=False)
    print(f'Outliers written to {args.outfile}')

    if args.plot_volcano:
        print('Plotting volcano plot...')
        plot_df = df.dropna(subset=['delta_PSI', 'padj'])
        plot_df = plot_df[plot_df['padj'] != 'n/a']
        plot_volcano(plot_df, args.outfile, args.padj_threshold,
                     args.delta_PSI_threshold, args.label_top_n_hits)
        print(f'Volcano plot written to {args.outfile}_volcano_plot.pdf')


if __name__ == '__main__':
    main()
