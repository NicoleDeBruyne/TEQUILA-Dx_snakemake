#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2024.01.17
# Adapted from Robert Wang (Xing Lab)
# Optimized: 2025

import argparse
import pandas as pd
import numpy as np
import os
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Assesses whether a given gene shows unusually strong levels of allelic imbalance.")
    parser.add_argument('--infile', type=str, required=True, help='Input file with genes to analyze and haplotype assignments. Should contain columns: \
                                                                sample, gene, hap1_read_count, hap2_read_count, unassigned_read_count, max_coverage, and max_phased_coverage')
    parser.add_argument("--phasing-threshold", type=float, default=0.5, help="Minimum proportion of max coverage that phased to haplotypes for analysis")
    parser.add_argument('--sample-coverage-threshold', type=int, default=20, help='Minimum coverage required over the gene in the sample of interest for analysis')
    parser.add_argument('--padj-threshold', default=0.05, type=float, help='Maximum adjusted p-value to be considered an outlier. (default: 0.05)')
    parser.add_argument('--haplotype-ratio-threshold', default=0.1, type=float, help='Minimum difference in haplotype ratio to be considered an outlier. (default: 0.1)')
    parser.add_argument('--plot-volcano', action='store_true', help='Plot a volcano plot of difference in haplotype ratio vs. -log10(p-value) for all junctions.')
    parser.add_argument('--label-top-n-hits', default=0, type=int, help='Label the top N hits in the volcano plot.')
    parser.add_argument('--outprefix', required=True, type=str, help='Prefix for final output files (i.e. prefix_ase_results.tsv and prefix_ase_outliers.tsv)')
    return parser.parse_args()

def binomial_test(hap1_count, hap2_count):
    """Binomial test to check if hap1 and hap2 counts differ from the expected 50/50 ratio."""
    n = int(hap1_count) + int(hap2_count)
    minor_count = int(min(hap1_count, hap2_count))
    return stats.binomtest(k=minor_count, n=n, p=0.5, alternative='less').pvalue

def process_gene(hap1_count, hap2_count, max_coverage, max_phased_coverage,
                 phasing_threshold, sample_coverage_threshold):
    """Perform allele-specific expression (ASE) analysis for a given gene."""

    if np.isnan(max_phased_coverage) or (max_phased_coverage < sample_coverage_threshold) or (max_phased_coverage / max_coverage < phasing_threshold):
        return "n/a", "n/a", "n/a"

    total_hap = hap1_count + hap2_count
    ratio = min(hap1_count, hap2_count) / total_hap
    p_value = binomial_test(hap1_count, hap2_count)
    return ratio, float(ratio - 0.5), p_value

def plot_volcano(df, outprefix, padj_threshold, haplotype_ratio_threshold, n=0):
    """Plot a volcano plot of difference in haplotype ratio vs. -log10(p-value) for all junctions."""

    df = df.copy()
    df['padj'] = df['padj'].replace(0, np.nextafter(0, 1)).astype('float64')
    df['-log10_padj'] = -np.log10(df['padj'])
    df['diff'] = df['diff'].astype(float)

    df['rank_diff'] = df['diff'].abs().rank(ascending=False)
    df['rank_padj'] = df['-log10_padj'].rank(ascending=False)
    df['combined_score'] = df['rank_diff'] + df['rank_padj']
    top_hits = df.nsmallest(min(n, len(df)), 'combined_score')

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(df['diff'], df['-log10_padj'], color='black', s=1)
    ax.set_xlabel('Difference in haplotype ratio')
    ax.set_ylabel('-log10(p-value)')
    ax.axhline(-np.log10(padj_threshold), color='red', linestyle='--')
    ax.axvline(haplotype_ratio_threshold, color='red', linestyle='--')
    ax.axvline(-haplotype_ratio_threshold, color='red', linestyle='--')
    for _, row in top_hits.iterrows():
        ax.text(row['diff'], row['-log10_padj'], row['gene'],
                fontsize=6, color='red', ha='left', va='bottom', zorder=100, rotation=45)
    plt.savefig(f'{outprefix}_volcano_plot.png', dpi=300)
    plt.savefig(f'{outprefix}_volcano_plot.pdf')

def main():
    """Main function."""

    print(f"\n\n\n******************************************************************************************")
    print(f"Detecting ASE outliers...")
    print(f"******************************************************************************************\n")

    args = parse_args()

    outfile = args.outprefix + "_ase_results.tsv"
    outlier_outfile = args.outprefix + "_ase_outliers.tsv"
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)

    # Read input
    infile = pd.read_table(args.infile, sep="\t")

    print('Performing ASE analysis...')

    # Vectorise numeric coercions once, before the loop.
    # Read by column name (not position) so this stays correct regardless of
    # column order, and fails loudly if an expected column is missing.
    required_cols = ['sample', 'gene', 'hap1_read_count',
                      'hap2_read_count', 'unassigned_read_count',
                      'max_coverage', 'max_phased_coverage']
    missing = [c for c in required_cols if c not in infile.columns]
    if missing:
        raise ValueError(f"Input file {args.infile} is missing expected column(s): {missing}. "
                          f"Found columns: {list(infile.columns)}")

    hap1_counts = pd.to_numeric(infile['hap1_read_count'], errors='coerce')
    hap2_counts = pd.to_numeric(infile['hap2_read_count'], errors='coerce')
    unassigned_counts = pd.to_numeric(infile['unassigned_read_count'], errors='coerce')
    max_coverages = pd.to_numeric(infile['max_coverage'], errors='coerce')
    max_phased_coverages = pd.to_numeric(infile['max_phased_coverage'], errors='coerce')

    results = []
    for i in range(len(infile)):
        row = infile.iloc[i]
        sample = row['sample']
        gene = row['gene'].split(".")[0]
        hap1_count = hap1_counts.iloc[i]
        hap2_count = hap2_counts.iloc[i]
        unassigned_count = unassigned_counts.iloc[i]
        max_coverage = max_coverages.iloc[i]
        max_phased_coverage = max_phased_coverages.iloc[i]

        ratio, diff, p_value = process_gene(
            hap1_count, hap2_count, max_coverage, max_phased_coverage,
            args.phasing_threshold, args.sample_coverage_threshold
        )
        results.append([sample, gene, hap1_count, hap2_count, unassigned_count,
                         max_coverage, max_phased_coverage, ratio, diff, p_value])

    df = pd.DataFrame(results, columns=[
        'sample', 'gene', 'hap1_count', 'hap2_count', 'unassigned_count',
        'max_coverage', 'max_phased_coverage', 'ratio', 'diff', 'p_value'
    ])

    # Multiple-testing correction
    print('Correcting for multiple testing...')
    df['padj'] = 'n/a'
    df['padj'] = df['padj'].astype(object)
    mask = pd.to_numeric(df['p_value'], errors='coerce').notna()
    if mask.sum() == 0:
        print("No valid p-values to correct.")
    else:
        padj_values = multipletests(pd.to_numeric(df.loc[mask, 'p_value'], errors='coerce').dropna(), method='fdr_by')[1]
        df.loc[mask, 'padj'] = padj_values

    df.to_csv(outfile, sep='\t', index=False)
    print(f'ASE analysis results written to {outfile}')

    print('Identifying outliers based on user-defined thresholds...')
    outliers = df[
        (pd.to_numeric(df['padj'], errors='coerce') < args.padj_threshold) &
        (abs(pd.to_numeric(df['diff'], errors='coerce')) >= args.haplotype_ratio_threshold)
    ]
    outliers.to_csv(outlier_outfile, sep='\t', index=False)
    print(f'ASE outliers written to {outlier_outfile}')

    if args.plot_volcano:
        print('Plotting volcano plot...')
        plot_df = df.dropna(subset=['diff', 'padj'])
        plot_df = plot_df[plot_df['padj'] != 'n/a']
        plot_volcano(plot_df, args.outprefix, args.padj_threshold,
                     args.haplotype_ratio_threshold, args.label_top_n_hits)
        print(f'Volcano plot written to {args.outprefix}_volcano_plot.png/.pdf')

if __name__ == "__main__":
    main()
