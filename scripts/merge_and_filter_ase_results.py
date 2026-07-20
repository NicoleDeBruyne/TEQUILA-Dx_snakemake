#!/usr/bin/env python

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2024.10.23

import argparse
import pandas as pd
import os
from math import ceil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

def parse_args():
    """ Parse command line arguments """
    
    parser = argparse.ArgumentParser(description="Remove genes with outlier ASE from long-read RNA-seq data that are present in a large number of samples (i.e. from the same TEQUILA panel and/or sequencing batch).")
    parser.add_argument("--infiles", nargs="+", help="Input TSV files containing the outlier genes to be filtered (output from run_ase_analysis.py).", required=True)
    parser.add_argument("--outprefix", help="Prefix for output files.", required=True)
    parser.add_argument("--title", default="Number of Genes with Allele-specific Expression by Sample", help="Title for the plot.")
    parser.add_argument('--min-haplotype-ratio', default=0, type=float, help='Exclude genes with haplotype ratio < min-haplotype-ratio, for example, if you believe these are false positives.')
    parser.add_argument('--delta-haplotype-ratio-threshold', default=0.1, type=float, help='Minimum difference in haplotype ratio to be considered an outlier. (default: 0.1)')
    parser.add_argument('--padj-threshold', default=0.05, type=float, help='Maximum adjusted p-value to be considered an outlier. (default: 0.05)')
    parser.add_argument("--sample-number-threshold", type=int, help="The maximum number of samples that a gene can be an outlier in to be retained.")
    parser.add_argument("--plot", action="store_true", help="If set, generates a plot of the number of genes with outlier ASE according to the haplotype-ratio and padj thresholds before and after sample number filtering.")
    
    return parser.parse_args()

def plot_outlier_counts(samples, dfs, legends, suptitle, outfile):
    """Generate bar plots and box plots."""

    # Count unique genes per sample for each dataset
    counts_list = []
    for df in dfs:
        counts_list.append(df.groupby('sample')['gene'].nunique().sort_index())

    # Ensure all samples are represented in each series in the same order
    sample_order = sorted(samples)
    for i, counts in enumerate(counts_list):
        counts_list[i] = counts.reindex(sample_order, fill_value=0)

    # ---------------- Bar plots ----------------
    barplot_file = os.path.splitext(outfile)[0] + "_barplot.pdf"
    fig, axes = plt.subplots(len(dfs), 1, figsize=(max(16, len(sample_order)*0.2), max(6, len(dfs)*4)), sharex=True)
    if len(dfs) == 1:
        axes = [axes]  # ensure iterable

    for ax, counts, (color, label) in zip(axes, counts_list, legends):
        ax.bar(counts.index, counts.values, color=color)
        ax.set_title(f'Average={counts.mean():.2f}, Median={counts.median()}, Range={counts.min()}–{counts.max()}')
        ax.set_ylabel('Count')
        ax.tick_params(axis='x', rotation=90)
        ax.legend(handles=[mpatches.Patch(color=color, label=label)],
                  loc='center left', bbox_to_anchor=(1.0, 0.5))

    axes[-1].set_xlabel('Sample')
    fig.suptitle(suptitle, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(barplot_file)
    plt.close()
    print(f"Bar plot saved to {barplot_file}")

    # ---------------- Box plot ----------------
    boxplot_file = os.path.splitext(outfile)[0] + "_boxplot.pdf"
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Data for boxplot: one box per dataset
    data = [counts.values for counts in counts_list]
    colors = [c for c, _ in legends]

    # Draw empty boxes with black edges and black median
    box = ax.boxplot(
        data,
        patch_artist=False,
        medianprops=dict(color='black'),
        boxprops=dict(color='black'),
        whiskerprops=dict(color='black'),
        capprops=dict(color='black'),
        showfliers=False
    )

    # Plot points with jitter and dataset colors
    jitter_strength = 0.08
    for i, counts in enumerate(data):
        x = np.full(len(counts), i + 1) + np.random.uniform(-jitter_strength, jitter_strength, len(counts))
        ax.scatter(x, counts, color=colors[i], alpha=0.7, s=20)

    # Clean x-axis (no labels, just spacing)
    ax.set_xticks(range(1, len(dfs) + 1))
    ax.set_xticklabels([])

    ax.set_ylabel('Number of genes')
    ax.set_title(suptitle)

    # Legend with stats included, with stats on new line
    legend_handles = []
    for counts, color, (_, label) in zip(counts_list, colors, legends):
        stats_text = f"{label}\nAvg={counts.mean():.2f}, Med={counts.median()}, Range={counts.min()}–{counts.max()}"
        legend_handles.append(mpatches.Patch(color=color, label=stats_text))
    ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1.05, 0.5))

    plt.tight_layout()
    plt.savefig(boxplot_file, bbox_inches='tight')
    plt.close()
    print(f"Box plot saved to {boxplot_file}")

def main():
    """ Main function """

    # Parse command line arguments
    args = parse_args()
    
    # Check that at least one input file was provided and that all input files are valid
    if len(args.infiles) == 0:
        raise ValueError("No input files provided.")
    for input_file in args.infiles:
        if not os.path.isfile(input_file):
            raise FileNotFoundError(f"Input file {input_file} not found.")

    # Read in the input TSV files
    print(f"Reading in {len(args.infiles)} input files...")
    merged_df = pd.DataFrame()
    samples = set()
    for input_file in args.infiles:
        df = pd.read_csv(input_file, sep="\t")
        if df.empty:
            print(f"WARNING: Input file {input_file} is empty.")
            continue
        samples.update(df['sample'].unique())
        df = df[(df['ratio'] >= args.min_haplotype_ratio) & (abs(df['diff']) >= args.delta_haplotype_ratio_threshold) & (df['padj'] <= args.padj_threshold)]
        df.insert(0, "input_file", input_file)
        merged_df = pd.concat([merged_df, df], ignore_index=True)

    # Check if the expected column is present
    if not 'gene' in merged_df.columns:
        raise ValueError(f"Input files must contain 'gene' column.")

    # Add a column to indicate the number of samples each gene is an outlier in
    merged_df['sample_count'] = merged_df.groupby('gene')['sample'].transform('nunique')

    # Remove genes that are outliers in >= sample_number_threshold samples
    if args.sample_number_threshold:
        filtered_df = merged_df[merged_df['sample_count'] <= args.sample_number_threshold]
    else:
        filtered_df = merged_df.copy()

    # Write the filtered DataFrame to file
    os.makedirs(os.path.dirname(args.outprefix), exist_ok=True)
    filtered_df.to_csv(args.outprefix + ".tsv", sep="\t", index=False)
    print(f"Merged ASE results and saved to {args.outprefix + '.tsv'}.")

    # Generate a plot of the number of genes before and after filtering
    if args.plot:
        legends = [
            ('#d95d5b', f'Genes with abs(Δ haplotype ratio) ≥ {args.delta_haplotype_ratio_threshold} and padj ≤ {args.padj_threshold}'), \
            ('#4c8fca', f'Genes with allelic imbalance in ≤ {args.sample_number_threshold} samples')
        ]
        plot_outlier_counts(samples, [merged_df, filtered_df], legends, args.title, f"{args.outprefix}_counts.pdf")

if __name__ == "__main__":
    main()