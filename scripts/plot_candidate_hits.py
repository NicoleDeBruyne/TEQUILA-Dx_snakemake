#!/usr/bin/env python3
 
# Author: Nicole DeBruyne (Lin Lab)
# Date: 2025.10.01

import argparse
import os, glob, pandas as pd, ast, matplotlib.pyplot as plt, numpy as np
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator

rcParams['pdf.fonttype'] = 42

def parse_args():
    """ Parse command line arguments """

    parser = argparse.ArgumentParser(description="Plot clustermaps for feature counts")
    parser.add_argument("--infile", required=True, help="Input file containing candidate genes. Expected columns: \
                                    sample, gene, variant, pathogenic_variant, ASE, outlier_junction, and optionally max_AMELIE_score and AMELIE_rank")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--amelie_score_threshold", type=float, help="AMELIE score threshold for plots.")
    parser.add_argument("--amelie_rank_threshold", type=float, help="AMELIE rank threshold for plots.")

    return parser.parse_args()

def plot(df, samples, catcol, categories, colors, outfile, title, legend_title):
    counts = df.groupby(["sample", catcol])["gene"].nunique().unstack(fill_value=0)
    counts = counts.reindex(samples, columns=categories, fill_value=0)

    # ---------- Barplot ----------
    fig, ax = plt.subplots(figsize=(12,6))
    counts.plot(kind="bar", stacked=True, ax=ax, color=colors)
    ax.set_xticklabels(samples, rotation=90)
    ax.set_xlabel("")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel("Gene Count")
    ax.set_title(title)
    ax.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    barplot_file = f"{os.path.splitext(outfile)[0]}_barplot{os.path.splitext(outfile)[1]}"
    plt.savefig(barplot_file)
    plt.close(fig)
    print(f"Bar plot saved to {barplot_file}")

    # ---------- Boxplot ----------
    total_counts = counts.sum(axis=1)
    fig, (ax_box, ax_bar) = plt.subplots(1, 2, figsize=(10,6), sharey=True, gridspec_kw={'width_ratios':[1,1.5]})
    ax_box.boxplot(total_counts, positions=[1], widths=0.4,
                patch_artist=True, boxprops=dict(facecolor='none', color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                flierprops=dict(marker=''))
    x_jitter = np.random.normal(1, 0.04, size=len(total_counts))
    ax_box.scatter(x_jitter, total_counts, color='#d95d5b', alpha=0.7, zorder=3)
    ax_box.set_xticks([1])
    ax_box.set_xticklabels(["Box Plot"])
    ax_box.set_ylabel("Gene Count")
    unique, freq = np.unique(total_counts, return_counts=True)
    ax_bar.barh(unique, freq, color='#d95d5b', alpha=0.7)
    ax_bar.set_xlabel("Number of Samples")
    plt.suptitle(f"Avg={total_counts.mean():.2f}, Med={total_counts.median()}, "
                f"Range={total_counts.min()}–{total_counts.max()}", fontsize=14)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    outfile2 = f"{os.path.splitext(outfile)[0]}_boxplot{os.path.splitext(outfile)[1]}"
    plt.savefig(outfile2)
    plt.close(fig)
    print(f"Box plot saved to {outfile2}")

def main():
    """ Main function """

    # Parse command line arguments
    args = parse_args()

    # Read in input file
    df = pd.read_csv(args.infile, sep='\t', keep_default_na=False)
    df = df.astype(object)
    df.fillna('.', inplace=True)

    # Add some columns
    df['RNA_dysregulation'] = (df['ASE'] | (df['outlier_junction'] != 'None'))
    df['final_grouping'] = df.apply(lambda row: 
        "RNA dysregulation with candidate variant" if row['RNA_dysregulation'] and row['variant'] else
        ("RNA dysregulation" if row['RNA_dysregulation'] else
        ("Candidate variant" if row['variant'] else "No findings")), axis=1
    )

    # Plot
    samples = sorted(df['sample'].unique())

    plot(df[df['pathogenic_variant']].copy(), samples, "pathogenic_variant", [True, False], ["#d95d5b", "#ea9a9c"], f'{args.outdir}/genes_with_pathogenic_variant.pdf', "Candidate Genes with Pathogenic Variant", "Pathogenic Variant")
    plot(df, samples, "ASE", [True], ["#57aa3e"], f'{args.outdir}/genes_with_ASE.pdf', "Candidate Genes", "Allelic Imbalance")
    plot(df, samples, "outlier_junction", ["Strong", "Moderate", "Weak"], ["#4c8fca", "#91c4e9", "#a9cce6"], f'{args.outdir}/genes_with_outlier_junction.pdf', "Candidate Genes", "Alternative Splicing")
    plot(df, samples, "final_grouping", ["RNA dysregulation", "RNA dysregulation with candidate variant"], ["#b271ab", "#cca0ca"], f'{args.outdir}/genes_with_RNA_dysregulation.pdf', "Candidate Genes", "Group")

    if args.amelie_score_threshold:
        filtered_df = df[pd.to_numeric(df['max_AMELIE_score'], errors='coerce') >= args.amelie_score_threshold]
        plot(filtered_df[filtered_df['pathogenic_variant']].copy(), samples, "pathogenic_variant", [True, False], ["#d95d5b", "#ea9a9c"], f'{args.outdir}/genes_with_pathogenic_variant_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes with Pathogenic Variant, AMELIE Score >= {args.amelie_score_threshold}", "Pathogenic Variant")
        plot(filtered_df, samples, "ASE", [True], ["#57aa3e"], f'{args.outdir}/genes_with_ASE_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes, AMELIE Score >= {args.amelie_score_threshold}", "Allelic Imbalance")
        plot(filtered_df, samples, "outlier_junction", ["Strong", "Moderate", "Weak"], ["#4c8fca", "#91c4e9", "#a9cce6"], f'{args.outdir}/genes_with_outlier_junction_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes, AMELIE Score >= {args.amelie_score_threshold}", "Outlier Junction")
        plot(filtered_df, samples, "final_grouping", ["RNA dysregulation", "RNA dysregulation with candidate variant"], ["#b271ab", "#cca0ca"], f'{args.outdir}/genes_with_RNA_dysregulation_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes, AMELIE Score >= {args.amelie_score_threshold}", "Group")

    if args.amelie_rank_threshold:
        filtered_df = df[pd.to_numeric(df['AMELIE_rank'], errors='coerce') <= args.amelie_rank_threshold]
        plot(filtered_df[filtered_df['pathogenic_variant']].copy(), samples, "pathogenic_variant", [True, False], ["#d95d5b", "#ea9a9c"], f'{args.outdir}/genes_with_pathogenic_variant_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes with Pathogenic Variant, AMELIE Rank <= {args.amelie_rank_threshold}", "Pathogenic Variant")
        plot(filtered_df, samples, "ASE", [True], ["#57aa3e"], f'{args.outdir}/genes_with_ASE_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes, AMELIE Rank <= {args.amelie_rank_threshold}", "Allelic Imbalance")
        plot(filtered_df, samples, "outlier_junction", ["Strong", "Moderate", "Weak"], ["#4c8fca", "#91c4e9", "#a9cce6"], f'{args.outdir}/genes_with_outlier_junction_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes, AMELIE Rank <= {args.amelie_rank_threshold}", "Outlier Junction")
        plot(filtered_df, samples, "final_grouping", ["RNA dysregulation", "RNA dysregulation with candidate variant"], ["#b271ab", "#cca0ca"], f'{args.outdir}/genes_with_RNA_dysregulation_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes, AMELIE Rank <= {args.amelie_rank_threshold}", "Group")

if __name__ == "__main__":
    main()
