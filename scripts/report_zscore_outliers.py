
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_BOX_COLOR = "black"
_POINT_COLOR = "#2c7fb8"
_OUTLIER_COLOR = "#c0392b"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize a gene x sample z-score matrix (e.g. *_zscores_cptm.tsv or "
                    "*_zscores_motr.tsv from quantify_gene_by_assignment.py/quantify_gene_expression.py/ "
                    "normalize_amalgam_matrix.py -- any normalization method or quantification source, "
                    "same gene-row/sample-column shape): for a given threshold, report how many genes "
                    "are below threshold in each sample and how many samples each gene is below "
                    "threshold in, plus a boxplot+strip of the raw z-scores (one box per gene, one "
                    "box per sample) with a companion bar chart of the below-threshold counts.")
    parser.add_argument("--zscore-matrix", required=True,
        help="Path to the z-score matrix TSV (first column 'gene', one column per sample).")
    parser.add_argument("--threshold", type=float, default=-3.5,
        help="Outlier cutoff; a value counts if it is <= this threshold (default: -3.5).")
    parser.add_argument("--title", default=None,
        help="Prefix for plot titles (e.g. 'cohort1 panelA blood, MOTR'). Defaults to the input filename.")
    parser.add_argument("--outprefix", required=True,
        help="Writes <outprefix>_by_sample.tsv, <outprefix>_by_gene.tsv, and "
             "<outprefix>_by_sample_boxplot.pdf, <outprefix>_by_gene_boxplot.pdf.")
    return parser.parse_args()


def _plot_boxplot_with_threshold(data, threshold, ylabel, count_label, title, outfile):
    """data: ordered dict/list of (category_name, 1-D array of z-scores). One box+jittered
    strip per category on top, a bar of how many of that category's points are <= threshold
    on the bottom, sharing the same x-axis/category order."""
    categories = [c for c, _ in data]
    n = len(categories)
    fig_w = max(6.0, min(40.0, n * 0.28 + 1.5))
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(fig_w, 8), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )

    box_data = [np.asarray(vals, dtype=float) for _, vals in data]
    _line = dict(color=_BOX_COLOR, linewidth=1)
    ax1.boxplot(
        box_data, showfliers=False, widths=0.6,
        boxprops=_line, whiskerprops=_line, capprops=_line, medianprops=_line,
    )

    rng = np.random.RandomState(0)
    below_counts = []
    for i, vals in enumerate(box_data, start=1):
        vals = vals[~np.isnan(vals)]
        jitter = (rng.rand(len(vals)) - 0.5) * 0.35
        is_below = vals <= threshold
        below_counts.append(int(is_below.sum()))
        if (~is_below).any():
            ax1.scatter(i + jitter[~is_below], vals[~is_below], color=_POINT_COLOR,
                        s=10, alpha=0.7, zorder=3, linewidths=0)
        if is_below.any():
            ax1.scatter(i + jitter[is_below], vals[is_below], color=_OUTLIER_COLOR,
                        s=22, zorder=4, linewidths=0)

    ax1.axhline(threshold, linestyle=":", color=_BOX_COLOR, linewidth=1.2, zorder=2)
    ax1.annotate(f"z = {threshold}", xy=(1, threshold), xycoords=("axes fraction", "data"),
                 xytext=(4, 3), textcoords="offset points", fontsize=8, va="bottom")
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)

    ax2.bar(range(1, n + 1), below_counts, color=_OUTLIER_COLOR, width=0.6)
    ax2.set_ylabel(count_label)
    ax2.set_xlim(0.3, n + 0.7)
    ax2.set_xticks(range(1, n + 1))
    tick_fontsize = 8 if n <= 40 else (6 if n <= 100 else 4)
    rotation = 45 if n <= 15 else 90
    ax2.set_xticklabels(categories, rotation=rotation, ha="right" if rotation == 45 else "center",
                         fontsize=tick_fontsize)

    fig.savefig(outfile)
    plt.close(fig)


def main():
    args = parse_args()
    threshold = args.threshold
    label = args.title or args.zscore_matrix

    df = pd.read_csv(args.zscore_matrix, sep="\t", index_col=0)
    df.index.name = "gene"
    is_outlier = df <= threshold

    by_sample = (
        is_outlier.sum(axis=0)
        .rename("n_genes_below_threshold")
        .rename_axis("sample")
        .reset_index()
        .sort_values("n_genes_below_threshold", ascending=False)
        .reset_index(drop=True)
    )
    by_gene = (
        is_outlier.sum(axis=1)
        .rename("n_samples_below_threshold")
        .rename_axis("gene")
        .reset_index()
        .sort_values("n_samples_below_threshold", ascending=False)
        .reset_index(drop=True)
    )

    out_by_sample = args.outprefix + "_by_sample.tsv"
    out_by_gene = args.outprefix + "_by_gene.tsv"
    by_sample.to_csv(out_by_sample, sep="\t", index=False)
    by_gene.to_csv(out_by_gene, sep="\t", index=False)

    # Per-gene boxplot: one box per gene, points = that gene's z-score in each sample.
    # Genes ordered by median z-score (ascending) so the most-dysregulated genes cluster on the left.
    gene_order = df.median(axis=1, skipna=True).sort_values(ascending=True).index
    gene_data = [(gene, df.loc[gene].dropna().to_numpy()) for gene in gene_order]
    out_gene_boxplot = args.outprefix + "_by_gene_boxplot.pdf"
    _plot_boxplot_with_threshold(
        gene_data, threshold,
        ylabel="z-score (across samples)",
        count_label=f"n samples\nz <= {threshold}",
        title=f"{label}\nper-gene z-score distribution ({len(gene_data)} gene(s))",
        outfile=out_gene_boxplot,
    )

    # Per-sample boxplot: one box per sample, points = that sample's z-score in each gene.
    # Samples ordered by median z-score (ascending) for the same reason.
    sample_order = df.median(axis=0, skipna=True).sort_values(ascending=True).index
    sample_data = [(sample, df[sample].dropna().to_numpy()) for sample in sample_order]
    out_sample_boxplot = args.outprefix + "_by_sample_boxplot.pdf"
    _plot_boxplot_with_threshold(
        sample_data, threshold,
        ylabel="z-score (across genes)",
        count_label=f"n genes\nz <= {threshold}",
        title=f"{label}\nper-sample z-score distribution ({len(sample_data)} sample(s))",
        outfile=out_sample_boxplot,
    )

    n_genes, n_samples = df.shape
    print(f"{n_genes} gene(s) x {n_samples} sample(s), threshold z <= {threshold}")
    print(f"Total outlier (gene, sample) pairs: {int(is_outlier.to_numpy().sum())}")
    print(f"Saved per-sample counts to {out_by_sample}")
    print(f"Saved per-gene counts to {out_by_gene}")
    print(f"Saved per-gene boxplot to {out_gene_boxplot}")
    print(f"Saved per-sample boxplot to {out_sample_boxplot}")


if __name__ == "__main__":
    main()
