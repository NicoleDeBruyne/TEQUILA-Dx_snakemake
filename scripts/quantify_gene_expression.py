#!/usr/bin/env python3
"""
scripts/quantify_gene_expression.py
Cohort-level merge step: combines every sample's own
{sample}_gene_count.tsv or {sample}_gene_coverage.tsv (per-gene raw value +
CPTM, written per-sample by scripts/quantify_gene_expression_sample.py via
rules/7_sample_gene_quantification.smk's _7A/_7B) into the cohort's CPTM +
raw-value matrices and per-gene boxplots. No BAM access here -- adding/
removing a sample from a cohort only reruns this cheap merge, not the
per-sample BAM scan. Invoked by rules/9_merge_results.smk (_9K for count,
_9L for coverage).

CPTM is read straight from the per-sample TSVs, not recomputed here: its
normalization (raw_value / that sample's own gene-sum * 1e6) is entirely a
per-sample computation, so pooling the cohort doesn't change any sample's
value -- see quantify_gene_expression_sample.py's module docstring.

The per-gene boxplot (one page per gene) always plots the CPTM value, with
the highest- and lowest-CPTM sample labeled directly on the plot.
"""

import argparse
import os

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
rcParams['pdf.fonttype'] = 42


def parse_args():
    parser = argparse.ArgumentParser(description='Merge per-sample gene count/coverage TSVs into a cohort matrix + boxplots')
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Every sample's {sample}_gene_count.tsv or {sample}_gene_coverage.tsv, in the desired column order")
    parser.add_argument(
        "--metric",
        required=True,
        choices=["count", "coverage"],
        help="Must match the --metric the per-sample TSVs were generated with -- picks the raw-value "
             "column name (raw_count vs raw_coverage) and the boxplot axis label.")
    parser.add_argument(
        "--outprefix",
        required=True,
        help="Prefix for output files: <outprefix>_matrix.tsv, <outprefix>_matrix_raw.tsv, and "
             "<outdir>/<gene>.pdf per gene (gene PDFs are written next to the matrix, not under the "
             "prefix's basename)")
    parser.add_argument('--title')
    add_alias_map_arg(parser)
    return parser.parse_args()


def make_gene_boxplots(cptm_df, outdir, metric_label):
    """One page per gene: a boxplot of every sample's CPTM value for that
    gene, with individual sample points overlaid and the highest- and
    lowest-CPTM sample labeled by name."""
    os.makedirs(outdir, exist_ok=True)
    for gene in cptm_df.index:
        vals = cptm_df.loc[gene].dropna()
        out_pdf = os.path.join(outdir, gene + ".pdf")
        if vals.empty:
            continue

        fig, ax = plt.subplots(figsize=(3.6, 4.2))
        _black = dict(color="black")
        ax.boxplot([vals.to_numpy()], showfliers=False, widths=0.5,
                   boxprops=_black, whiskerprops=_black, capprops=_black, medianprops=_black)

        jitter = (np.random.RandomState(0).rand(len(vals)) - 0.5) * 0.15
        ax.scatter(1 + jitter, vals.to_numpy(), color="#2c7fb8", zorder=3, s=18)

        max_sample = vals.idxmax()
        min_sample = vals.idxmin()
        ax.scatter([1 + jitter[list(vals.index).index(max_sample)]], [vals[max_sample]],
                   color="#c0392b", zorder=4, s=30)
        ax.annotate(max_sample, (1 + jitter[list(vals.index).index(max_sample)], vals[max_sample]),
                    textcoords="offset points", xytext=(6, 0), fontsize=7, color="#c0392b", va="center")
        if min_sample != max_sample:
            ax.scatter([1 + jitter[list(vals.index).index(min_sample)]], [vals[min_sample]],
                       color="#c0392b", zorder=4, s=30)
            ax.annotate(min_sample, (1 + jitter[list(vals.index).index(min_sample)], vals[min_sample]),
                        textcoords="offset points", xytext=(6, 0), fontsize=7, color="#c0392b", va="center")

        ax.set_xticks([])
        ax.set_ylabel("CPTM (" + metric_label + ")")
        ax.set_title(gene, fontsize=10)
        fig.tight_layout()
        fig.savefig(out_pdf)
        plt.close(fig)


def main():
    args = parse_args()

    raw_col = "raw_count" if args.metric == "count" else "raw_coverage"
    metric_label = "read count" if args.metric == "count" else "max coverage"

    per_sample_dfs = [pd.read_csv(f, sep="\t") for f in args.infiles]
    samples = [d['sample'].iloc[0] for d in per_sample_dfs if not d.empty]
    long_df = pd.concat(per_sample_dfs, ignore_index=True)
    genes = list(long_df['gene'].drop_duplicates())

    cptm_df = long_df.pivot(index='gene', columns='sample', values='cptm').reindex(index=genes, columns=samples)
    cptm_df.index.name = "gene"

    raw_df = long_df.pivot(index='gene', columns='sample', values=raw_col).reindex(index=genes, columns=samples)
    raw_df.index.name = "gene"

    # Primary matrix: clean genes x samples layout of the relative-expression
    # (CPTM) value -- this is the file named exactly by the calling rule
    # (e.g. gene_count_matrix.tsv / gene_coverage_matrix.tsv).
    out_matrix = args.outprefix + "_matrix.tsv"
    cptm_df.to_csv(out_matrix, sep="\t")
    print("Saved CPTM matrix: " + out_matrix)

    # Alias-labeled copy: always produced (mirrors the real-ID matrix
    # verbatim when --alias-map is empty), so the rule's declared output
    # exists regardless of whether this cohort actually has any aliases.
    alias_map = parse_alias_map(args.alias_map)
    alias_cptm_df = cptm_df.copy()
    alias_cptm_df.columns = resolve_all(alias_cptm_df.columns, alias_map)
    out_matrix_alias = args.outprefix + "_matrix_alias.tsv"
    alias_cptm_df.to_csv(out_matrix_alias, sep="\t")
    print("Saved alias-labeled CPTM matrix: " + out_matrix_alias)

    # Secondary matrix: the same layout with raw (un-normalized) values, for
    # reference/debugging -- e.g. distinguishing a true zero-expression gene
    # from a panel-design dropout, which CPTM alone can't tell apart.
    out_raw = args.outprefix + "_matrix_raw.tsv"
    raw_df.to_csv(out_raw, sep="\t")
    print("Saved raw-value matrix: " + out_raw)

    plot_outdir = os.path.dirname(args.outprefix)
    make_gene_boxplots(cptm_df, plot_outdir, metric_label)
    print("Saved per-gene boxplots to: " + plot_outdir)


if __name__ == "__main__":
    main()
