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

Also computes a second normalization of the same raw counts: MOTR
("median of target ratios") -- DESeq2's "poscounts" median-of-ratios
size-factor normalization -- see scripts/motr.py's module docstring for
the algorithm and its per-sample exclusion rule. Written to
<outprefix>_matrix_motr.tsv (+ alias copy); the existing CPTM matrix is
<outprefix>_matrix_cptm.tsv.

Two per-gene boxplots are written per gene (<gene>_cptm.pdf,
<gene>_motr.pdf) -- see scripts/gene_boxplots.py.

Also writes a low-expression outlier score for every (gene, sample): a
leave-one-out, shrinkage-based robust z-score in log2 space -- see
scripts/expression_outliers.py's module docstring for the full algorithm
and why it's used instead of a parametric (e.g. negative-binomial) fit.
"""

import argparse
import os

import pandas as pd
import numpy as np

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from expression_outliers import compute_outlier_scores
from motr import add_motr_args, compute_size_factors
from gene_boxplots import make_gene_boxplots


def add_outlier_args(parser):
    """Shared CLI options for the low-expression outlier score -- same
    defaults/semantics across all four quantification methods. See
    expression_outliers.py's module docstring for the algorithm."""
    parser.add_argument("--outlier-pseudocount", type=float, default=1.0,
        help="Added to CPTM before log2-transforming. Default: 1.0")
    parser.add_argument("--outlier-shrinkage-k", type=float, default=10.0,
        help="Shrinkage constant: weight on a gene's own leave-one-out MAD is "
             "n/(n+k) vs. the cohort-wide trend. Default: 10.0")
    parser.add_argument("--outlier-min-mad", type=float, default=0.1,
        help="Floor on the shrunk MAD (log2 units), guarding against a near-zero-MAD "
             "gene producing an absurd z-score. Default: 0.1")
    parser.add_argument("--outlier-zscore-threshold", type=float, default=3.0,
        help="A sample is flagged as a low-expression outlier for a gene when "
             "z <= -this value. Default: 3.0")


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
        help="Prefix for output files: <outprefix>_matrix_cptm.tsv, <outprefix>_matrix_motr.tsv, "
             "<outprefix>_matrix_raw.tsv, <outprefix>_zscores_cptm.tsv, <outprefix>_zscores_motr.tsv, "
             "and <outdir>/<gene>_cptm.pdf + <outdir>/<gene>_motr.pdf per gene (gene PDFs are written "
             "next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    add_motr_args(parser)
    add_alias_map_arg(parser)
    add_outlier_args(parser)
    return parser.parse_args()


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

    alias_map = parse_alias_map(args.alias_map)

    # Primary matrix: clean genes x samples layout of the relative-expression
    # (CPTM) value.
    out_matrix = args.outprefix + "_matrix_cptm.tsv"
    cptm_df.to_csv(out_matrix, sep="\t")
    print("Saved CPTM matrix: " + out_matrix)

    # Alias-labeled copy: always produced (mirrors the real-ID matrix
    # verbatim when --alias-map is empty), so the rule's declared output
    # exists regardless of whether this cohort actually has any aliases.
    alias_cptm_df = cptm_df.copy()
    alias_cptm_df.columns = resolve_all(alias_cptm_df.columns, alias_map)
    out_matrix_alias = args.outprefix + "_matrix_cptm_alias.tsv"
    alias_cptm_df.to_csv(out_matrix_alias, sep="\t")
    print("Saved alias-labeled CPTM matrix: " + out_matrix_alias)

    # Secondary matrix: the same layout with raw (un-normalized) values, for
    # reference/debugging -- e.g. distinguishing a true zero-expression gene
    # from a panel-design dropout, which CPTM alone can't tell apart. Every
    # gene here is already BED-panel-restricted (done per-sample, upstream,
    # by scripts/quantify_gene_expression_sample.py) -- there's no separate
    # genome-wide matrix for this method.
    out_raw = args.outprefix + "_matrix_raw.tsv"
    raw_df.to_csv(out_raw, sep="\t")
    print("Saved raw-value matrix: " + out_raw)

    # MOTR ("median of target ratios"): DESeq2-style median-of-ratios
    # normalization -- see scripts/motr.py's module docstring.
    size_factors = compute_size_factors(raw_df, max_zero_fraction=args.motr_max_zero_fraction)
    motr_df = raw_df.div(size_factors, axis=1)
    motr_df.index.name = "gene"

    out_matrix_motr = args.outprefix + "_matrix_motr.tsv"
    motr_df.to_csv(out_matrix_motr, sep="\t")
    print("Saved MOTR matrix: " + out_matrix_motr)

    alias_motr_df = motr_df.copy()
    alias_motr_df.columns = resolve_all(alias_motr_df.columns, alias_map)
    out_matrix_motr_alias = args.outprefix + "_matrix_motr_alias.tsv"
    alias_motr_df.to_csv(out_matrix_motr_alias, sep="\t")
    print("Saved alias-labeled MOTR matrix: " + out_matrix_motr_alias)

    plot_outdir = os.path.dirname(args.outprefix)
    make_gene_boxplots(cptm_df, plot_outdir, "CPTM (" + metric_label + ")", "_cptm")
    make_gene_boxplots(motr_df, plot_outdir, "MOTR (" + metric_label + ")", "_motr")
    print("Saved per-gene CPTM/MOTR boxplots to: " + plot_outdir)

    # Low-expression outlier score (see expression_outliers.py's module
    # docstring for the algorithm), computed once per normalization -- a
    # gene/sample flagged as low-expression on CPTM may not be on MOTR
    # (or vice versa), since MOTR's per-sample size factor can shift a
    # gene's relative rank within its own sample, so both are kept rather
    # than picking one.
    z_cptm_df, _ = compute_outlier_scores(
        cptm_df,
        pseudocount=args.outlier_pseudocount,
        shrinkage_k=args.outlier_shrinkage_k,
        min_mad=args.outlier_min_mad,
        z_threshold=args.outlier_zscore_threshold,
    )
    out_zscores_cptm = args.outprefix + "_zscores_cptm.tsv"
    z_cptm_df.to_csv(out_zscores_cptm, sep="\t")
    print("Saved CPTM z-score matrix: " + out_zscores_cptm)

    z_motr_df, _ = compute_outlier_scores(
        motr_df,
        pseudocount=args.outlier_pseudocount,
        shrinkage_k=args.outlier_shrinkage_k,
        min_mad=args.outlier_min_mad,
        z_threshold=args.outlier_zscore_threshold,
    )
    out_zscores_motr = args.outprefix + "_zscores_motr.tsv"
    z_motr_df.to_csv(out_zscores_motr, sep="\t")
    print("Saved MOTR z-score matrix: " + out_zscores_motr)


if __name__ == "__main__":
    main()
