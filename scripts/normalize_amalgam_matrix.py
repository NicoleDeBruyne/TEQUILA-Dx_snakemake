
import argparse
import os

import pandas as pd
import numpy as np

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from expression_outliers import compute_outlier_scores
from motr import add_motr_args, compute_size_factors
from gene_boxplots import make_gene_boxplots


def parse_args():
    parser = argparse.ArgumentParser(
        description="CPTM/MOTR normalization, outlier z-scores, and per-gene boxplots for the "
                    "AMALGAM gene expression matrix -- the same outputs the other three "
                    "quantification methods produce. Expects a symbol-keyed, BED-panel-filtered "
                    "raw count matrix (gene_amalgam_matrix_raw.tsv from aggregate_amalgam_matrices.py).")
    parser.add_argument('--raw-matrix', required=True,
        help="by_amalgam/gene_amalgam_matrix_raw.tsv from scripts/aggregate_amalgam_matrices.py "
             "(raw counts, gene-symbol-keyed, BED-panel genes only).")
    parser.add_argument('--outprefix', required=True,
        help="Prefix for output files: <outprefix>_matrix_cptm.tsv, "
             "<outprefix>_matrix_motr.tsv, <outprefix>_zscores_cptm.tsv, "
             "<outprefix>_zscores_motr.tsv, and "
             "<outdir>/<gene>_cptm.pdf + <outdir>/<gene>_motr.pdf per targeted-panel gene "
             "(gene PDFs are written next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    add_motr_args(parser)
    add_alias_map_arg(parser)
    add_outlier_args_local(parser)
    return parser.parse_args()


def add_outlier_args_local(parser):
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


def main():
    args = parse_args()

    # N5 (aggregate_amalgam_matrices.py) has already translated gene_id -> symbol,
    # collapsed duplicates, and filtered to the BED panel.
    raw_df = pd.read_csv(args.raw_matrix, sep="\t", index_col=0)
    raw_df.index.name = "gene"

    alias_map = parse_alias_map(args.alias_map)

    col_sums = raw_df.sum(axis=0)
    cptm_df = raw_df.div(col_sums.replace(0, np.nan), axis=1).fillna(0) * 1e6
    cptm_df.index.name = "gene"

    out_matrix = args.outprefix + "_matrix_cptm.tsv"
    cptm_df.to_csv(out_matrix, sep="\t")
    print("Saved targeted-panel CPTM matrix: " + out_matrix)

    alias_cptm_df = cptm_df.copy()
    alias_cptm_df.columns = resolve_all(alias_cptm_df.columns, alias_map)
    out_matrix_alias = args.outprefix + "_matrix_cptm_alias.tsv"
    alias_cptm_df.to_csv(out_matrix_alias, sep="\t")
    print("Saved alias-labeled targeted-panel CPTM matrix: " + out_matrix_alias)

    size_factors = compute_size_factors(raw_df, max_zero_fraction=args.motr_max_zero_fraction)
    motr_df = raw_df.div(size_factors, axis=1)
    motr_df.index.name = "gene"

    out_matrix_motr = args.outprefix + "_matrix_motr.tsv"
    motr_df.to_csv(out_matrix_motr, sep="\t")
    print("Saved targeted-panel MOTR matrix: " + out_matrix_motr)

    alias_motr_df = motr_df.copy()
    alias_motr_df.columns = resolve_all(alias_motr_df.columns, alias_map)
    out_matrix_motr_alias = args.outprefix + "_matrix_motr_alias.tsv"
    alias_motr_df.to_csv(out_matrix_motr_alias, sep="\t")
    print("Saved alias-labeled targeted-panel MOTR matrix: " + out_matrix_motr_alias)

    plot_outdir = os.path.dirname(args.outprefix)
    make_gene_boxplots(cptm_df, plot_outdir, "CPTM (AMALGAM)", "_cptm")
    make_gene_boxplots(motr_df, plot_outdir, "MOTR (AMALGAM)", "_motr")
    print("Saved per-gene CPTM/MOTR boxplots to: " + plot_outdir)

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
