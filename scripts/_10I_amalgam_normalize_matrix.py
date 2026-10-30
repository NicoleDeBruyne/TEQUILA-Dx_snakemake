# 10I - CPTM and MOTR normalized versions of the AMALGAM gene matrix

import argparse

import pandas as pd
import numpy as np

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from motr import add_motr_args, compute_size_factors


def parse_args():
    parser = argparse.ArgumentParser(
        description="CPTM/MOTR normalization and per-gene boxplots for the "
                    "AMALGAM gene expression matrix -- the same outputs the other three "
                    "quantification methods produce. Expects a symbol-keyed, BED-panel-filtered "
                    "raw count matrix (gene_amalgam_matrix_raw.tsv from _9E_amalgam_aggregate_matrices.py).")
    parser.add_argument('--raw-matrix', required=True,
        help="by_amalgam/gene_amalgam_matrix_raw.tsv from _9E_amalgam_aggregate_matrices.py "
             "(raw counts, gene-symbol-keyed, BED-panel genes only).")
    parser.add_argument('--outprefix', required=True,
        help="Prefix for output files: <outprefix>_matrix_cptm.tsv, "
             "<outprefix>_matrix_motr.tsv, and "
             "<outdir>/<gene>_cptm.pdf + <outdir>/<gene>_motr.pdf per targeted-panel gene "
             "(gene PDFs are written next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    add_motr_args(parser)
    add_alias_map_arg(parser)
    return parser.parse_args()


def main():
    args = parse_args()

    raw_df = pd.read_csv(args.raw_matrix, sep="\t", index_col=0)
    raw_df.index.name = "gene"

    alias_map = parse_alias_map(args.alias_map)

    # CPTM = counts per ten million panel reads
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


if __name__ == "__main__":
    main()
