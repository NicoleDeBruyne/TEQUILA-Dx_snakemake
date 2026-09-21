
import argparse
import gzip
import os
import re

import pandas as pd
import numpy as np

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from expression_outliers import compute_outlier_scores
from motr import add_motr_args, compute_size_factors
from gene_boxplots import make_gene_boxplots


_ATTR_RE_CACHE = {}


def _attr(attr_str, key):
    rx = _ATTR_RE_CACHE.get(key)
    if rx is None:
        rx = re.compile(key + r'\s+"([^"]+)"')
        _ATTR_RE_CACHE[key] = rx
    m = rx.search(attr_str)
    return m.group(1) if m else None


def _strip_ver(s):
    return re.sub(r'\.\d+$', '', s) if s else s


def load_gene_id_to_symbol(gtf_path):
    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    mapping = {}
    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            gid = _strip_ver(_attr(attrs, "gene_id"))
            gname = _attr(attrs, "gene_name") or _attr(attrs, "gene_symbol")
            if gid and gname:
                mapping[gid] = gname
    return mapping


def load_targeted_genes(bed_path):
    genes = []
    seen = set()
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            gene = cols[3]
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return sorted(genes)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Restrict AMALGAM's genome-wide gene_matrix.tsv to the run's BED panel "
                    "(translating gene_id -> gene symbol via the reference GTF) and add "
                    "CPTM/MOTR normalization, z-scores, and per-gene boxplots -- the same "
                    "outputs the other three quantification methods produce.")
    parser.add_argument('--gene-matrix', required=True,
        help="by_amalgam/quantification/*_gene_matrix.tsv from scripts/aggregate_amalgam_matrices.py "
             "(raw counts, gene_id-keyed, genome-wide).")
    parser.add_argument('--bed', required=True,
        help="This group's BED panel (column 4 = target gene symbols).")
    parser.add_argument('--gtf', required=True,
        help="Reference annotation GTF, for gene_id -> gene_name translation.")
    parser.add_argument('--outprefix', required=True,
        help="Prefix for output files: <outprefix>_matrix_raw.tsv, <outprefix>_matrix_cptm.tsv, "
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

    raw_all_df = pd.read_csv(args.gene_matrix, sep="\t", index_col=0)

    gid_to_symbol = load_gene_id_to_symbol(args.gtf)
    symbol_index = [gid_to_symbol.get(_strip_ver(gid), gid) for gid in raw_all_df.index]

    raw_all_df.index = symbol_index
    raw_all_df.index.name = "gene"
    raw_by_symbol = raw_all_df.groupby(level=0).sum()

    targeted_genes = load_targeted_genes(args.bed)
    raw_df = raw_by_symbol.reindex(targeted_genes).fillna(0).astype(int)
    raw_df.index.name = "gene"

    out_raw = args.outprefix + "_matrix_raw.tsv"
    raw_df.to_csv(out_raw, sep="\t")
    print("Saved targeted-panel raw-value matrix: " + out_raw)

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
