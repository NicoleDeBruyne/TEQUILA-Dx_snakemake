# 10F/10G - merge gene counts (or coverage) into cohort matrices: raw, CPTM and MOTR

import argparse

import pandas as pd

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from motr import add_motr_args, compute_size_factors


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
             "<outprefix>_matrix_raw.tsv, "
             "and <outdir>/<gene>_cptm.pdf + <outdir>/<gene>_motr.pdf per gene (gene PDFs are written "
             "next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    add_motr_args(parser)
    add_alias_map_arg(parser)
    return parser.parse_args()


def main():
    args = parse_args()

    raw_col = "raw_count" if args.metric == "count" else "raw_coverage"

    # long per-sample tables -> gene x sample
    per_sample_dfs = [pd.read_csv(f, sep="\t") for f in args.infiles]
    samples = [d['sample'].iloc[0] for d in per_sample_dfs if not d.empty]
    long_df = pd.concat(per_sample_dfs, ignore_index=True)
    genes = list(long_df['gene'].drop_duplicates())

    cptm_df = long_df.pivot(index='gene', columns='sample', values='cptm').reindex(index=genes, columns=samples)
    cptm_df.index.name = "gene"

    raw_df = long_df.pivot(index='gene', columns='sample', values=raw_col).reindex(index=genes, columns=samples)
    raw_df.index.name = "gene"

    alias_map = parse_alias_map(args.alias_map)

    out_matrix = args.outprefix + "_matrix_cptm.tsv"
    cptm_df.to_csv(out_matrix, sep="\t")
    print("Saved CPTM matrix: " + out_matrix)

    alias_cptm_df = cptm_df.copy()
    alias_cptm_df.columns = resolve_all(alias_cptm_df.columns, alias_map)
    out_matrix_alias = args.outprefix + "_matrix_cptm_alias.tsv"
    alias_cptm_df.to_csv(out_matrix_alias, sep="\t")
    print("Saved alias-labeled CPTM matrix: " + out_matrix_alias)

    out_raw = args.outprefix + "_matrix_raw.tsv"
    raw_df.to_csv(out_raw, sep="\t")
    print("Saved raw-value matrix: " + out_raw)

    # MOTR has to come from the raw counts, not CPTM
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


if __name__ == "__main__":
    main()
