# 10H - same as 10F but for the splice-site assignment counts, plus a summary of assignment outcomes

import argparse

import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import rcParams

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from motr import add_motr_args, compute_size_factors
rcParams['pdf.fonttype'] = 42


def parse_args():
    parser = argparse.ArgumentParser(
        description='Merge per-sample splice-site-assignment gene TSVs into a cohort matrix + boxplots')
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Every sample's {sample}_gene_assignment.tsv, in the desired column order")
    parser.add_argument("--stats-infiles", nargs="+", required=True,
        help="Every sample's {sample}_read_outcomes.tsv, same order as --infiles")
    parser.add_argument(
        "--outprefix",
        required=True,
        help="Prefix for output files: <outprefix>_matrix_cptm.tsv (targeted-panel CPTM), "
             "<outprefix>_matrix_motr.tsv (targeted-panel MOTR -- median-of-target-ratios, "
             "DESeq2-style size-factor normalization), "
             "<outprefix>_matrix_raw.tsv (targeted-panel raw counts), "
             "<outprefix>_matrix_raw_all_genes.tsv (every gene with >=1 assigned read cohort-wide), "
             "and "
             "<outdir>/<gene>_cptm.pdf + <outdir>/<gene>_motr.pdf per targeted-panel gene "
             "(gene PDFs are written next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    add_motr_args(parser)
    add_alias_map_arg(parser)
    return parser.parse_args()


_STATS_CATEGORIES = [
    ("spliced_assigned",                   "Spliced: assigned"),
    ("spliced_unassigned_zero_shared",     "Spliced: 0 shared splice sites"),
    ("spliced_unassigned_tied",            "Spliced: tied genes"),
    ("unspliced_assigned",                 "Unspliced: assigned"),
    ("unspliced_unassigned_zero_overlap",  "Unspliced: 0 overlapping gene exons"),
    ("unspliced_unassigned_multi_overlap", "Unspliced: multiple overlapping genes"),
    ("unspliced_unassigned_not_contained", "Unspliced: not entirely within exon"),
]


def main():
    args = parse_args()

    per_sample_dfs = [pd.read_csv(f, sep="\t") for f in args.infiles]
    samples = [d['sample'].iloc[0] for d in per_sample_dfs if not d.empty]
    long_df = pd.concat(per_sample_dfs, ignore_index=True)

    all_genes = sorted(long_df['gene'].drop_duplicates())
    raw_all_df = (
        long_df.pivot(index='gene', columns='sample', values='raw_count')
        .reindex(index=all_genes, columns=samples).fillna(0).astype(int)
    )
    raw_all_df.index.name = "gene"

    out_raw_all = args.outprefix + "_matrix_raw_all_genes.tsv"
    raw_all_df.to_csv(out_raw_all, sep="\t")
    print("Saved all-gene raw-value matrix: " + out_raw_all)

    # only panel genes have a CPTM value, so use that to restrict
    targeted_genes = sorted(long_df.loc[long_df['cptm'].notna(), 'gene'].drop_duplicates())
    raw_df = raw_all_df.reindex(targeted_genes).fillna(0).astype(int)
    raw_df.index.name = "gene"

    out_raw = args.outprefix + "_matrix_raw.tsv"
    raw_df.to_csv(out_raw, sep="\t")
    print("Saved targeted-panel raw-value matrix: " + out_raw)

    cptm_df = (
        long_df.pivot(index='gene', columns='sample', values='cptm')
        .reindex(index=targeted_genes, columns=samples)
    )
    cptm_df.index.name = "gene"

    out_matrix = args.outprefix + "_matrix_cptm.tsv"
    cptm_df.to_csv(out_matrix, sep="\t")
    print("Saved targeted-panel CPTM matrix: " + out_matrix)

    alias_map = parse_alias_map(args.alias_map)
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


    # how reads got assigned in each sample (useful for spotting weird samples)
    stats_df = pd.concat([pd.read_csv(f, sep="\t") for f in args.stats_infiles], ignore_index=True)
    stats_df = stats_df.set_index("sample").reindex(samples)
    stats_df = stats_df[["n_total"] + [col for col, _ in _STATS_CATEGORIES]]

    out_stats = args.outprefix + "_read_outcomes.tsv"
    stats_df.to_csv(out_stats, sep="\t")
    print("Saved per-sample assignment-outcome stats: " + out_stats)


if __name__ == "__main__":
    main()
