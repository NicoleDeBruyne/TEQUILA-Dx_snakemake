#!/usr/bin/env python3
"""
scripts/quantify_gene_by_assignment.py
Cohort-level merge step: combines every sample's own
{sample}_gene_assignment.tsv (per-gene raw_count + CPTM, written per-sample
by scripts/quantify_gene_by_assignment_sample.py via
rules/7_sample_gene_quantification.smk's _7C) and {sample}_read_outcomes.tsv
into the cohort's genome-wide raw matrix, BED-panel raw + CPTM matrices,
per-gene boxplots, and the read-assignment-outcome summary plot. No BAM or
GTF access here -- adding/removing a sample from a cohort only reruns this
cheap merge, not the per-sample BAM scan + GTF parse. Invoked by
rules/9_merge_results.smk's _9M.

CPTM is read straight from the per-sample TSVs, not recomputed here (same
reasoning as scripts/quantify_gene_expression.py). Which genes are
"BED-panel genes" for the CPTM/boxplot subset is inferred from which rows
have a non-null cptm value (populated only for panel genes by the
per-sample step) -- so this merge step doesn't need its own BED file.

Also writes a low-expression outlier score for every (gene, sample), on
the targeted-panel CPTM matrix -- see scripts/expression_outliers.py's
module docstring for the algorithm. This is the outlier score
scripts/merge_hits.py annotates onto merged_all_hits.tsv (the
splice-site-assignment method was chosen for that annotation since,
unlike --metric count/coverage, it isn't restricted to simple BED-region
overlap/pileup and so isn't as easily confounded by a neighboring gene's
reads spilling into a tightly-packed panel region).

Also computes a second normalization of the same targeted-panel raw
counts: MOTR ("median of target ratios") -- DESeq2's "poscounts"
median-of-ratios size-factor normalization (their own fix for the classic
method's all-or-nothing zero-count handling), restricted to genes on this
run's BED panel rather than the whole transcriptome, since that's the
only gene set this pipeline ever reports on. See
compute_size_factors()'s docstring for the algorithm and its per-sample
fallback. Written to a separate <outprefix>_matrix_motr.tsv (+ _alias
copy) -- the existing CPTM matrix is renamed to <outprefix>_matrix_cptm.tsv
(from plain _matrix.tsv) so the two normalizations are unambiguous side
by side.
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
from expression_outliers import compute_outlier_scores
from motr import add_motr_args, compute_size_factors
from gene_boxplots import make_gene_boxplots
rcParams['pdf.fonttype'] = 42


# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

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
             "<outprefix>_zscores_cptm.tsv + <outprefix>_zscores_motr.tsv (low-expression outlier "
             "z-scores, one per normalization), and "
             "<outdir>/<gene>_cptm.pdf + <outdir>/<gene>_motr.pdf per targeted-panel gene "
             "(gene PDFs are written next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    add_motr_args(parser)
    add_alias_map_arg(parser)
    add_outlier_args(parser)
    return parser.parse_args()


# Category display order/labels/colors for the assignment-outcome plot below.
# Greens for assigned, reds/oranges for the specific unassigned reasons --
# grouped spliced-then-unspliced so the two read types are visually adjacent.
_STATS_CATEGORIES = [
    ("spliced_assigned",                   "Spliced: assigned",                        "#2ca25f"),
    ("spliced_unassigned_zero_shared",     "Spliced: 0 shared splice sites",           "#fc9272"),
    ("spliced_unassigned_tied",            "Spliced: tied genes",                      "#de2d26"),
    ("unspliced_assigned",                 "Unspliced: assigned",                      "#66c2a4"),
    ("unspliced_unassigned_zero_overlap",  "Unspliced: 0 overlapping gene exons",      "#fdae6b"),
    ("unspliced_unassigned_multi_overlap", "Unspliced: multiple overlapping genes",    "#e6550d"),
    ("unspliced_unassigned_not_contained", "Unspliced: not entirely within exon",      "#a63603"),
]


def make_assignment_summary_plot(stats_df, out_pdf, title):
    """One stacked horizontal bar per sample: how many of that sample's
    alignments landed in each assigned/unassigned-reason category (see
    _STATS_CATEGORIES). Read counts, not gene counts -- this is about
    assignment outcome, not the per-gene quantification matrices."""
    samples = list(stats_df.index)
    n_samples = len(samples)

    fig_height = max(3, 0.35 * n_samples + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    y_pos = np.arange(n_samples)
    left = np.zeros(n_samples)
    for col, label, color in _STATS_CATEGORIES:
        vals = stats_df[col].to_numpy(dtype=float)
        ax.barh(y_pos, vals, left=left, height=0.7, color=color, label=label)
        left += vals

    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples)
    ax.invert_yaxis()  # first sample at the top
    ax.set_xlabel("Alignments")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    per_sample_dfs = [pd.read_csv(f, sep="\t") for f in args.infiles]
    samples = [d['sample'].iloc[0] for d in per_sample_dfs if not d.empty]
    long_df = pd.concat(per_sample_dfs, ignore_index=True)

    # Every gene that received >=1 assigned read anywhere in the cohort
    # (union across samples' per-sample files -- each already lists only
    # its own nonzero genes plus its own panel genes).
    all_genes = sorted(long_df['gene'].drop_duplicates())
    raw_all_df = (
        long_df.pivot(index='gene', columns='sample', values='raw_count')
        .reindex(index=all_genes, columns=samples).fillna(0).astype(int)
    )
    raw_all_df.index.name = "gene"

    out_raw_all = args.outprefix + "_matrix_raw_all_genes.tsv"
    raw_all_df.to_csv(out_raw_all, sep="\t")
    print("Saved all-gene raw-value matrix: " + out_raw_all)

    # BED-panel subset: rows with a non-null cptm anywhere are panel genes
    # (populated only for panel genes by the per-sample step) -- explicit
    # zero row for any panel gene with no assigned reads at all, matching
    # --metric count/coverage's convention of always listing every panel gene.
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

    # Alias-labeled copy: always produced (mirrors the real-ID matrix
    # verbatim when --alias-map is empty), so the rule's declared output
    # exists regardless of whether this group actually has any aliases.
    alias_map = parse_alias_map(args.alias_map)
    alias_cptm_df = cptm_df.copy()
    alias_cptm_df.columns = resolve_all(alias_cptm_df.columns, alias_map)
    out_matrix_alias = args.outprefix + "_matrix_cptm_alias.tsv"
    alias_cptm_df.to_csv(out_matrix_alias, sep="\t")
    print("Saved alias-labeled targeted-panel CPTM matrix: " + out_matrix_alias)

    # MOTR ("median of target ratios"): DESeq2-style median-of-ratios
    # normalization, restricted to targeted-panel genes only (raw_df, not
    # raw_all_df) -- see compute_size_factors()'s docstring.
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

    metric_label = "assigned reads"
    plot_outdir = os.path.dirname(args.outprefix)
    make_gene_boxplots(cptm_df, plot_outdir, "CPTM (" + metric_label + ")", "_cptm")
    make_gene_boxplots(motr_df, plot_outdir, "MOTR (" + metric_label + ")", "_motr")
    print("Saved per-gene CPTM/MOTR boxplots to: " + plot_outdir)

    # Low-expression outlier score (see expression_outliers.py's module
    # docstring for the algorithm), computed once per normalization -- the
    # CPTM version is the file rules/9_merge_results.smk's _9D2 passes to
    # merge_hits.py to annotate merged_all_hits.tsv (unchanged from
    # before); the MOTR version is new, for anyone comparing the two
    # normalizations' outlier calls directly.
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

    # Per-sample breakdown of assigned vs. unassigned-and-why (see
    # _STATS_CATEGORIES), one row per sample, columns in the same order the
    # plot stacks them in -- ordered by the sample order --infiles was given in.
    stats_df = pd.concat([pd.read_csv(f, sep="\t") for f in args.stats_infiles], ignore_index=True)
    stats_df = stats_df.set_index("sample").reindex(samples)
    stats_df = stats_df[["n_total"] + [col for col, _, _ in _STATS_CATEGORIES]]

    out_stats = args.outprefix + "_read_outcomes.tsv"
    stats_df.to_csv(out_stats, sep="\t")
    print("Saved per-sample assignment-outcome stats: " + out_stats)

    out_stats_pdf = args.outprefix + "_read_outcomes.pdf"
    make_assignment_summary_plot(stats_df, out_stats_pdf, args.title + " - read assignment outcomes")
    print("Saved read-assignment outcome plot: " + out_stats_pdf)


if __name__ == "__main__":
    main()
