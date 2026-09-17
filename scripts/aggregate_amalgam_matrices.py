#!/usr/bin/env python3
"""
scripts/aggregate_amalgam_matrices.py

Combines every sample's AMALGAM Quantify_Transcripts.py output
(<sample>_transcript_quantification.tsv, columns: transcript_id, gene_id,
count, ...) in a (bed_id, sample_type) group into two cohort-wide
matrices: one per-transcript, one per-gene (gene-level = sum of that
gene's transcripts' counts). Extracted from the aggregation step of the
group's original manual sbatch pipeline into its own script, matching
this repo's convention of a dedicated scripts/*.py file per pipeline
step rather than inline Python in a rules/*.smk shell block.

The gene matrix's raw counts are also normalized into a CPTM ("counts per
target million") value here -- unlike --metric count/coverage/assignment
(scripts/quantify_gene_expression_sample.py,
scripts/quantify_gene_by_assignment_sample.py), which each compute their
own CPTM per-sample before this cohort-level merge even happens, AMALGAM's
per-sample output is raw transcript counts with no per-sample normalization
step of its own -- so it has to happen here, on the assembled cohort
matrix, the one time this script has every sample's raw values in hand at
once. Same normalization convention as the other three methods: raw_value /
(sum of every gene's raw value for that sample) * 1e6, so all four methods'
CPTM values are on a comparable scale.

Also writes a low-expression outlier score for every (gene, sample), on
this CPTM-normalized gene matrix -- see scripts/expression_outliers.py's
module docstring for the algorithm.
"""

import argparse

import pandas as pd
import numpy as np

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from expression_outliers import compute_outlier_scores, outliers_long_format


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
        description="Combine per-sample AMALGAM transcript quantification files into "
                    "cohort-wide transcript_matrix.tsv and gene_matrix.tsv.")
    parser.add_argument('--infiles', required=True, nargs='+',
        help="Every sample's <sample>_transcript_quantification.tsv in this group.")
    parser.add_argument('--samples', required=True, nargs='+',
        help="Sample name for each --infiles entry, same order/length.")
    parser.add_argument('--outprefix', required=True,
        help="Writes <outprefix>_transcript_matrix.tsv, <outprefix>_gene_matrix.tsv (raw), "
             "<outprefix>_gene_matrix_cptm.tsv (normalized), and outlier score files -- see module docstring.")
    add_alias_map_arg(parser)
    add_outlier_args(parser)
    args = parser.parse_args()
    if len(args.infiles) != len(args.samples):
        parser.error("--infiles and --samples must have the same number of entries")
    return args


def main():
    args = parse_args()

    # Two-pass design. Pass 1 collects the union of every (transcript_id,
    # gene_id) pair across all samples, reading only those two (cheap,
    # low-cardinality-per-file) columns, folding each file's ids into a
    # single running index one at a time and discarding the file's own
    # frame immediately. Pass 2 then reads each sample's counts and
    # aligns them to that ONE shared index.
    #
    # This matters because earlier versions of this fix (see git history)
    # held every sample's own id/count data in memory simultaneously via
    # a growing Python list before doing anything with it -- first with
    # full (transcript_id, gene_id, count) tables in a single-pass design,
    # then, even after splitting into two passes, by collecting pass 1's
    # per-file (transcript_id, gene_id) frames into a list before
    # deduplicating them (same bug, just 2 columns instead of 3 -- caught
    # because the job died with zero "Processed <sample>" lines printed,
    # meaning it never even reached pass 2). transcript_id/gene_id are
    # strings duplicated identically across nearly every sample (same
    # shared reference transcriptome from _9D3's filtered.gtf), so
    # holding N samples' worth of them at once is what actually exceeded
    # the job's memory limit -- not the final matrix-assembly step.
    #
    # Building one shared tx_index up front, incrementally, means every
    # sample's data resident in memory afterward is just a numeric count
    # array aligned to that single shared index, rather than its own full
    # copy of every transcript/gene name string.
    tx_index = None
    for f in args.infiles:
        ids = pd.read_csv(f, sep='\t', usecols=['transcript_id', 'gene_id'])
        idx = pd.MultiIndex.from_frame(ids)
        tx_index = idx.unique() if tx_index is None else tx_index.union(idx, sort=False)

    transcript_matrix = pd.DataFrame(index=tx_index)
    gene_tables = []
    for sample, f in zip(args.samples, args.infiles):
        df = pd.read_csv(f, sep='\t', usecols=['transcript_id', 'gene_id', 'count'])
        transcript_matrix[sample] = (
            df.set_index(['transcript_id', 'gene_id'])['count'].reindex(tx_index)
        )
        gene = df.groupby('gene_id')['count'].sum().to_frame(name=sample)
        gene_tables.append(gene)
        print(f"Processed {sample}", flush=True)

    transcript_matrix = transcript_matrix.fillna(0).reset_index()

    gene_matrix = pd.concat(gene_tables, axis=1).fillna(0)

    out_transcript = args.outprefix + '_transcript_matrix.tsv'
    out_gene = args.outprefix + '_gene_matrix.tsv'
    transcript_matrix.to_csv(out_transcript, sep='\t', index=False)
    gene_matrix.to_csv(out_gene, sep='\t')
    print(f"Done: {out_transcript} and {out_gene} written.", flush=True)

    # Alias-labeled copy of the gene matrix: always produced (mirrors the
    # real-ID matrix verbatim when --alias-map is empty), so the rule's
    # declared output exists regardless of whether this group actually has
    # any aliases.
    alias_map = parse_alias_map(args.alias_map)
    alias_gene_matrix = gene_matrix.copy()
    alias_gene_matrix.columns = resolve_all(alias_gene_matrix.columns, alias_map)
    out_gene_alias = args.outprefix + '_gene_matrix_alias.tsv'
    alias_gene_matrix.to_csv(out_gene_alias, sep='\t')
    print(f"Done: {out_gene_alias} written.", flush=True)

    # CPTM normalization (raw_value / that sample's own gene-sum * 1e6) --
    # see this module's docstring for why it happens here rather than
    # per-sample, unlike --metric count/coverage/assignment.
    col_sums = gene_matrix.sum(axis=0)
    cptm_df = gene_matrix.div(col_sums.replace(0, np.nan), axis=1).fillna(0) * 1e6
    out_cptm = args.outprefix + '_gene_matrix_cptm.tsv'
    cptm_df.to_csv(out_cptm, sep='\t')
    print(f"Done: {out_cptm} written.", flush=True)

    alias_cptm_df = cptm_df.copy()
    alias_cptm_df.columns = resolve_all(alias_cptm_df.columns, alias_map)
    out_cptm_alias = args.outprefix + '_gene_matrix_cptm_alias.tsv'
    alias_cptm_df.to_csv(out_cptm_alias, sep='\t')
    print(f"Done: {out_cptm_alias} written.", flush=True)

    # Low-expression outlier score (see expression_outliers.py's module
    # docstring for the algorithm), on the CPTM-normalized gene matrix.
    z_df, is_outlier_df = compute_outlier_scores(
        cptm_df,
        pseudocount=args.outlier_pseudocount,
        shrinkage_k=args.outlier_shrinkage_k,
        min_mad=args.outlier_min_mad,
        z_threshold=args.outlier_zscore_threshold,
    )
    out_zscores = args.outprefix + '_outlier_zscores.tsv'
    z_df.to_csv(out_zscores, sep='\t')
    print(f"Done: {out_zscores} written.", flush=True)

    out_outliers = args.outprefix + '_outliers.tsv'
    outliers_long_format(z_df, is_outlier_df).to_csv(out_outliers, sep='\t', index=False)
    print(f"Done: {out_outliers} written.", flush=True)


if __name__ == '__main__':
    main()