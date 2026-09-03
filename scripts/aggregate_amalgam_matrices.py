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
"""

import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-sample AMALGAM transcript quantification files into "
                    "cohort-wide transcript_matrix.tsv and gene_matrix.tsv.")
    parser.add_argument('--infiles', required=True, nargs='+',
        help="Every sample's <sample>_transcript_quantification.tsv in this group.")
    parser.add_argument('--samples', required=True, nargs='+',
        help="Sample name for each --infiles entry, same order/length.")
    parser.add_argument('--outprefix', required=True,
        help="Writes <outprefix>_transcript_matrix.tsv and <outprefix>_gene_matrix.tsv.")
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
    # shared reference transcriptome from _10C3's filtered.gtf), so
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


if __name__ == '__main__':
    main()