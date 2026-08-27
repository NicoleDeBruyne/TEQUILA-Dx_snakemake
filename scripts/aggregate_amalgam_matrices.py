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

    transcript_tables = []
    gene_tables = []
    for sample, f in zip(args.samples, args.infiles):
        df = pd.read_csv(f, sep='\t')
        tx = df[['transcript_id', 'gene_id', 'count']].rename(columns={'count': sample})
        transcript_tables.append(tx)
        gene = df.groupby('gene_id')['count'].sum().to_frame(name=sample)
        gene_tables.append(gene)
        print(f"Processed {sample}", flush=True)

    transcript_matrix = transcript_tables[0]
    for tx in transcript_tables[1:]:
        transcript_matrix = transcript_matrix.merge(tx, on=['transcript_id', 'gene_id'], how='outer')
    transcript_matrix = transcript_matrix.fillna(0)

    gene_matrix = pd.concat(gene_tables, axis=1).fillna(0)

    out_transcript = args.outprefix + '_transcript_matrix.tsv'
    out_gene = args.outprefix + '_gene_matrix.tsv'
    transcript_matrix.to_csv(out_transcript, sep='\t', index=False)
    gene_matrix.to_csv(out_gene, sep='\t')
    print(f"Done: {out_transcript} and {out_gene} written.", flush=True)


if __name__ == '__main__':
    main()
