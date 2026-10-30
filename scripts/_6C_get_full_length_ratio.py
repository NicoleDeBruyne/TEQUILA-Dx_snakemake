# 6C - full-length ratio per gene: on average, how much of the canonical transcript's exons each read covers

import argparse
import os

import pandas as pd
import numpy as np
import pysam
from gtf_utils import ref_blocks, load_canonical_transcripts, load_gene_list


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute one sample's per-gene full-length transcript coverage ratio (FLR) from its own BAM")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bed", required=True,
        help="BED file for this sample's panel. One row per gene (chrom, start, end, gene, ...); "
             "column 4 is the gene symbol. Only used to get the gene list -- the actual region "
             "fetched per gene is its canonical transcript's own span from --gtf.")
    parser.add_argument("--gtf", required=True, help="Path to the reference annotation GTF (plain or .gz).")
    parser.add_argument("--outfile", required=True)
    return parser.parse_args()


def _overlap_length(blocks_a, blocks_b):
    i = j = 0
    total = 0
    while i < len(blocks_a) and j < len(blocks_b):
        a_s, a_e = blocks_a[i]
        b_s, b_e = blocks_b[j]
        lo, hi = max(a_s, b_s), min(a_e, b_e)
        if lo < hi:
            total += hi - lo
        if a_e < b_e:
            i += 1
        else:
            j += 1
    return total


# only primary alignments, so every read counts once
def compute_avg_flr(bam_handle, chrom, span_start, span_end, exons, exonic_length):
    if exonic_length <= 0:
        return np.nan, 0
    ratios = []
    for read in bam_handle.fetch(chrom, span_start, span_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        read_blocks = ref_blocks(read)
        if not read_blocks:
            continue
        overlap = _overlap_length(read_blocks, exons)
        ratios.append(overlap / exonic_length)
    if not ratios:
        return np.nan, 0
    return float(np.mean(ratios)), len(ratios)


def main():
    args = parse_args()

    genes = load_gene_list(args.bed)
    if not genes:
        raise ValueError("No genes found in BED file: " + args.bed)

    print("Loading canonical transcripts for " + str(len(genes)) + " gene(s) from " + args.gtf + "...")
    transcripts = load_canonical_transcripts(args.gtf, genes)
    missing = [g for g in genes if g not in transcripts]
    if missing:
        print("[WARNING] No canonical transcript found in GTF for " + str(len(missing)) +
              " gene(s) (will be NaN in output): " + ", ".join(missing[:20]) +
              (" ..." if len(missing) > 20 else ""))

    rows = []
    if os.path.exists(args.bam):
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            for gene in genes:
                if gene not in transcripts:
                    rows.append(dict(sample=args.sample, gene=gene, avgFLR=np.nan, read_count=0))
                    continue
                chrom, span_start, span_end, exons, exonic_length = transcripts[gene]
                avg_flr, n_reads = compute_avg_flr(bam, chrom, span_start, span_end, exons, exonic_length)
                rows.append(dict(sample=args.sample, gene=gene, avgFLR=avg_flr, read_count=n_reads))
    else:
        for gene in genes:
            rows.append(dict(sample=args.sample, gene=gene, avgFLR=np.nan, read_count=0))

    df = pd.DataFrame(rows, columns=["sample", "gene", "avgFLR", "read_count"])
    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"{args.sample}: FLR computed for {len(genes)} gene(s). Saved: {args.outfile}")


if __name__ == "__main__":
    main()
