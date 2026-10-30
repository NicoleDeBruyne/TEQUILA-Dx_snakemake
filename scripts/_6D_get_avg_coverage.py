# 6D - depth summary per gene (median/mean/SD/breadth) over the canonical exons

import argparse
import os

import numpy as np
import pandas as pd
import pysam
from gtf_utils import ref_blocks, load_canonical_transcripts, load_gene_list


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute one sample's per-gene depth summary over each gene's canonical transcript")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bed", required=True,
        help="BED file for this sample's panel; column 4 is the gene symbol. Only used for the gene list.")
    parser.add_argument("--gtf", required=True, help="Path to the reference annotation GTF (plain or .gz).")
    parser.add_argument("--min-depth", type=int, default=50,
        help="Depth a base must reach to count toward breadth_coverage. Default: 50")
    parser.add_argument("--outfile", required=True)
    return parser.parse_args()


# difference array instead of a real pileup, way faster and good enough for depth
def exonic_depth(bam_handle, chrom, span_start, span_end, exons):
    span_len = span_end - span_start
    diff = np.zeros(span_len + 1, dtype=np.int64)
    for read in bam_handle.fetch(chrom, span_start, span_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        for b_start, b_end in ref_blocks(read):
            lo = max(b_start, span_start) - span_start
            hi = min(b_end, span_end) - span_start
            if lo < hi:
                diff[lo] += 1
                diff[hi] -= 1
    depth = np.cumsum(diff[:-1])
    mask = np.zeros(span_len, dtype=bool)
    for e_start, e_end in exons:
        mask[e_start - span_start:e_end - span_start] = True
    return depth[mask]


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

    def empty_row(gene, length=np.nan):
        return dict(gene=gene, length=length, median_coverage=np.nan,
                    average_coverage=np.nan, std_coverage=np.nan, breadth_coverage=np.nan)

    rows = []
    bam = pysam.AlignmentFile(args.bam, "rb") if os.path.exists(args.bam) else None
    if bam is None:
        print("[WARNING] BAM not found, writing NaN coverage for every gene: " + args.bam)
    try:
        for gene in genes:
            if gene not in transcripts:
                rows.append(empty_row(gene))
                continue
            chrom, span_start, span_end, exons, exonic_length = transcripts[gene]
            if bam is None or exonic_length <= 0:
                rows.append(empty_row(gene, exonic_length))
                continue
            depth = exonic_depth(bam, chrom, span_start, span_end, exons)
            rows.append(dict(
                gene=gene,
                length=int(exonic_length),
                median_coverage=float(np.median(depth)),
                average_coverage=float(np.mean(depth)),
                std_coverage=float(np.std(depth)),
                breadth_coverage=float(np.mean(depth >= args.min_depth)),
            ))
    finally:
        if bam is not None:
            bam.close()

    df = pd.DataFrame(rows, columns=["gene", "length", "median_coverage", "average_coverage", "std_coverage",
                                     "breadth_coverage"])
    df["length"] = df["length"].astype("Int64")
    df.to_csv(args.outfile, sep="\t", index=False, float_format="%.4f")
    print(f"{args.sample}: coverage computed for {len(genes)} gene(s). Saved: {args.outfile}")


if __name__ == "__main__":
    main()
