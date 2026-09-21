
import argparse

import pandas as pd
import numpy as np
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description="Approximate one sample's relative gene expression from its own BAM")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument(
        "--bed",
        required=True,
        help="BED file for this sample's panel. One row per gene (chrom, start, end, gene, ...); "
             "column 4 is the gene symbol.")
    parser.add_argument(
        "--metric",
        required=True,
        choices=["count", "coverage"],
        help="count = number of distinct reads overlapping the gene region; "
             "coverage = max per-base pileup depth in the gene region")
    parser.add_argument("--outfile", required=True)
    return parser.parse_args()


def load_gene_regions(bed):
    gene_regions = {}
    with open(bed) as b:
        for line in b:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            gene = fields[3]
            gene_regions[gene] = (fields[0], int(fields[1]), int(fields[2]))
    return gene_regions


def _keep_read(r):
    return not r.is_unmapped and not r.is_secondary


def gene_read_count(bam, chrom, start, end):
    ids = set()
    with pysam.AlignmentFile(bam, "rb") as f:
        for r in f.fetch(chrom, start, end):
            if _keep_read(r):
                ids.add(r.query_name)
    return len(ids)


def gene_max_coverage(bam, chrom, start, end):
    with pysam.AlignmentFile(bam, "rb") as f:
        per_base = f.count_coverage(chrom, start, end, quality_threshold=0, read_callback=_keep_read)
    if end <= start:
        return 0
    depth = np.array(per_base).sum(axis=0)
    return int(depth.max()) if depth.size else 0


def main():
    args = parse_args()

    gene_regions = load_gene_regions(args.bed)
    if not gene_regions:
        raise ValueError("No gene regions found in BED file: " + args.bed)

    print(str(len(gene_regions)) + " gene(s), metric=" + args.metric + ", sample=" + args.sample)

    fn = gene_read_count if args.metric == "count" else gene_max_coverage
    raw = {gene: fn(args.bam, chrom, start, end) for gene, (chrom, start, end) in gene_regions.items()}

    raw_col = "raw_count" if args.metric == "count" else "raw_coverage"

    total = sum(raw.values())
    rows = [
        dict(sample=args.sample, gene=gene, **{raw_col: v}, cptm=(v / total * 1e6 if total else 0.0))
        for gene, v in raw.items()
    ]

    df = pd.DataFrame(rows, columns=["sample", "gene", raw_col, "cptm"])
    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"Saved: {args.outfile}")


if __name__ == "__main__":
    main()
