#!/usr/bin/env python3
"""
scripts/quantify_gene_expression_sample.py
Approximates one sample's relative gene expression across the panel's genes,
using one of two lightweight proxies (no external quantification tool):

  --metric count    -- number of distinct reads overlapping each gene's BED
                        region (a reasonable proxy for transcript abundance
                        with full-length long reads, where each read is
                        roughly one transcript molecule).
  --metric coverage -- max per-base pileup depth anywhere in each gene's BED
                        region (much more sensitive to exactly where reads
                        pile up -- e.g. one probe/amplicon-covered exon --
                        than to overall transcript abundance, but included
                        as an alternative/sanity-check view).

Invoked per-sample by rules/7_sample_gene_quantification.smk (_7A for
count, _7B for coverage).

Writes both the raw value (read count, or max depth) and "counts per target
million" (CPTM): raw_value / (sum of every gene's raw value for THIS
sample) * 1e6 -- i.e. relative to the total signal across every gene *on
this panel* for this one sample, not to the sample's total sequencing
depth. This makes values comparable across samples regardless of depth,
while staying meaningful for a targeted panel (where "fraction of all
on-target reads" would be diluted by off-gene-body panel regions like
flanking/intronic probes). CPTM's denominator is entirely this sample's own
data, so it's exactly as valid computed here, per-sample, as it was when
computed after pooling a whole cohort's raw values into one matrix.
"""

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
    """{gene: (chrom, start, end)} -- one row per gene, column 4 is the gene
    symbol. Matches the BED convention used throughout this pipeline (e.g.
    scripts/phase_reads.py's extract_gene_regions): if a gene appears more
    than once, the last row wins."""
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
    """Primary, mapped alignments only -- matches the read-filtering
    convention already used for on-target counting elsewhere in this
    pipeline (scripts/get_on_target_rate.py), so gene-level counts here are
    consistent with the cohort's on-target-rate numbers."""
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
    depth = np.array(per_base).sum(axis=0)  # sum A/C/G/T arrays -> per-base depth
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

    # CPTM = value / (sum of every gene's value for THIS sample) * 1e6. A
    # sample where every gene is 0 (e.g. a failed/empty BAM) would divide by
    # zero -- leave CPTM as 0 for that sample rather than NaN/inf.
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
