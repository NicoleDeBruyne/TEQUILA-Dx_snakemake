#!/usr/bin/env python3
"""
scripts/quantify_gene_expression.py
Approximates relative gene expression across a cohort of BAMs sharing one
BED panel, using one of two lightweight proxies (no external quantification
tool):

  --metric count    -- number of distinct reads overlapping each gene's BED
                        region (a reasonable proxy for transcript abundance
                        with full-length long reads, where each read is
                        roughly one transcript molecule).
  --metric coverage -- max per-base pileup depth anywhere in each gene's BED
                        region (much more sensitive to exactly where reads
                        pile up -- e.g. one probe/amplicon-covered exon --
                        than to overall transcript abundance, but included
                        as an alternative/sanity-check view).

Invoked by rules/10_quantify_genes.smk (_10A for count, _10B for coverage).

Both metrics are reported two ways in the output matrix:
  - the raw value (read count, or max depth)
  - "counts per target million" (CPTM): raw_value / (sum of every gene's
    raw value for that sample) * 1e6 -- i.e. relative to the total signal
    across every gene *on this panel* for that sample, not to the sample's
    total sequencing depth. This makes values comparable across samples
    regardless of depth, while staying meaningful for a targeted panel
    (where "fraction of all on-target reads" would be diluted by
    off-gene-body panel regions like flanking/intronic probes).

The per-gene boxplot (one page per gene) always plots the CPTM value, with
the highest- and lowest-CPTM sample labeled directly on the plot.
"""

import argparse
import os
import concurrent.futures

import pandas as pd
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42


def parse_args():
    parser = argparse.ArgumentParser(description='Approximate relative gene expression across a cohort from BAM read counts or coverage')
    parser.add_argument(
        "--mapping-file",
        required=True,
        help="TSV file without a header and with the columns: name, bam")
    parser.add_argument(
        "--bed",
        required=True,
        help="BED file for the panel shared by every sample in --mapping-file. "
             "One row per gene (chrom, start, end, gene, ...); column 4 is the gene symbol.")
    parser.add_argument(
        "--metric",
        required=True,
        choices=["count", "coverage"],
        help="count = number of distinct reads overlapping the gene region; "
             "coverage = max per-base pileup depth in the gene region")
    parser.add_argument(
        "--outprefix",
        required=True,
        help="Prefix for output files: <outprefix>_matrix.tsv, and <outdir>/<gene>.pdf per gene "
             "(gene PDFs are written next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads")
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
    pipeline (scripts/plot_on_target_rates.py), so gene-level counts here
    are consistent with the cohort's on-target-rate numbers."""
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


def quantify_sample(sample, bam, gene_regions, metric):
    fn = gene_read_count if metric == "count" else gene_max_coverage
    return sample, {gene: fn(bam, chrom, start, end) for gene, (chrom, start, end) in gene_regions.items()}


def make_gene_boxplots(cptm_df, outdir, metric_label):
    """One page per gene: a boxplot of every sample's CPTM value for that
    gene, with individual sample points overlaid and the highest- and
    lowest-CPTM sample labeled by name."""
    os.makedirs(outdir, exist_ok=True)
    for gene in cptm_df.index:
        vals = cptm_df.loc[gene].dropna()
        out_pdf = os.path.join(outdir, gene + ".pdf")
        if vals.empty:
            continue

        fig, ax = plt.subplots(figsize=(3.6, 4.2))
        _black = dict(color="black")
        ax.boxplot([vals.to_numpy()], showfliers=False, widths=0.5,
                   boxprops=_black, whiskerprops=_black, capprops=_black, medianprops=_black)

        jitter = (np.random.RandomState(0).rand(len(vals)) - 0.5) * 0.15
        ax.scatter(1 + jitter, vals.to_numpy(), color="#2c7fb8", zorder=3, s=18)

        max_sample = vals.idxmax()
        min_sample = vals.idxmin()
        ax.scatter([1 + jitter[list(vals.index).index(max_sample)]], [vals[max_sample]],
                   color="#c0392b", zorder=4, s=30)
        ax.annotate(max_sample, (1 + jitter[list(vals.index).index(max_sample)], vals[max_sample]),
                    textcoords="offset points", xytext=(6, 0), fontsize=7, color="#c0392b", va="center")
        if min_sample != max_sample:
            ax.scatter([1 + jitter[list(vals.index).index(min_sample)]], [vals[min_sample]],
                       color="#c0392b", zorder=4, s=30)
            ax.annotate(min_sample, (1 + jitter[list(vals.index).index(min_sample)], vals[min_sample]),
                        textcoords="offset points", xytext=(6, 0), fontsize=7, color="#c0392b", va="center")

        ax.set_xticks([])
        ax.set_ylabel("CPTM (" + metric_label + ")")
        ax.set_title(gene, fontsize=10)
        fig.tight_layout()
        fig.savefig(out_pdf)
        plt.close(fig)


def main():
    args = parse_args()

    df = pd.read_csv(args.mapping_file, sep="\t", header=None, names=["sample", "bam"])
    gene_regions = load_gene_regions(args.bed)
    if not gene_regions:
        raise ValueError("No gene regions found in BED file: " + args.bed)

    print("Quantifying " + str(len(gene_regions)) + " gene(s) across " + str(len(df)) +
          " sample(s) using metric=" + args.metric + " with " + str(args.threads) + " thread(s)...")

    raw = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(quantify_sample, r["sample"], r["bam"], gene_regions, args.metric): r["sample"]
                   for _, r in df.iterrows()}
        for f in concurrent.futures.as_completed(futures):
            sample, gene_vals = f.result()
            raw[sample] = gene_vals
            print("Finished: " + sample)

    # genes as rows, samples as columns
    raw_df = pd.DataFrame(raw)[df["sample"].tolist()]
    raw_df.index.name = "gene"

    # CPTM = value / (sum of every gene's value for that sample) * 1e6.
    # Samples where every gene is 0 (e.g. a failed/empty BAM) would divide
    # by zero -- leave those columns as 0 CPTM rather than NaN/inf.
    col_sums = raw_df.sum(axis=0)
    cptm_df = raw_df.div(col_sums.replace(0, np.nan), axis=1).fillna(0) * 1e6

    metric_label = "read count" if args.metric == "count" else "max coverage"

    # Primary matrix: clean genes x samples layout of the relative-expression
    # (CPTM) value -- this is the file named exactly by the calling rule
    # (e.g. gene_count_matrix.tsv / gene_coverage_matrix.tsv).
    out_matrix = args.outprefix + "_matrix.tsv"
    cptm_df.to_csv(out_matrix, sep="\t")
    print("Saved CPTM matrix: " + out_matrix)

    # Secondary matrix: the same layout with raw (un-normalized) values, for
    # reference/debugging -- e.g. distinguishing a true zero-expression gene
    # from a panel-design dropout, which CPTM alone can't tell apart.
    out_raw = args.outprefix + "_matrix_raw.tsv"
    raw_df.to_csv(out_raw, sep="\t")
    print("Saved raw-value matrix: " + out_raw)

    plot_outdir = os.path.dirname(args.outprefix)
    make_gene_boxplots(cptm_df, plot_outdir, metric_label)
    print("Saved per-gene boxplots to: " + plot_outdir)


if __name__ == "__main__":
    main()
