#!/usr/bin/env python3
"""
scripts/quantify_gene_by_assignment.py
A fourth relative-gene-expression proxy (alongside quantify_gene_expression.py's
--metric count/coverage): assigns each alignment to the single gene, among
every gene annotated in the GTF (not just genes on the run's BED panel), that
it shares the most annotated splice sites with.

For every gene in the GTF, this parses:
  - its genomic span (for a genome-wide overlap index, so a read is only ever
    compared against genes it could plausibly belong to, not all ~60,000+
    genes in the GTF)
  - the union of exon intervals across all its transcripts (used only for
    unspliced reads, see below)
  - the *set* of individual splice-site positions (not junction pairs) across
    every intron of every transcript of that gene: ss1 = the first base of
    the intron, ss2 = the last base of the intron, both 1-based -- the same
    convention used for junction strings elsewhere in this pipeline (see
    scripts/identify_cohort_junction_outliers.py's parse_gtf_junctions).
    Splice sites are pooled from every transcript of a gene and compared as
    individual positions (not donor/acceptor pairs), since two transcripts of
    the same gene can share a donor while differing at the acceptor, and a
    read should be able to "vote" for a gene based on either site
    independently.

Each primary or supplementary alignment (secondary and unmapped alignments
are skipped) in each BAM is then assigned as follows:
  - Spliced alignment (>=1 'N' CIGAR op): its own splice-site positions
    (derived from its CIGAR, same ss1/ss2 convention) are intersected against
    every GTF gene whose span it overlaps. It's assigned to the gene with the
    largest number of shared sites, PROVIDED that gene is unique (i.e. not
    tied with another gene) and shares at least one site. Ties, and reads
    whose sites match no overlapping gene at all, are left unassigned.
  - Unspliced alignment (no 'N' CIGAR op): assigned to a gene only if its
    aligned span overlaps that gene's exons and no *other* gene's exons.
    Overlapping the exons of more than one gene, or no gene's exons at all,
    leaves it unassigned.

Unlike --metric count/coverage in quantify_gene_expression.py (which only
ever see genes on the run's BED panel, since they're computed directly from
BED regions), assignment happens against the full GTF gene set. Two raw
matrices are written: one across every GTF gene that received >=1 assigned
read anywhere in the cohort, and one restricted to the run's BED panel genes
(with an explicit zero row for any panel gene that received none). CPTM
("counts per target million") and the per-gene boxplots are both computed
only from the BED-panel matrix, normalized against the panel's own total
assigned-read signal -- the same convention --metric count/coverage use --
since a genome-wide CPTM denominator (~every GTF gene) would not be
meaningful for a targeted panel.

Invoked by rules/10_quantify_genes.smk (_10D).
"""

import argparse
import gzip
import os
import re
from collections import defaultdict
import concurrent.futures

import pandas as pd
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42


# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Assign each read to the GTF gene it shares the most annotated splice sites with, '
                     'and approximate relative gene expression across a cohort from the resulting counts')
    parser.add_argument(
        "--mapping-file",
        required=True,
        help="TSV file without a header and with the columns: name, bam")
    parser.add_argument(
        "--bed",
        required=True,
        help="BED file for the panel shared by every sample in --mapping-file. "
             "One row per gene (chrom, start, end, gene, ...); column 4 is the gene symbol. "
             "Only used to pick the CPTM/boxplot subset -- reads are assigned against every gene in --gtf.")
    parser.add_argument(
        "--gtf",
        required=True,
        help="GTF/GTF.gz annotation. Every gene in this file (not just BED-panel genes) is a candidate "
             "assignment target.")
    parser.add_argument(
        "--outprefix",
        required=True,
        help="Prefix for output files: <outprefix>_matrix.tsv (targeted-panel CPTM), "
             "<outprefix>_matrix_raw.tsv (targeted-panel raw counts), "
             "<outprefix>_matrix_raw_all_genes.tsv (every GTF gene with >=1 assigned read cohort-wide), "
             "and <outdir>/<gene>.pdf per targeted-panel gene "
             "(gene PDFs are written next to the matrix, not under the prefix's basename)")
    parser.add_argument('--title')
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# BED (targeted panel -- only used to select the CPTM/boxplot subset)
# ---------------------------------------------------------------------------

def load_targeted_genes(bed):
    """Ordered list of gene symbols from BED column 4 -- matches the BED
    convention used throughout this pipeline (e.g. quantify_gene_expression.py's
    load_gene_regions): if a gene appears more than once, the last row wins,
    but its position in the returned order is its *first* occurrence."""
    genes = []
    seen = set()
    with open(bed) as b:
        for line in b:
            if not line.strip() or line.startswith('#'):
                continue
            gene = line.strip().split('\t')[3]
            if gene not in seen:
                genes.append(gene)
                seen.add(gene)
    return genes


# ---------------------------------------------------------------------------
# GTF parsing -- every gene's span, merged exons, and splice-site positions
# ---------------------------------------------------------------------------

_ATTR_RE_CACHE = {}


def _attr(attr_str, key):
    rx = _ATTR_RE_CACHE.get(key)
    if rx is None:
        rx = re.compile(key + r'\s+"([^"]+)"')
        _ATTR_RE_CACHE[key] = rx
    m = rx.search(attr_str)
    return m.group(1) if m else None


def _strip_ver(s):
    return re.sub(r'\.\d+$', '', s) if s else s


def _merge_intervals(intervals):
    """Sorted, merged list of (start, end) half-open intervals."""
    if not intervals:
        return []
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(iv) for iv in merged]


def parse_gtf(gtf_path):
    """Returns {gene: {"chrom": str, "start": int, "end": int,
                        "exons": [(start, end), ...] (merged, sorted),
                        "splice_sites": set(int)}}
    over every gene with >=1 exon in the GTF. start/end are 0-based
    half-open (matching the merged exon intervals); splice site positions
    are 1-based (ss1 = first intron base, ss2 = last intron base), matching
    the junction-string convention used elsewhere in this pipeline."""

    open_fn = gzip.open if gtf_path.endswith(".gz") else open

    # gene -> chrom, and gene -> transcript_id -> [(start0, end), ...]
    gene_chrom = {}
    gene_tx_exons = defaultdict(lambda: defaultdict(list))

    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue

            chrom = parts[0]
            start0 = int(parts[3]) - 1   # 0-based
            end = int(parts[4])          # half-open
            attrs = parts[8]

            gname = _attr(attrs, "gene_name") or _attr(attrs, "gene_symbol")
            gid = _strip_ver(_attr(attrs, "gene_id"))
            gene_key = gname or gid
            if gene_key is None:
                continue

            tx_id = _attr(attrs, "transcript_id")
            if tx_id is None:
                continue

            gene_chrom[gene_key] = chrom
            gene_tx_exons[gene_key][tx_id].append((start0, end))

    genes = {}
    for gene_key, tx_dict in gene_tx_exons.items():
        all_exons = [iv for exons in tx_dict.values() for iv in exons]
        merged_exons = _merge_intervals(all_exons)
        if not merged_exons:
            continue

        splice_sites = set()
        for exons in tx_dict.values():
            if len(exons) < 2:
                continue
            exons_sorted = sorted(exons, key=lambda x: x[0])
            for i in range(len(exons_sorted) - 1):
                _, end_e = exons_sorted[i]
                start_n, _ = exons_sorted[i + 1]
                ss1 = end_e + 1     # first intron base, 1-based
                ss2 = start_n       # last intron base, 1-based
                splice_sites.add(ss1)
                splice_sites.add(ss2)

        genes[gene_key] = {
            "chrom": gene_chrom[gene_key],
            "start": merged_exons[0][0],
            "end": max(e for _, e in merged_exons),
            "exons": merged_exons,
            "splice_sites": splice_sites,
        }

    return genes


# ---------------------------------------------------------------------------
# Genome-wide gene index -- so a read is only ever compared against the
# handful of genes it overlaps, not all ~60,000+ genes in the GTF.
# ---------------------------------------------------------------------------

_BIN_SIZE = 100_000


def build_gene_bins(genes, bin_size=_BIN_SIZE):
    bins = defaultdict(list)
    for gene_key, info in genes.items():
        b_start = info["start"] // bin_size
        b_end = (max(info["end"] - 1, info["start"])) // bin_size
        for b in range(b_start, b_end + 1):
            bins[(info["chrom"], b)].append(gene_key)
    return bins


def candidate_genes(chrom, start, end, genes, gene_bins, bin_size=_BIN_SIZE):
    """Gene keys whose span actually overlaps [start, end) on chrom,
    restricted up front to the bins that span overlaps."""
    b_start = start // bin_size
    b_end = max(end - 1, start) // bin_size
    seen = set()
    out = []
    for b in range(b_start, b_end + 1):
        for gene_key in gene_bins.get((chrom, b), ()):
            if gene_key in seen:
                continue
            seen.add(gene_key)
            info = genes[gene_key]
            if info["start"] < end and info["end"] > start:
                out.append(gene_key)
    return out


def _overlaps_exons(exons, start, end):
    """True if [start, end) overlaps any interval in the sorted, merged
    exons list."""
    for e_start, e_end in exons:
        if e_start >= end:
            break
        if e_end > start:
            return True
    return False


# ---------------------------------------------------------------------------
# Per-read splice-site extraction
# ---------------------------------------------------------------------------

_CIGAR_CONSUMES_REF = {0, 2, 3, 7, 8}  # M, D, N, =, X


def read_splice_sites(read):
    """Set of this alignment's own splice-site positions (1-based ss1/ss2
    per intron, pooled), derived from its CIGAR. Empty set for an unspliced
    (no 'N' op) alignment."""
    sites = set()
    ref_pos = read.reference_start  # 0-based
    for op, length in read.cigartuples:
        if op == 3:  # N -- intron
            sites.add(ref_pos + 1)        # ss1: first intron base, 1-based
            sites.add(ref_pos + length)   # ss2: last intron base, 1-based
            ref_pos += length
        elif op in _CIGAR_CONSUMES_REF:  # M, D, =, X
            ref_pos += length
        # I, S, H, P: consume query only, not reference -- skip
    return sites


def _keep_read(r):
    """Primary and supplementary alignments, excluding unmapped and
    secondary -- per this method's read-filtering convention (each
    alignment record, primary or supplementary, is evaluated and assigned
    independently)."""
    return not r.is_unmapped and not r.is_secondary


# ---------------------------------------------------------------------------
# Per-sample assignment
# ---------------------------------------------------------------------------

def assign_sample(sample, bam, genes, gene_bins):
    """Returns (sample, {gene: assigned_read_count}, stats_dict)."""

    counts = defaultdict(int)
    n_total = n_spliced_assigned = n_spliced_unassigned = 0
    n_unspliced_assigned = n_unspliced_unassigned = 0

    with pysam.AlignmentFile(bam, "rb") as f:
        for read in f.fetch(until_eof=True):
            if not _keep_read(read):
                continue
            if read.cigartuples is None:
                continue
            n_total += 1

            sites = read_splice_sites(read)
            cands = candidate_genes(read.reference_name, read.reference_start,
                                     read.reference_end, genes, gene_bins)

            if sites:
                # Spliced: assign to the unique gene with the most shared
                # splice-site positions, provided it shares >=1 and isn't tied.
                best_gene, best_count, n_at_best = None, 0, 0
                for gene_key in cands:
                    shared = len(sites & genes[gene_key]["splice_sites"])
                    if shared > best_count:
                        best_gene, best_count, n_at_best = gene_key, shared, 1
                    elif shared == best_count and shared > 0:
                        n_at_best += 1
                if best_gene is not None and best_count > 0 and n_at_best == 1:
                    counts[best_gene] += 1
                    n_spliced_assigned += 1
                else:
                    n_spliced_unassigned += 1
            else:
                # Unspliced: assign only if the read's span overlaps exactly
                # one candidate gene's exons.
                overlapping = [
                    gene_key for gene_key in cands
                    if _overlaps_exons(genes[gene_key]["exons"], read.reference_start, read.reference_end)
                ]
                if len(overlapping) == 1:
                    counts[overlapping[0]] += 1
                    n_unspliced_assigned += 1
                else:
                    n_unspliced_unassigned += 1

    stats = {
        "n_total": n_total,
        "n_spliced_assigned": n_spliced_assigned,
        "n_spliced_unassigned": n_spliced_unassigned,
        "n_unspliced_assigned": n_unspliced_assigned,
        "n_unspliced_unassigned": n_unspliced_unassigned,
    }
    return sample, dict(counts), stats


# ---------------------------------------------------------------------------
# Plotting (same convention as quantify_gene_expression.py's per-gene boxplots)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    df = pd.read_csv(args.mapping_file, sep="\t", header=None, names=["sample", "bam"])
    targeted_genes = load_targeted_genes(args.bed)
    if not targeted_genes:
        raise ValueError("No gene regions found in BED file: " + args.bed)

    print("Parsing GTF " + args.gtf + " for every annotated gene's splice sites...")
    genes = parse_gtf(args.gtf)
    if not genes:
        raise ValueError("No genes with exons found in GTF file: " + args.gtf)
    gene_bins = build_gene_bins(genes)
    print("Parsed " + str(len(genes)) + " gene(s) from the GTF (assignment target universe).")

    print("Assigning reads across " + str(len(df)) + " sample(s), " + str(len(targeted_genes)) +
          " targeted gene(s) on the BED panel, using " + str(args.threads) + " thread(s)...")

    raw = {}
    all_stats = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(assign_sample, r["sample"], r["bam"], genes, gene_bins): r["sample"]
                   for _, r in df.iterrows()}
        for f in concurrent.futures.as_completed(futures):
            sample, gene_counts, stats = f.result()
            raw[sample] = gene_counts
            all_stats[sample] = stats
            print("Finished: " + sample + " (" + str(stats["n_spliced_assigned"] + stats["n_unspliced_assigned"]) +
                  "/" + str(stats["n_total"]) + " alignments assigned)")

    # Every GTF gene that received >=1 assigned read anywhere in the cohort.
    raw_all_df = pd.DataFrame(raw).fillna(0).astype(int)
    raw_all_df = raw_all_df[df["sample"].tolist()]
    raw_all_df.index.name = "gene"
    raw_all_df = raw_all_df.sort_index()

    out_raw_all = args.outprefix + "_matrix_raw_all_genes.tsv"
    raw_all_df.to_csv(out_raw_all, sep="\t")
    print("Saved all-gene raw-value matrix: " + out_raw_all)

    # Targeted (BED-panel) subset -- explicit zero row for any panel gene
    # with no assigned reads at all, matching --metric count/coverage's
    # convention of always listing every panel gene.
    raw_df = raw_all_df.reindex(targeted_genes).fillna(0).astype(int)
    raw_df.index.name = "gene"

    out_raw = args.outprefix + "_matrix_raw.tsv"
    raw_df.to_csv(out_raw, sep="\t")
    print("Saved targeted-panel raw-value matrix: " + out_raw)

    # CPTM = value / (sum of every *targeted* gene's value for that sample) * 1e6.
    # Samples where every targeted gene is 0 would divide by zero -- leave
    # those columns as 0 CPTM rather than NaN/inf.
    col_sums = raw_df.sum(axis=0)
    cptm_df = raw_df.div(col_sums.replace(0, np.nan), axis=1).fillna(0) * 1e6

    out_matrix = args.outprefix + "_matrix.tsv"
    cptm_df.to_csv(out_matrix, sep="\t")
    print("Saved targeted-panel CPTM matrix: " + out_matrix)

    metric_label = "assigned reads"
    plot_outdir = os.path.dirname(args.outprefix)
    make_gene_boxplots(cptm_df, plot_outdir, metric_label)
    print("Saved per-gene boxplots to: " + plot_outdir)


if __name__ == "__main__":
    main()
