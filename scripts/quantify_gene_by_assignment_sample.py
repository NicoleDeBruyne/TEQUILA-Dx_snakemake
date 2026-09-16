#!/usr/bin/env python3
"""
scripts/quantify_gene_by_assignment_sample.py
A fourth relative-gene-expression proxy (alongside quantify_gene_expression_sample.py's
--metric count/coverage): assigns each of this sample's alignments to the
single gene, among every gene annotated in the GTF (not just genes on the
run's BED panel), that it shares the most annotated splice sites with.

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
are skipped) in this sample's BAM is then assigned as follows:
  - Spliced alignment (>=1 'N' CIGAR op): its own splice-site positions
    (derived from its CIGAR, same ss1/ss2 convention) are intersected against
    every GTF gene whose span it overlaps. It's assigned to the gene with the
    largest number of shared sites, PROVIDED that gene is unique (i.e. not
    tied with another gene) and shares at least one site. Ties, and reads
    whose sites match no overlapping gene at all, are left unassigned.
  - Unspliced alignment (no 'N' CIGAR op): assigned to a gene only if its
    aligned span overlaps that gene's (merged) exons and no *other* gene's
    exons, AND is entirely CONTAINED within that one gene's merged exon
    set -- e.g. a transcript with an exon 100-200 and another with an exon
    190-250 merge into one 100-250 block, so a read spanning 120-230
    counts as contained even though no single annotated exon covers that
    whole range. If the gene is monoexonic (its merged exon set is a
    single interval), containment is relaxed to plain overlap, since a
    single-exon gene has no internal intron for a partial overlap to fall
    into. Overlapping the exons of more than one gene, overlapping none,
    or overlapping exactly one gene's exons but not being contained within
    them, all leave the read unassigned.

Unlike --metric count/coverage in quantify_gene_expression_sample.py (which
only ever see genes on the run's BED panel, since they're computed directly
from BED regions), assignment happens against the full GTF gene set. The
output TSV includes a raw_count row for every gene that got >=1 assigned
read in this sample, UNION every BED-panel gene (with an explicit 0 row for
a panel gene that got none) -- writing every one of the ~60,000+ GTF genes
regardless would bloat this file for no benefit. CPTM ("counts per target
million") is populated only for BED-panel genes, normalized against the
panel's own total assigned-read signal for this sample -- the same
convention --metric count/coverage use -- since a genome-wide CPTM
denominator (~every GTF gene) would not be meaningful for a targeted panel.
The cohort-level merge step (scripts/quantify_gene_by_assignment.py) infers
which genes are BED-panel genes from which rows have a non-null cptm,
rather than needing its own BED file.

Invoked per-sample by rules/7_sample_gene_quantification.smk (_7C).

Note: unlike --metric count/coverage, this script re-parses the GTF into
its genome-wide splice-site/gene-bin index on every invocation -- one per
sample now, rather than once per cohort submission as before the
per-sample/merge split. GTF parsing is normally far cheaper than the BAM
scan itself, but for a very large annotation this is a real (if usually
small) added cost per sample, traded for the ability to cache each
sample's own assignment result independently of cohort membership.
"""

import argparse
import gzip
import re
from collections import defaultdict

import pandas as pd
import pysam


# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Assign one sample\'s reads to the GTF gene each shares the most annotated splice '
                     'sites with, and compute that sample\'s raw counts + CPTM for its BED-panel genes')
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument(
        "--bed",
        required=True,
        help="BED file for this sample's panel. One row per gene (chrom, start, end, gene, ...); "
             "column 4 is the gene symbol. Only used to pick the CPTM subset -- reads are assigned "
             "against every gene in --gtf.")
    parser.add_argument(
        "--gtf",
        required=True,
        help="GTF/GTF.gz annotation. Every gene in this file (not just BED-panel genes) is a candidate "
             "assignment target.")
    parser.add_argument("--outfile", required=True, help="Per-gene raw_count/cptm TSV.")
    parser.add_argument("--stats-outfile", required=True, help="One-row assignment-outcome stats TSV.")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# BED (targeted panel -- only used to select the CPTM subset)
# ---------------------------------------------------------------------------

def load_targeted_genes(bed):
    """Ordered list of gene symbols from BED column 4 -- matches the BED
    convention used throughout this pipeline (e.g.
    quantify_gene_expression_sample.py's load_gene_regions): if a gene
    appears more than once, the last row wins, but its position in the
    returned order is its *first* occurrence."""
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


def _contained_in_exons(exons, start, end):
    """True if [start, end) is entirely contained within a SINGLE interval
    in the sorted, merged exons list -- e.g. a transcript with an exon
    100-200 and another with an exon 190-250 merge into one 100-250
    interval, so a read spanning 120-230 counts as contained even though
    no single annotated exon covers that whole range. A read that starts
    inside one merged interval and extends past its end (into a genuine
    gap between merged exon blocks) is NOT contained, even if it's fully
    covered by exon sequence from some other transcript not merged into
    this same block."""
    for e_start, e_end in exons:
        if e_start <= start < e_end:
            return end <= e_end
        if e_start > start:
            break
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

def assign_sample(bam, genes, gene_bins):
    """Returns ({gene: assigned_read_count}, stats_dict). stats_dict breaks
    unassigned reads down by the specific reason they were left unassigned
    (see module docstring for the assignment rules each of these
    corresponds to), not just a spliced/unspliced assigned/unassigned
    total."""

    counts = defaultdict(int)
    stats = {
        "n_total": 0,
        "spliced_assigned": 0,
        "spliced_unassigned_zero_shared": 0,      # >=1 candidate gene, but none shares a splice site
        "spliced_unassigned_tied": 0,              # >=2 genes tied for the most shared splice sites
        "unspliced_assigned": 0,
        "unspliced_unassigned_zero_overlap": 0,    # no candidate gene's exons overlapped
        "unspliced_unassigned_multi_overlap": 0,   # >=2 candidate genes' exons overlapped
        "unspliced_unassigned_not_contained": 0,   # exactly 1 gene overlapped, but read not contained in it
    }

    with pysam.AlignmentFile(bam, "rb") as f:
        for read in f.fetch(until_eof=True):
            if not _keep_read(read):
                continue
            if read.cigartuples is None:
                continue
            stats["n_total"] += 1

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
                    stats["spliced_assigned"] += 1
                elif best_count == 0:
                    stats["spliced_unassigned_zero_shared"] += 1
                else:
                    stats["spliced_unassigned_tied"] += 1
            else:
                # Unspliced: assign only if the read's span overlaps
                # exactly one candidate gene's (merged) exons -- the
                # uniqueness gate -- AND is entirely CONTAINED within that
                # one gene's merged exon set. A monoexonic gene (its merged
                # exon set collapses to a single interval) relaxes the
                # second requirement to plain overlap, since a single-exon
                # gene has no internal intron to fall outside of the way a
                # partial-exon overlap on a multiexonic gene would.
                overlapping = [
                    gene_key for gene_key in cands
                    if _overlaps_exons(genes[gene_key]["exons"], read.reference_start, read.reference_end)
                ]
                if len(overlapping) == 0:
                    stats["unspliced_unassigned_zero_overlap"] += 1
                elif len(overlapping) > 1:
                    stats["unspliced_unassigned_multi_overlap"] += 1
                else:
                    gene_key = overlapping[0]
                    gene_exons = genes[gene_key]["exons"]
                    is_monoexonic = len(gene_exons) == 1
                    if is_monoexonic or _contained_in_exons(gene_exons, read.reference_start, read.reference_end):
                        counts[gene_key] += 1
                        stats["unspliced_assigned"] += 1
                    else:
                        stats["unspliced_unassigned_not_contained"] += 1

    return dict(counts), stats


# Column order for the per-sample stats TSV -- same categories the
# cohort-level plot (scripts/quantify_gene_by_assignment.py) stacks, in the
# same order.
_STATS_COLUMNS = [
    "n_total",
    "spliced_assigned", "spliced_unassigned_zero_shared", "spliced_unassigned_tied",
    "unspliced_assigned", "unspliced_unassigned_zero_overlap",
    "unspliced_unassigned_multi_overlap", "unspliced_unassigned_not_contained",
]


def main():
    args = parse_args()

    targeted_genes = load_targeted_genes(args.bed)
    if not targeted_genes:
        raise ValueError("No gene regions found in BED file: " + args.bed)

    print("Parsing GTF " + args.gtf + " for every annotated gene's splice sites...")
    genes = parse_gtf(args.gtf)
    if not genes:
        raise ValueError("No genes with exons found in GTF file: " + args.gtf)
    gene_bins = build_gene_bins(genes)
    print("Parsed " + str(len(genes)) + " gene(s) from the GTF (assignment target universe).")

    print("Assigning reads for " + args.sample + ", " + str(len(targeted_genes)) + " targeted gene(s) on the BED panel...")
    counts, stats = assign_sample(args.bam, genes, gene_bins)
    print(args.sample + ": " + str(stats["spliced_assigned"] + stats["unspliced_assigned"]) +
          "/" + str(stats["n_total"]) + " alignments assigned")

    # Row set: every gene that got >=1 assigned read UNION every BED-panel
    # gene (so a panel gene with 0 assigned reads still gets an explicit
    # zero row, matching --metric count/coverage's convention).
    panel_total = sum(counts.get(g, 0) for g in targeted_genes)
    all_genes = sorted(set(counts) | set(targeted_genes))
    rows = []
    for gene in all_genes:
        raw_count = counts.get(gene, 0)
        is_panel_gene = gene in targeted_genes
        cptm = (raw_count / panel_total * 1e6 if panel_total else 0.0) if is_panel_gene else None
        rows.append(dict(sample=args.sample, gene=gene, raw_count=raw_count, cptm=cptm))

    df = pd.DataFrame(rows, columns=["sample", "gene", "raw_count", "cptm"])
    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"Saved: {args.outfile}")

    stats_row = {"sample": args.sample, **{c: stats[c] for c in _STATS_COLUMNS}}
    stats_df = pd.DataFrame([stats_row], columns=["sample"] + _STATS_COLUMNS)
    stats_df.to_csv(args.stats_outfile, sep="\t", index=False)
    print(f"Saved: {args.stats_outfile}")


if __name__ == "__main__":
    main()
