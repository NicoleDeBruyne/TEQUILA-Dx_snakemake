#!/usr/bin/env python3
"""
scripts/get_full_length_ratio_sample.py
For each gene on the panel, computes one sample's "full-length ratio" (FLR):
the fraction of that gene's canonical transcript's exonic reference
positions that each overlapping read's alignment actually spans via a
non-N (non-intron-skip) CIGAR operation (M/D/=/X), averaged across every
primary, mapped read overlapping the transcript's span. A read that fully
spans every annotated exon of the canonical transcript (no dropped exons,
no truncation) scores close to 1; a read covering only part of the
transcript body (a truncated cDNA/library artifact, a partially-covered
amplicon, an alternate/partial isoform, etc.) scores lower.

Reads directly from this sample's own BAM (the same file _6A/_6B/
quantify_gene_count_sample.py use) with one open file handle, fetching each
gene's canonical-transcript span with pysam's region-indexed fetch() --
NOT phase_reads.py's per-gene bulk BAM files. Those per-gene files are
themselves nothing more than `samtools view <region> <this same original
BAM>` (see phase_reads.py's filter_bam_by_region()), so reading them here
would mean opening one small file per gene -- tens of thousands of file
opens for a large panel -- for data that's one indexed fetch() away in a
file that's opened once anyway. This also means this rule doesn't need to
wait on phase_reads.py at all.

The "canonical transcript" per gene is picked from --gtf: the transcript
tagged "Ensembl_canonical" (GENCODE/Ensembl convention) if one exists for
that gene, otherwise the transcript with the largest total exonic length.
Its genomic span (min exon start to max exon end) is the fetch() region.

Every gene from --bed is written, even ones with no canonical transcript
found in --gtf (NaN avgFLR / 0 read_count) -- the cohort-level merge step
(scripts/get_full_length_ratio.py) assumes every sample's TSV covers the
same gene set for the shared bed panel, so this keeps that true even when
a gene is missing from the GTF for one reason or another.

Invoked per-sample by rules/6_sample_qc.smk (_6C).
"""

import argparse
import gzip
import os

import pandas as pd
import numpy as np
import pysam


# Reference-consuming, non-N CIGAR ops: M=0, D=2, ==7, X=8. I=1/S=4/H=5/P=6
# don't consume reference; N=3 (intron skip) is explicitly excluded.
_REF_CONSUMING_NON_N = {0, 2, 7, 8}


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


def _open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _gtf_attr(attr_str, key):
    # GTF attribute fields look like: gene_name "FOO"; transcript_id "ENST...";
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) != 2:
            continue
        k, v = parts
        if k == key:
            return v.strip().strip('"')
    return None


def load_gene_list(bed):
    """Gene symbols from BED column 4, in file order, de-duplicated -- same
    convention as scripts/quantify_gene_count_sample.py's load_gene_regions()
    and scripts/phase_reads.py's extract_gene_regions()."""
    genes = []
    seen = set()
    with open(bed) as b:
        for line in b:
            if not line.strip() or line.startswith('#'):
                continue
            gene = line.strip().split('\t')[3]
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return genes


def load_canonical_transcripts(gtf_path, genes):
    """Returns {gene: (chrom, span_start, span_end, [(exon_start, exon_end), ...], exonic_length)}
    for the canonical transcript of every gene in `genes` found in the GTF.
    All coordinates are 0-based half-open, matching pysam's reference
    coordinate convention. (chrom, span_start, span_end) is the transcript's
    overall genomic footprint (min exon start to max exon end) -- the region
    passed to AlignmentFile.fetch()."""
    genes = set(genes)

    # transcript_id -> dict(gene, chrom, strand, canonical, exons=[(start,end),...])
    transcripts = {}

    with _open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields[:9]
            if feature not in ("transcript", "exon"):
                continue
            gene_name = _gtf_attr(attrs, "gene_name") or _gtf_attr(attrs, "gene_id")
            if gene_name not in genes:
                continue
            transcript_id = _gtf_attr(attrs, "transcript_id")
            if transcript_id is None:
                continue

            if feature == "transcript":
                is_canonical = 'tag "Ensembl_canonical"' in attrs
                t = transcripts.setdefault(transcript_id, {
                    "gene": gene_name, "chrom": chrom, "strand": strand,
                    "canonical": False, "exons": [],
                })
                t["canonical"] = t["canonical"] or is_canonical
            else:  # exon
                t = transcripts.setdefault(transcript_id, {
                    "gene": gene_name, "chrom": chrom, "strand": strand,
                    "canonical": False, "exons": [],
                })
                # GTF is 1-based inclusive -> convert to 0-based half-open.
                t["exons"].append((int(start) - 1, int(end)))

    # Pick, per gene: the Ensembl_canonical-tagged transcript if any, else
    # the transcript with the largest total exonic length.
    by_gene = {}
    for tid, t in transcripts.items():
        if not t["exons"]:
            continue
        by_gene.setdefault(t["gene"], []).append(t)

    result = {}
    for gene, cands in by_gene.items():
        canonical_cands = [c for c in cands if c["canonical"]]
        pool = canonical_cands if canonical_cands else cands
        best = max(pool, key=lambda c: sum(e - s for s, e in c["exons"]))
        exons = sorted(best["exons"])
        # Merge any overlapping/adjacent exon records (defensive -- GENCODE
        # exons are normally already disjoint per transcript).
        merged = []
        for s, e in exons:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        exonic_length = sum(e - s for s, e in merged)
        span_start = merged[0][0]
        span_end = merged[-1][1]
        result[gene] = (best["chrom"], span_start, span_end, merged, exonic_length)

    return result


def _read_ref_blocks(read):
    """Reference-coordinate (0-based half-open) blocks covered by this
    read's non-N, reference-consuming CIGAR ops (M/D/=/X)."""
    blocks = []
    pos = read.reference_start
    block_start = None
    for op, length in read.cigartuples or []:
        if op in _REF_CONSUMING_NON_N:
            if block_start is None:
                block_start = pos
            pos += length
        else:
            if block_start is not None:
                blocks.append((block_start, pos))
                block_start = None
            if op == 3:  # N: intron skip, still consumes reference
                pos += length
    if block_start is not None:
        blocks.append((block_start, pos))
    return blocks


def _overlap_length(blocks_a, blocks_b):
    """Total overlap length between two lists of disjoint, sorted (start,end)
    interval tuples."""
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


def compute_avg_flr(bam_handle, chrom, span_start, span_end, exons, exonic_length):
    """Returns (avg_flr, n_reads), fetching from an already-open
    AlignmentFile via an indexed region query. avg_flr is NaN when there's
    no usable data; n_reads is always a real count (0 when there's
    nothing)."""
    if exonic_length <= 0:
        return np.nan, 0
    ratios = []
    for read in bam_handle.fetch(chrom, span_start, span_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        read_blocks = _read_ref_blocks(read)
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
