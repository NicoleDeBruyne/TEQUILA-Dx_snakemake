
import argparse
import gzip
import os

import pandas as pd
import numpy as np
import pysam


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
    genes = set(genes)

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
            else:
                t = transcripts.setdefault(transcript_id, {
                    "gene": gene_name, "chrom": chrom, "strand": strand,
                    "canonical": False, "exons": [],
                })
                t["exons"].append((int(start) - 1, int(end)))

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
            if op == 3:
                pos += length
    if block_start is not None:
        blocks.append((block_start, pos))
    return blocks


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


def compute_avg_flr(bam_handle, chrom, span_start, span_end, exons, exonic_length):
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
