
import argparse
import gzip
import re
from collections import defaultdict

import pandas as pd
import pysam



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



def load_targeted_genes(bed):
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

    open_fn = gzip.open if gtf_path.endswith(".gz") else open

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
            start0 = int(parts[3]) - 1
            end = int(parts[4])
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
                ss1 = end_e + 1
                ss2 = start_n
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
    for e_start, e_end in exons:
        if e_start >= end:
            break
        if e_end > start:
            return True
    return False


def _contained_in_exons(exons, start, end):
    for e_start, e_end in exons:
        if e_start <= start < e_end:
            return end <= e_end
        if e_start > start:
            break
    return False



_CIGAR_CONSUMES_REF = {0, 2, 3, 7, 8}


def read_splice_sites(read):
    sites = set()
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op == 3:
            sites.add(ref_pos + 1)
            sites.add(ref_pos + length)
            ref_pos += length
        elif op in _CIGAR_CONSUMES_REF:
            ref_pos += length
    return sites


def _keep_read(r):
    return not r.is_unmapped and not r.is_secondary



def assign_sample(bam, genes, gene_bins):

    counts = defaultdict(int)
    stats = {
        "n_total": 0,
        "spliced_assigned": 0,
        "spliced_unassigned_zero_shared": 0,
        "spliced_unassigned_tied": 0,
        "unspliced_assigned": 0,
        "unspliced_unassigned_zero_overlap": 0,
        "unspliced_unassigned_multi_overlap": 0,
        "unspliced_unassigned_not_contained": 0,
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
