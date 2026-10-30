# 10P - sanity check that each RNA sample matches its WGS, using high-confidence NanoTS SNVs in the panel

import argparse
import os
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pysam
from palette import use_fonts

use_fonts()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check that each sample's RNA-seq comes from the same person as its WGS. High-confidence "
                    "NanoTS SNVs from every sample are compared with every WGS VCF in the cohort; the score is "
                    "the fraction of a sample's NanoTS SNVs whose alt allele is also present in the WGS VCF.")
    parser.add_argument("--samples", nargs="+", required=True, help="All samples on the panel.")
    parser.add_argument("--nanoTS-vcfs", nargs="+", required=True,
        help="Normalized NanoTS VCFs, in the same order as --samples.")
    parser.add_argument("--wgs-samples", nargs="+", required=True, help="Samples that have WGS.")
    parser.add_argument("--wgs-vcfs", nargs="+", required=True,
        help="Panel-restricted WGS VCFs from _1E_prepare_wgs_vcf.py, in the same order as --wgs-samples.")
    parser.add_argument("--bed", required=True, help="Panel BED (chrom, start, end, gene).")
    parser.add_argument("--min-dp", type=int, default=20,
        help="Minimum NanoTS depth (DP) for an SNV to count as high-confidence.")
    parser.add_argument("--min-sites", type=int, default=20,
        help="Samples with fewer high-confidence NanoTS SNVs than this are flagged too_few_sites.")
    parser.add_argument("--threshold", type=float, default=0.8,
        help="Concordance with the sample's own WGS below this is flagged low; a sample without WGS whose best "
             "match reaches this is flagged matches_other_wgs.")
    parser.add_argument("--outprefix", required=True,
        help="Output prefix for _matrix.tsv, _summary.tsv and _heatmap.pdf.")
    return parser.parse_args()


def load_regions(bed):
    by_chrom = {}
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            by_chrom.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    merged = []
    for chrom, intervals in by_chrom.items():
        cur = None
        for start, end in sorted(intervals):
            if cur and start <= cur[1]:
                cur[1] = max(cur[1], end)
                continue
            if cur:
                merged.append((chrom, cur[0], cur[1]))
            cur = [start, end]
        if cur:
            merged.append((chrom, cur[0], cur[1]))
    return merged


def fetch_records(vcf, regions):
    contigs = set(vcf.header.contigs)
    for chrom, start, end in regions:
        if chrom not in contigs:
            continue
        for rec in vcf.fetch(chrom, start, end):
            if start < rec.pos <= end and rec.alts:
                yield rec


def has_alt(sample):
    gt = sample.get("GT")
    return bool(gt) and None not in gt and any(a > 0 for a in gt)


def nanoTS_snvs(path, regions, min_dp):
    kept = []
    with pysam.VariantFile(path) as vcf:
        for rec in fetch_records(vcf, regions):
            if list(rec.filter.keys()) != ["PASS"] or len(rec.alts) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue
            s = rec.samples[0]
            if not has_alt(s) or (s.get("DP") or 0) < min_dp:
                continue
            kept.append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
    per_pos = Counter((c, p) for c, p, _, _ in kept)
    return {k for k in kept if per_pos[(k[0], k[1])] == 1}


def wgs_alleles(path, regions):
    alleles = set()
    with pysam.VariantFile(path) as vcf:
        for rec in fetch_records(vcf, regions):
            if has_alt(rec.samples[0]):
                alleles.add((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
    return alleles


# ok / low / mismatch (some other WGS fits better) / matches a different WGS
def summarize(samples, wgs_samples, rna_sites, scores, args):
    rows = []
    for s in samples:
        n = len(rna_sites[s])
        row_scores = scores.loc[s]
        best = row_scores.idxmax() if n else None
        best_score = row_scores.max() if n else None
        own = row_scores.get(s) if s in wgs_samples and n else None
        if n < args.min_sites:
            flag = "too_few_sites"
        elif s in wgs_samples:
            if best != s and best_score > own:
                flag = "mismatch"
            elif own < args.threshold:
                flag = "low"
            else:
                flag = "ok"
        else:
            flag = "matches_other_wgs" if best_score >= args.threshold else "no_wgs"
        rows.append(dict(sample=s, has_wgs=s in wgs_samples, n_rna_snvs=n,
                         own_wgs_concordance=round(own, 4) if own is not None else ".",
                         best_wgs_match=best if best is not None else ".",
                         best_wgs_concordance=round(best_score, 4) if best_score is not None else ".",
                         flag=flag))
    return pd.DataFrame(rows)


def plot_heatmap(scores, outfile):
    n_rows, n_cols = scores.shape
    fig, ax = plt.subplots(figsize=(max(5.0, 1.6 + 0.32 * n_cols), max(3.6, 1.2 + 0.22 * n_rows)))
    im = ax.imshow(scores.values.astype(float), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(scores.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(scores.index, fontsize=6)
    ax.set_xlabel("WGS VCF")
    ax.set_ylabel("RNA-seq sample (NanoTS SNVs)")
    col_index = {w: j for j, w in enumerate(scores.columns)}
    for i, s in enumerate(scores.index):
        if s in col_index:
            ax.add_patch(plt.Rectangle((col_index[s] - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="#d62728", lw=1.2))
    if n_rows * n_cols <= 900:
        for i in range(n_rows):
            for j in range(n_cols):
                v = scores.iat[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.5,
                        color="black" if v > 0.6 else "white")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Fraction of NanoTS SNVs found in WGS", fontsize=6)
    cbar.ax.tick_params(labelsize=5)
    ax.set_title("RNA vs WGS concordance (red box = sample's own WGS)", fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)


def main():
    print("\n\n\n******************************************************************************************")
    print("Checking RNA vs WGS concordance...")
    print("******************************************************************************************\n")

    args = parse_args()
    if len(args.samples) != len(args.nanoTS_vcfs):
        raise ValueError("--samples and --nanoTS-vcfs must have the same length.")
    if len(args.wgs_samples) != len(args.wgs_vcfs):
        raise ValueError("--wgs-samples and --wgs-vcfs must have the same length.")

    regions = load_regions(args.bed)
    rna_sites = {}
    for s, path in zip(args.samples, args.nanoTS_vcfs):
        rna_sites[s] = nanoTS_snvs(path, regions, args.min_dp)
        print(f"{s}: {len(rna_sites[s])} high-confidence NanoTS SNV(s)")
    wgs = {}
    for s, path in zip(args.wgs_samples, args.wgs_vcfs):
        wgs[s] = wgs_alleles(path, regions)
        print(f"{s} WGS: {len(wgs[s])} variant allele(s) in the panel gene regions")

    # concordance = fraction of the sample's RNA SNVs that show up in each WGS
    scores = pd.DataFrame(
        [[len(rna_sites[s] & wgs[w]) / len(rna_sites[s]) if rna_sites[s] else 0.0 for w in args.wgs_samples]
         for s in args.samples],
        index=args.samples, columns=args.wgs_samples)

    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)), exist_ok=True)
    scores.round(4).to_csv(args.outprefix + "_matrix.tsv", sep="\t", index_label="sample")
    summary = summarize(args.samples, set(args.wgs_samples), rna_sites, scores, args)
    summary.to_csv(args.outprefix + "_summary.tsv", sep="\t", index=False)
    plot_heatmap(scores, args.outprefix + "_heatmap.pdf")

    print("\nSummary:")
    print(summary.to_string(index=False))
    flagged = summary[summary["flag"].isin(["mismatch", "low", "matches_other_wgs"])]
    if not flagged.empty:
        print(f"\nWARNING: {len(flagged)} sample(s) flagged: "
              + ", ".join(f"{r.sample} ({r.flag})" for r in flagged.itertuples()))
    print(f"\nWrote {args.outprefix}_matrix.tsv, {args.outprefix}_summary.tsv and {args.outprefix}_heatmap.pdf")


if __name__ == "__main__":
    main()
