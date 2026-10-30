# 3C - check how well the RNA read phasing agrees with the sample's own WGS phasing
# (only for samples that have a WGS VCF)

import argparse
import itertools
import os

import pandas as pd
import pysam

BASES = "ACGT"
SITE_COLS = ["sample", "mode", "gene", "chrom", "pos", "ref", "alt", "wgs_phased", "wgs_block", "wgs_alt_hap",
             "hap1_ref", "hap1_alt", "hap2_ref", "hap2_alt", "rna_alt_hap"]
GENE_COLS = ["sample", "mode", "gene", "wgs_het_sites", "wgs_phased_sites", "assessed_sites", "pairs",
             "concordant_pairs", "adjacent_pairs", "adjacent_discordant"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare the haplotype BAMs from read phasing with the phase of heterozygous SNVs in the "
                    "sample's WGS VCF. For each WGS het site the alt allele is placed on hap1 or hap2 from the "
                    "reads, and every pair of sites in the same WGS phase block is scored as agreeing or not "
                    "on cis vs trans. Hap labels are arbitrary per gene, so only relative phase is compared.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--wgs-vcf", required=True,
        help="Panel WGS VCF from _1E_prepare_wgs_vcf.py. Use this one, not the 1F output: 1F already "
             "mixes in the NanoTS phase blocks, so it isn't an independent truth set.")
    parser.add_argument("--mode", nargs="+", required=True, metavar="NAME=MAPPING",
        help="One or more phasing runs to check, as name=gene_bam_mapping_file.tsv "
             "(e.g. wgs_blinded=... wgs_aided=...).")
    parser.add_argument("--min-depth", type=int, default=5,
        help="Reads needed on each haplotype at a site before we trust where its alt allele went.")
    parser.add_argument("--min-frac", type=float, default=0.8,
        help="Alt fraction needed on the alt haplotype (and at most 1 - this on the other one).")
    parser.add_argument("--outprefix", required=True,
        help="Writes <outprefix>_sites.tsv and <outprefix>_genes.tsv.")
    return parser.parse_args()


def load_wgs_hets(path):
    sites = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            if len(rec.ref) != 1 or len(rec.alts or ()) != 1 or len(rec.alts[0]) != 1:
                continue  # SNVs only, pileups over indels are too fiddly to trust
            call = rec.samples[0]
            gt = call.get("GT")
            if gt is None or sorted(gt) != [0, 1]:
                continue
            ps = call.get("PS") if "PS" in call else None
            sites.append(dict(chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=rec.alts[0],
                              wgs_phased=bool(call.phased),
                              # no PS on a phased call means the whole chromosome is one block
                              wgs_block=f"{rec.chrom}:{ps}" if ps is not None else rec.chrom,
                              wgs_alt_hap=gt.index(1)))
    return pd.DataFrame(sites)


def allele_counts(bam, chrom, pos, ref, alt):
    cov = bam.count_coverage(chrom, pos - 1, pos, quality_threshold=0, read_callback="all")
    counts = {b: int(cov[i][0]) for i, b in enumerate(BASES)}
    return counts.get(ref, 0), counts.get(alt, 0)


def rna_alt_hap(r1, a1, r2, a2, min_depth, min_frac):
    if r1 + a1 < min_depth or r2 + a2 < min_depth:
        return None
    f1, f2 = a1 / (r1 + a1), a2 / (r2 + a2)
    if f1 >= min_frac and f2 <= 1 - min_frac:
        return 0
    if f2 >= min_frac and f1 <= 1 - min_frac:
        return 1
    return None  # alt on both or neither, can't call it


def check_mode(mode, mapping_path, hets, args):
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    site_rows, gene_rows = [], []
    for _, m in mapping.iterrows():
        chrom, span = m["region"].split(":")
        start, end = (int(x) for x in span.split("-"))
        in_gene = hets[(hets["chrom"] == chrom) & (hets["pos"] > start) & (hets["pos"] <= end)]
        assessed = []
        if len(in_gene):
            with pysam.AlignmentFile(m["hap1_bam"]) as b1, pysam.AlignmentFile(m["hap2_bam"]) as b2:
                for _, s in in_gene.iterrows():
                    r1, a1 = allele_counts(b1, chrom, s["pos"], s["ref"], s["alt"])
                    r2, a2 = allele_counts(b2, chrom, s["pos"], s["ref"], s["alt"])
                    hap = rna_alt_hap(r1, a1, r2, a2, args.min_depth, args.min_frac)
                    site_rows.append(dict(sample=args.sample, mode=mode, gene=m["gene"], chrom=chrom, pos=s["pos"],
                                          ref=s["ref"], alt=s["alt"], wgs_phased=s["wgs_phased"],
                                          wgs_block=s["wgs_block"], wgs_alt_hap=s["wgs_alt_hap"] + 1,
                                          hap1_ref=r1, hap1_alt=a1, hap2_ref=r2, hap2_alt=a2,
                                          rna_alt_hap="" if hap is None else hap + 1))
                    if hap is not None and s["wgs_phased"]:
                        assessed.append((s["wgs_block"], s["pos"], s["wgs_alt_hap"], hap))

        # every pair in the same WGS block: do RNA and WGS agree on cis (same hap) vs trans?
        pairs = agree = adjacent = adjacent_bad = 0
        for _, block in itertools.groupby(sorted(assessed), key=lambda x: x[0]):
            block = list(block)
            for a, b in itertools.combinations(block, 2):
                pairs += 1
                agree += (a[2] == b[2]) == (a[3] == b[3])
            # neighbouring sites only, closer to the usual switch error idea
            for a, b in zip(block, block[1:]):
                adjacent += 1
                adjacent_bad += (a[2] == b[2]) != (a[3] == b[3])
        gene_rows.append(dict(sample=args.sample, mode=mode, gene=m["gene"], wgs_het_sites=len(in_gene),
                              wgs_phased_sites=int(in_gene["wgs_phased"].sum()) if len(in_gene) else 0,
                              assessed_sites=len(assessed), pairs=pairs, concordant_pairs=agree,
                              adjacent_pairs=adjacent, adjacent_discordant=adjacent_bad))
    return site_rows, gene_rows


def main():
    print("\n\n\n******************************************************************************************")
    print("Checking read phasing against WGS...")
    print("******************************************************************************************\n")
    args = parse_args()

    hets = load_wgs_hets(args.wgs_vcf)
    if hets.empty:
        print(f"WARNING: no heterozygous SNVs in {args.wgs_vcf}, nothing to compare")
        hets = pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "wgs_phased", "wgs_block", "wgs_alt_hap"])
    else:
        print(f"{len(hets)} WGS het SNVs, {int(hets['wgs_phased'].sum())} of them phased")
        if not hets["wgs_phased"].any():
            print("WARNING: the WGS VCF has no phased calls, so pairs can't be scored")

    all_sites, all_genes = [], []
    for spec in args.mode:
        mode, _, path = spec.partition("=")
        if not path or not os.path.isfile(path):
            raise SystemExit(f"ERROR: --mode {spec!r}: expected name=existing mapping file")
        sites, genes = check_mode(mode, path, hets, args)
        all_sites += sites
        all_genes += genes
        g = pd.DataFrame(genes, columns=GENE_COLS)
        n_pairs = int(g["pairs"].sum()) if len(g) else 0
        rate = g["concordant_pairs"].sum() / n_pairs if n_pairs else float("nan")
        print(f"{mode}: {len(g)} phased genes, {int(g['assessed_sites'].sum()) if len(g) else 0} sites placed on a "
              f"haplotype, {n_pairs} pairs, {rate:.3f} concordant")

    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)), exist_ok=True)
    pd.DataFrame(all_sites, columns=SITE_COLS).to_csv(f"{args.outprefix}_sites.tsv", sep="\t", index=False)
    pd.DataFrame(all_genes, columns=GENE_COLS).to_csv(f"{args.outprefix}_genes.tsv", sep="\t", index=False)
    print(f"\nWrote {args.outprefix}_sites.tsv and {args.outprefix}_genes.tsv\nDone.")


if __name__ == "__main__":
    main()
