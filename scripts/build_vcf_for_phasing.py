#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2025.10.30
# Optimized: 2025

import argparse
import os
import pysam
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description="Extracts biallelic PASS SNVs from nanoTS (DP≥20, AF≥0.1) \
                                                    and PASS indels called by both Clair3 and DeepVariant (DP≥20, AF≥0.1).")
    parser.add_argument("--longcallR-vcf", required=True, help="LongcallR VCF (SNVs only).")
    parser.add_argument("--nanoTS-vcf", required=True, help="NanoTS VCF (SNVs only).")
    parser.add_argument("--clair3-vcf", required=True, help="Clair3 VCF.")
    parser.add_argument("--deepvariant-vcf", required=True, help="DeepVariant VCF.")
    parser.add_argument("--outfile", required=True, help="Output merged VCF file.")
    parser.add_argument("--sample-name", default="SAMPLE", help="Sample name to use in output VCF.")
    parser.add_argument("--bcftools-exec", default="bcftools", help="Path to bcftools.")
    parser.add_argument("--tabix-exec", default="tabix", help="Path to tabix.")
    return parser.parse_args()

def get_variants(vcf_path, indels_only=False):
    """Return a dict {(chrom,pos,ref,alt): (GT,DP,AF)} filtered for PASS, DP>=20, AF>=0.1."""

    variants = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.fetch():
            # Fast filter: PASS only, biallelic only
            if rec.filter.keys() != ["PASS"]:
                continue
            if len(rec.alts) != 1:
                continue
            ref, alt = rec.ref, rec.alts[0]
            # Skip SNVs when indels_only is requested
            if indels_only and len(ref) == 1 and len(alt) == 1:
                continue

            # Access first (and only) sample directly by index — faster than list()
            s = rec.samples[0]
            DP = s.get("DP") or 0
            AF = s.get("AF") or s.get("VAF") or (0,)
            if isinstance(AF, (list, tuple)):
                AF = AF[0]
            if DP < 20 or AF < 0.1:
                continue
            variants[(rec.chrom, rec.pos, ref, alt)] = (s.get("GT"), DP, AF)

    return variants

def write_vcf(outfile, variants, sample_name):
    """Write minimal VCF."""

    # Single pass over sorted items, collect chroms on the fly
    sorted_variants = sorted(variants.items())
    # Collect unique chroms in sorted order without a second pass
    seen = set()
    chroms = []
    for (chrom, _, _, _), _ in sorted_variants:
        if chrom not in seen:
            seen.add(chrom)
            chroms.append(chrom)

    lines = [
        "##fileformat=VCFv4.2\n",
        *(f"##contig=<ID={c}>\n" for c in chroms),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n",
    ]
    for (chrom, pos, ref, alt), (GT, DP, AF) in sorted_variants:
        gt_str = "/".join(map(str, GT)) if GT else "./."
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP:AF\t{gt_str}:{DP}:{AF:.3f}\n")

    with open(outfile, "w") as f:
        f.writelines(lines)

def main():
    """ Main function """

    # Parse arguments
    args = parse_args()

    # Load longcallR and NanoTS SNVs
    print("Loading biallelic PASS SNVs from longcallR...")
    longcallR_SNVs = get_variants(args.longcallR_vcf)
    nanoTS_SNVs = get_variants(args.nanoTS_vcf)
    print(f"  → Found {len(longcallR_SNVs)} SNVs from longcallR VCF that pass filters")
    print(f"  → Found {len(nanoTS_SNVs)} SNVs from nanoTS VCF that pass filters")

    # Load Clair3 and DeepVariant indels
    print("Loading biallelic PASS indels from Clair3 and DeepVariant...")
    clair3_indels = get_variants(args.clair3_vcf, indels_only=True)
    deepvariant_indels = get_variants(args.deepvariant_vcf, indels_only=True)
    # Dict comprehension is faster than iterating keys when deepvariant_indels is also a dict
    shared_indels = {k: clair3_indels[k] for k in clair3_indels.keys() & deepvariant_indels.keys()}
    print(f"  → Found {len(shared_indels)} indels shared by Clair3 and DeepVariant that pass filters")

    # Merge variants (nanoTS SNVs + shared indels)
    merged = {**nanoTS_SNVs, **shared_indels}

    # Write output VCF
    outdir = os.path.dirname(args.outfile)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = args.outfile
    if outfile.endswith(".gz"):
        outfile = outfile[:-3]
    write_vcf(outfile, merged, args.sample_name)
    print(f"\nTotal variants written: {len(merged)}")

    # bgzip + tabix
    subprocess.run([args.bcftools_exec, "sort", "-Oz", "-o", f"{outfile}.gz", outfile], check=True)
    os.remove(outfile)
    subprocess.run([args.tabix_exec, "-f", "-p", "vcf", f"{outfile}.gz"], check=True)
    print(f"\n✅ Output written to {outfile}.gz and indexed.")

if __name__ == "__main__":
    main()
