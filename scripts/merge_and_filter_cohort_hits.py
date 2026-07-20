#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2026.07.06
#
# Orchestrates the cross-sample merge-and-filter stage for a single (bed, sample_type)
# group: merges & filters variant calls, ASE results, and outlier splice junctions across
# all samples in the group, then merges those three hit types per sample and produces one
# cohort-level candidate-hits table. This replaces the old merge_candidate_hits.sh, which
# did the same thing via bash + file globbing; here the sample/tissue/file associations
# are supplied directly by Snakemake instead of being re-derived by globbing directories.
#
# Reuses (via subprocess) the existing, unmodified:
#   merge_and_filter_variants.py, merge_and_filter_ase_results.py,
#   merge_and_filter_junction_results.py, merge_hits.py, plot_candidate_hits.py

import argparse
import os
import shutil
import subprocess
import sys
from math import ceil

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge and filter variant, ASE, and splice-junction outlier hits "
                    "across all samples in a (bed, sample_type) group.")
    parser.add_argument("--group-id", required=True, help="Identifier for this group (used in plot titles).")
    parser.add_argument("--outdir", required=True, help="Output directory for this group's merged results.")
    parser.add_argument("--omim", required=True, help="Path to OMIM data (columns: approved_gene_symbol, phenotypes, inheritance_patterns).")

    parser.add_argument("--samples", nargs="+", required=True, help="Sample names in this group.")
    parser.add_argument("--variant-files", nargs="+", required=True, help="Per-sample *_filtered_variants.tsv paths, aligned with --samples.")
    parser.add_argument("--ase-files", nargs="+", required=True, help="Per-sample *_binomial_ase_results.tsv paths, aligned with --samples.")
    parser.add_argument("--junction-manifest", nargs="*", default=[], help="Entries of 'tissue:::sample:::filepath' for each (sample, tissue) all_junctions.tsv.")

    # Variant merge thresholds
    parser.add_argument("--num-callers-threshold-snv", type=int, default=2)
    parser.add_argument("--num-callers-threshold-indel", type=int, default=2)
    parser.add_argument("--min-dp-snv", type=int, default=20)
    parser.add_argument("--min-dp-indel", type=int, default=20)
    parser.add_argument("--variant-sample-fraction", type=float, default=0.05)

    # ASE merge thresholds
    parser.add_argument("--min-haplotype-ratio", type=float, default=0.1)
    parser.add_argument("--delta-haplotype-ratio-threshold", type=float, default=0.2)
    parser.add_argument("--ase-padj-threshold", type=float, default=0.001)
    parser.add_argument("--ase-sample-fraction", type=float, default=0.05)

    # Junction merge thresholds
    parser.add_argument("--jxn-coverage-threshold", type=int, default=50)
    parser.add_argument("--jxn-padj-threshold", type=float, default=0.001)
    parser.add_argument("--delta-psi-threshold", type=float, default=0.2)
    parser.add_argument("--jxn-sample-fraction", type=float, default=0.01)

    return parser.parse_args()


def run(cmd):
    print(f"\n+ {' '.join(str(c) for c in cmd)}\n", flush=True)
    subprocess.run(cmd, check=True)


def merge_variants(args):
    """Merge and filter variant calls across all samples. Returns path to the merged
    hits tsv, or None if no input files were available."""

    os.makedirs(f"{args.outdir}/variant_calling", exist_ok=True)
    input_files = []
    for sample, variant_file in zip(args.samples, args.variant_files):
        if not os.path.isfile(variant_file):
            print(f"WARNING: No filtered variants file found for sample {sample} at {variant_file}", file=sys.stderr)
            continue
        input_files.append(variant_file)

    if not input_files:
        print("\nNo input files found for variant calling analysis")
        return None

    print("\nMerging and filtering variant hits...\n")
    n = ceil(len(input_files) * args.variant_sample_fraction)
    outprefix = f"{args.outdir}/variant_calling/all_candidate_variants_{n}samples"
    run([
        sys.executable, f"{SCRIPT_DIR}/merge_and_filter_variants.py",
        "--infiles", *input_files,
        "--outprefix", outprefix,
        "--num-callers-threshold-SNV", str(args.num_callers_threshold_snv),
        "--num-callers-threshold-indel", str(args.num_callers_threshold_indel),
        "--min-DP-SNV", str(args.min_dp_snv),
        "--min-DP-indel", str(args.min_dp_indel),
        "--sample-number-threshold", str(n),
        "--plot",
        "--plot-variant-type", "SNV", "indel",
        "--title", f"{args.group_id} Variant Counts",
    ])
    return f"{outprefix}.tsv"


def merge_ase(args):
    """Merge and filter ASE results across all samples. Returns path to the merged
    hits tsv, or None if no input files were available."""

    os.makedirs(f"{args.outdir}/ase_analysis", exist_ok=True)
    input_files = [f for f in args.ase_files if os.path.isfile(f)]
    missing = set(args.ase_files) - set(input_files)
    for f in missing:
        print(f"WARNING: No ASE results file found at {f}", file=sys.stderr)

    if not input_files:
        print("\nNo input files found for ASE analysis")
        return None

    print("\n\n\nMerging and filtering ASE hits...\n")
    n = ceil(len(input_files) * args.ase_sample_fraction)
    outprefix = (f"{args.outdir}/ase_analysis/outlier_ase_{args.ase_padj_threshold}padj_"
                 f"{args.delta_haplotype_ratio_threshold}diff_{n}samples")
    run([
        sys.executable, f"{SCRIPT_DIR}/merge_and_filter_ase_results.py",
        "--infiles", *input_files,
        "--outprefix", outprefix,
        "--min-haplotype-ratio", str(args.min_haplotype_ratio),
        "--delta-haplotype-ratio-threshold", str(args.delta_haplotype_ratio_threshold),
        "--padj-threshold", str(args.ase_padj_threshold),
        "--plot",
        "--title", f"{args.group_id}: Number of Genes with Allele-specific Expression by Sample",
        "--sample-number-threshold", str(n),
    ])
    return f"{outprefix}.tsv"


def merge_junctions(args):
    """Merge and filter outlier junctions across all samples, once per tissue.
    Returns {tissue: merged_hits_tsv_path} for tissues that had input files."""

    manifest = []
    for entry in args.junction_manifest:
        tissue, sample, filepath = entry.split(":::", 2)
        manifest.append((tissue, sample, filepath))

    tissues = sorted({t for t, _, _ in manifest})
    results = {}

    for tissue in tissues:
        input_files = [fp for t, _, fp in manifest if t == tissue and os.path.isfile(fp)]
        if not input_files:
            print(f"\n\n\nNo input files found for junction analysis in tissue: {tissue}")
            continue

        print(f"\n\n\nMerging and filtering junction hits in tissue: {tissue}...\n")
        outdir = f"{args.outdir}/junction_analysis/gtex_{tissue}"
        os.makedirs(outdir, exist_ok=True)
        outprefix = f"{outdir}/outlier_junctions_gtex_{tissue}"
        n = ceil(len(input_files) * args.jxn_sample_fraction)

        run([
            sys.executable, f"{SCRIPT_DIR}/merge_and_filter_junction_results.py",
            "--infiles", *input_files,
            "--outprefix", outprefix,
            "--jxn-coverage-threshold", str(args.jxn_coverage_threshold),
            "--padj-threshold", str(args.jxn_padj_threshold),
            "--delta-PSI-threshold", str(args.delta_psi_threshold),
            "--event-types", "exon_skipping", "exon_inclusion", "alt_ss1", "alt_ss2",
            "--sample-number-threshold", str(n),
            "--filter-by-cohort-IQR",
            "--plot",
            "--title", f"{args.group_id}: Number of Genes with Outlier Junctions by Sample",
        ])

        # merge_and_filter_junction_results.py only appends a "_cohortIQR" suffix if it
        # actually computed cohort statistics, which requires >= 8 samples (see that
        # script's internal check) -- mirror that here when predicting the output name.
        base = f"{outprefix}_{args.jxn_coverage_threshold}jxncov_{args.jxn_padj_threshold}padj_{args.delta_psi_threshold}deltaPSI_event"
        if len(input_files) >= 8:
            src = f"{base}_cohortIQR_{n}samples.tsv"
        else:
            src = f"{base}_{n}samples.tsv"
        final = f"{outprefix}_final.tsv"
        shutil.copyfile(src, final)
        results[tissue] = final

    return results


def split_by_sample(tsv_path, sample_col, samples, default_header, outdir, suffix):
    """Split a merged hits tsv into one file per sample. Samples with no rows (or if
    tsv_path is None) get a header-only stub with default_header."""

    os.makedirs(outdir, exist_ok=True)
    per_sample = {}

    if tsv_path is not None and os.path.isfile(tsv_path):
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    else:
        df = pd.DataFrame(columns=default_header)

    for sample in samples:
        out_path = f"{outdir}/{sample}_{suffix}.tsv"
        if sample_col in df.columns:
            sub = df[df[sample_col].astype(str) == str(sample)]
        else:
            sub = pd.DataFrame(columns=default_header)
        if sub.empty:
            pd.DataFrame(columns=default_header).to_csv(out_path, sep="\t", index=False)
        else:
            sub.to_csv(out_path, sep="\t", index=False)
        per_sample[sample] = out_path

    return per_sample


def split_junctions_by_sample(tissue_results, samples, default_header, outdir):
    """Concatenate each sample's rows across all tissues into one junction-hits file per sample."""

    os.makedirs(outdir, exist_ok=True)
    per_tissue_dfs = []
    for tissue, tsv_path in tissue_results.items():
        if os.path.isfile(tsv_path):
            per_tissue_dfs.append(pd.read_csv(tsv_path, sep="\t", dtype=str))

    combined = pd.concat(per_tissue_dfs, ignore_index=True) if per_tissue_dfs else pd.DataFrame(columns=default_header)

    per_sample = {}
    for sample in samples:
        out_path = f"{outdir}/{sample}_junction_hits.tsv"
        if "sample" in combined.columns:
            sub = combined[combined["sample"].astype(str) == str(sample)]
        else:
            sub = pd.DataFrame(columns=default_header)
        if sub.empty:
            pd.DataFrame(columns=default_header).to_csv(out_path, sep="\t", index=False)
        else:
            sub.to_csv(out_path, sep="\t", index=False)
        per_sample[sample] = out_path

    return per_sample


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"\n\n\n******************************************************************************************")
    print(f"Merging hits for group: {args.group_id}")
    print(f"******************************************************************************************\n")

    all_variant_hits = merge_variants(args)
    all_ase_hits = merge_ase(args)
    tissue_junction_hits = merge_junctions(args)

    print("\n\n\nMerging hits for each sample, using the following files:")
    print(f"Variant hits: {all_variant_hits}")
    print(f"ASE hits: {all_ase_hits}")
    for tissue, f in tissue_junction_hits.items():
        print(f"Junction hits ({tissue}): {f}")

    merged_hits_dir = f"{args.outdir}/merged_hits"
    by_sample_dir = f"{merged_hits_dir}/by_sample"
    os.makedirs(by_sample_dir, exist_ok=True)

    print("\nSplitting variant hits by sample...")
    variant_default_header = ['sample', 'chrom', 'pos', 'ref', 'alt', 'GT', 'gnomAD_AF',
                               'CLNSIG', 'gene', 'CADD_PHRED', 'SpliceAI',
                               'num_callers', 'sample_count']
    variant_by_sample = split_by_sample(all_variant_hits, "sample", args.samples,
                                         variant_default_header, by_sample_dir, "variant_hits")

    print("Splitting ASE hits by sample...")
    ase_default_header = ['sample', 'gene', 'ratio', 'sample_count']
    ase_by_sample = split_by_sample(all_ase_hits, "sample", args.samples,
                                     ase_default_header, by_sample_dir, "ase_hits")

    print("Splitting junction hits by sample...")
    junction_default_header = ['sample', 'phasing', 'gene', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count']
    junction_by_sample = split_junctions_by_sample(tissue_junction_hits, args.samples,
                                                    junction_default_header, by_sample_dir)

    print("Merging hits per sample...")
    for sample in args.samples:
        run([
            sys.executable, f"{SCRIPT_DIR}/merge_hits.py",
            "--outfile", f"{by_sample_dir}/{sample}_hits.tsv",
            "--sample-name", sample,
            "--variant-hits", variant_by_sample[sample],
            "--ase-hits", ase_by_sample[sample],
            "--junction-hits", junction_by_sample[sample],
            "--omim", args.omim,
        ])

    print("\nMerging all hits into one file...")
    dfs = [pd.read_csv(f"{by_sample_dir}/{sample}_hits.tsv", sep="\t", dtype=str) for sample in args.samples]
    all_hits = pd.concat(dfs, ignore_index=True)
    all_hits_path = f"{merged_hits_dir}/all_hits.tsv"
    all_hits.to_csv(all_hits_path, sep="\t", index=False)
    print(f"Wrote {all_hits_path}")

    run([
        sys.executable, f"{SCRIPT_DIR}/plot_candidate_hits.py",
        "--infile", all_hits_path,
        "--outdir", merged_hits_dir,
    ])


if __name__ == "__main__":
    main()
