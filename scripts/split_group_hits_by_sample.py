#!/usr/bin/env python3

# Splits a (bed, sample_type) group's merged variant, ASE, and junction hit
# tables into one file per sample, ready for merge_hits.py to combine
# per-sample. Fully standalone -- reads its inputs from disk, doesn't shell
# out to or import any other pipeline script.
#
# Samples with no matching rows in a given merged table (e.g. no variant
# hits at all) get a header-only stub file, using the appropriate default
# header, so downstream per-sample merging always has a well-formed file to
# read regardless of whether that sample had any hits.

import argparse
import os

import pandas as pd

VARIANT_HEADER = ['sample', 'chrom', 'pos', 'ref', 'alt', 'GT', 'gnomAD_AF',
                   'CLNSIG', 'gene', 'CADD_PHRED', 'SpliceAI',
                   'num_callers', 'sample_count']
ASE_HEADER = ['sample', 'gene', 'ratio', 'sample_count']
JUNCTION_HEADER = ['sample', 'phasing', 'gene', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count']


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split a group's merged variant/ASE/junction hit tables into per-sample files.")
    parser.add_argument("--variant-tsv", required=True,
        help="Group-level merged variant hits (output of merge_and_filter_variants.py).")
    parser.add_argument("--ase-tsv", required=True,
        help="Group-level merged ASE hits (output of merge_and_filter_ase_results.py).")
    parser.add_argument("--tissues", nargs="*", default=[],
        help="Tissues present in this group (same order as --junction-files).")
    parser.add_argument("--junction-files", nargs="*", default=[],
        help="Per-tissue merged junction hit files (output of merge_and_filter_junction_results.py, "
             "one per --tissues entry). Concatenated across tissues before splitting by sample.")
    parser.add_argument("--samples", nargs="+", required=True, help="Sample names in this group.")
    parser.add_argument("--outdir", required=True, help="Output directory for per-sample stub files.")
    args = parser.parse_args()
    if len(args.tissues) != len(args.junction_files):
        parser.error("--tissues and --junction-files must have the same number of entries")
    return args


def split_by_sample(df, sample_col, samples, default_header, outdir, suffix):
    """Write one {outdir}/{sample}_{suffix}.tsv per sample. Samples with no
    matching rows get a header-only stub with default_header."""
    os.makedirs(outdir, exist_ok=True)
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
        print(f"  {sample}: {len(sub)} row(s) -> {out_path}")


def main():
    args = parse_args()

    print("\nSplitting variant hits by sample...")
    variant_df = (pd.read_csv(args.variant_tsv, sep="\t", dtype=str)
                  if os.path.isfile(args.variant_tsv) else pd.DataFrame(columns=VARIANT_HEADER))
    split_by_sample(variant_df, "sample", args.samples, VARIANT_HEADER, args.outdir, "variant_hits")

    print("\nSplitting ASE hits by sample...")
    ase_df = (pd.read_csv(args.ase_tsv, sep="\t", dtype=str)
              if os.path.isfile(args.ase_tsv) else pd.DataFrame(columns=ASE_HEADER))
    split_by_sample(ase_df, "sample", args.samples, ASE_HEADER, args.outdir, "ase_hits")

    print("\nSplitting junction hits by sample (concatenated across tissues)...")
    per_tissue_dfs = []
    for tissue, jxn_file in zip(args.tissues, args.junction_files):
        if os.path.isfile(jxn_file):
            per_tissue_dfs.append(pd.read_csv(jxn_file, sep="\t", dtype=str))
        else:
            print(f"WARNING: No junction hits file found for tissue {tissue} at {jxn_file}")
    junction_df = (pd.concat(per_tissue_dfs, ignore_index=True)
                   if per_tissue_dfs else pd.DataFrame(columns=JUNCTION_HEADER))
    split_by_sample(junction_df, "sample", args.samples, JUNCTION_HEADER, args.outdir, "junction_hits")

    print("\nDone.")


if __name__ == "__main__":
    main()
