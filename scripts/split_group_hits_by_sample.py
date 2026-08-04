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
                   'ANNOVAR_AAChange.refGene', 'ANNOVAR_GeneDetail.refGene',
                   'num_callers', 'sample_count']
ASE_HEADER = ['sample', 'gene', 'ratio', 'sample_count']
JUNCTION_HEADER = ['sample', 'phasing', 'gene', 'junction', 'jxn_coverage', 'delta_PSI',
                    'annotation', 'event', 'sample_count']


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
    parser.add_argument("--cohort-junction-tsv", required=False, default=None,
        help="Group-level cohort-comparison junction outliers (rules/7_cohort_junction_analysis.smk's "
             "*_outliers_filtered.tsv). Raw column names differ from --junction-files (and vary by "
             "whether that group used the beta_binomial or modified_zscore method) -- normalized to "
             "the same schema as --junction-files here so merge_hits.py can treat them identically.")
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


def normalize_cohort_junction_df(df):
    """rules/7_cohort_junction_analysis.smk's *_outliers_filtered.tsv has a
    much richer, differently-named column set than the GTEx-comparison
    junction files (and its effect-size column's name/units depend on
    whether that group used the beta_binomial or modified_zscore method --
    delta_junction_PSI (a PSI-scale value in [-1, 1]) for the former,
    modz_junction_PSI (an unbounded z-score) for the latter). Normalize
    down to the same schema as JUNCTION_HEADER so merge_hits.py can
    aggregate both with identical logic."""
    out = pd.DataFrame()
    out['sample'] = df['sample']
    out['phasing'] = df['phasing']
    out['gene'] = df['gene']
    out['junction'] = df['junction']
    out['jxn_coverage'] = df['junction_coverage'] if 'junction_coverage' in df.columns else pd.NA
    if 'delta_junction_PSI' in df.columns:
        out['delta_PSI'] = df['delta_junction_PSI']
    elif 'modz_junction_PSI' in df.columns:
        out['delta_PSI'] = df['modz_junction_PSI']
    else:
        out['delta_PSI'] = pd.NA
    out['annotation'] = df['junction_type'] if 'junction_type' in df.columns else '.'
    out['event'] = df['event_type'] if 'event_type' in df.columns else '.'
    out['sample_count'] = df['n_sample_outlier_junction_PSI'] if 'n_sample_outlier_junction_PSI' in df.columns else pd.NA
    return out


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

    print("\nSplitting cohort-comparison junction hits by sample...")
    if args.cohort_junction_tsv:
        _required_cols = {'sample', 'phasing', 'gene', 'junction'}
        raw_cohort_df = (pd.read_csv(args.cohort_junction_tsv, sep="\t", dtype=str)
                          if os.path.isfile(args.cohort_junction_tsv) else pd.DataFrame())
        if _required_cols.issubset(raw_cohort_df.columns) and not raw_cohort_df.empty:
            cohort_junction_df = normalize_cohort_junction_df(raw_cohort_df)
        else:
            # Missing file, empty file, or a skipped group (rule 7 writes a
            # near-empty outliers_filtered.tsv for groups below its
            # min-samples threshold) -- fall through to header-only stubs
            # below rather than treating this as an error.
            cohort_junction_df = pd.DataFrame(columns=JUNCTION_HEADER)
    else:
        cohort_junction_df = pd.DataFrame(columns=JUNCTION_HEADER)
    split_by_sample(cohort_junction_df, "sample", args.samples, JUNCTION_HEADER, args.outdir, "cohort_junction_hits")

    print("\nDone.")


if __name__ == "__main__":
    main()
