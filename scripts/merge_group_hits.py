#!/usr/bin/env python3

# scripts/merge_group_hits.py
# Group-level candidate-hits builder -- one job per (bed_id, sample_type)
# group. Replaces what used to be three separate pipeline stages:
#   1. splitting the group's already-merged variant/ASE/junction tables out
#      per sample (formerly scripts/split_group_hits_by_sample.py),
#   2. building each sample's ranked hit table (one job per sample, formerly
#      calling scripts/merge_hits.py's CLI once per sample),
#   3. concatenating every sample's table back into one group-level
#      all_hits.tsv (formerly rule _6E's awk concat).
# All three now happen in a single job/script, using merge_hits.py as a
# library for the per-sample build_hit_table logic -- see
# rules/6_merge_hits.smk's _6D1/_6D2 for why this replaced N per-sample
# split+merge jobs plus a separate group-level concat job.

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_hits


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build every sample's ranked candidate-hits table for a group, directly from the "
                    "group's already-merged variant/ASE/junction tables, and write one concatenated "
                    "all_hits.tsv -- replaces split_group_hits_by_sample.py + per-sample merge_hits.py + "
                    "the old awk concat.")
    parser.add_argument("--outfile", required=True, help="Path to write the group-level all_hits.tsv to.")
    parser.add_argument("--variant-tsv", required=True,
        help="Group-level merged variant hits (output of merge_and_filter_variants.py).")
    parser.add_argument("--ase-tsv", required=True,
        help="Group-level merged ASE hits (output of merge_and_filter_ase_results.py).")
    parser.add_argument("--tissues", nargs="*", default=[],
        help="Tissues present in this group (same order as --junction-files).")
    parser.add_argument("--junction-files", nargs="*", default=[],
        help="Per-tissue merged junction hit files (output of merge_and_filter_junction_results.py, "
             "one per --tissues entry). Concatenated across tissues before filtering by sample.")
    parser.add_argument("--cohort-junction-tsv", required=False, default=None,
        help="Group-level cohort-comparison junction outliers (rules/7_cohort_junction_analysis.smk's "
             "*_outliers_filtered.tsv). Raw column names differ from --junction-files (and vary by "
             "whether this group used the beta_binomial or modified_zscore method) -- normalized to "
             "merge_hits.load_cohort_junction_df's schema here. If omitted/missing, every sample's "
             "cohort_* columns are filled with '.' (merge_hits.build_hit_table's existing fallback).")
    parser.add_argument("--samples", nargs="+", required=True, help="Every sample name in this group.")
    parser.add_argument("--omim", required=False, default=None,
        help="Path to OMIM data. If omitted, phenotypes/inheritance_patterns are filled with '.'.")
    args = parser.parse_args()
    if len(args.tissues) != len(args.junction_files):
        parser.error("--tissues and --junction-files must have the same number of entries")
    return args


def normalize_cohort_junction_df(df):
    """rules/7_cohort_junction_analysis.smk's *_outliers_filtered.tsv has a
    much richer, differently-named column set than the GTEx-comparison
    junction files: six independent metrics per junction (junction_PSI,
    junction_PSI_approx, 5ss_IR_ratio, 3ss_IR_ratio, junction_full_IR_ratio,
    junction_IPA_ratio -- see identify_cohort_junction_outliers.py's
    _METRIC_EVENTS), and each metric's effect-size column's name/units
    depend on whether that group used the beta_binomial method
    (delta_{metric}, a ratio-scale value in [-1, 1]) or modified_zscore
    (modz_{metric}, an unbounded z-score). Normalize down to one delta
    column per metric -- whichever of delta_{metric}/modz_{metric} is
    present -- matching merge_hits.load_cohort_junction_df's expected schema."""
    out = pd.DataFrame()
    out['sample'] = df['sample']
    out['phasing'] = df['phasing']
    out['gene'] = df['gene']
    out['junction'] = df['junction']
    out['jxn_coverage'] = df['junction_coverage'] if 'junction_coverage' in df.columns else pd.NA

    _METRIC_TO_OUTCOL = {
        'junction_PSI':           'delta_PSI',
        'junction_PSI_approx':    'delta_PSI_approx',
        '5ss_IR_ratio':           'delta_5ss_IR',
        '3ss_IR_ratio':           'delta_3ss_IR',
        'junction_full_IR_ratio': 'delta_full_IR',
        'junction_IPA_ratio':     'delta_IPA',
    }
    for metric, outcol in _METRIC_TO_OUTCOL.items():
        if f'delta_{metric}' in df.columns:
            out[outcol] = df[f'delta_{metric}']
        elif f'modz_{metric}' in df.columns:
            out[outcol] = df[f'modz_{metric}']
        else:
            out[outcol] = pd.NA

    out['annotation'] = df['junction_type'] if 'junction_type' in df.columns else '.'
    out['event'] = df['event_type'] if 'event_type' in df.columns else '.'
    out['sample_count'] = df['n_sample_outlier_junction_PSI'] if 'n_sample_outlier_junction_PSI' in df.columns else pd.NA
    return out


_EMPTY_JUNCTION_COLS = ['sample', 'gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event']


def main():
    args = parse_args()

    print(f"\nBuilding group-level hits for {len(args.samples)} sample(s)...")

    variant_df = merge_hits.load_variant_df(args.variant_tsv)
    ase_df = merge_hits.load_ase_df(args.ase_tsv)

    print("Reading junction hits (concatenated across tissues)...")
    per_tissue_dfs = []
    for tissue, jxn_file in zip(args.tissues, args.junction_files):
        if os.path.isfile(jxn_file):
            per_tissue_dfs.append(merge_hits.load_junction_df(jxn_file))
        else:
            print(f"WARNING: No junction hits file found for tissue {tissue} at {jxn_file}")
    junction_df = (pd.concat(per_tissue_dfs, ignore_index=True) if per_tissue_dfs
                   else pd.DataFrame(columns=_EMPTY_JUNCTION_COLS))

    print("Reading cohort-comparison junction hits...")
    cohort_junction_df = None
    if args.cohort_junction_tsv and os.path.isfile(args.cohort_junction_tsv):
        _required_cols = {'sample', 'phasing', 'gene', 'junction'}
        raw_cohort_df = pd.read_csv(args.cohort_junction_tsv, sep='\t', dtype=str)
        if _required_cols.issubset(raw_cohort_df.columns) and not raw_cohort_df.empty:
            cohort_junction_df = normalize_cohort_junction_df(raw_cohort_df)
        # else: missing required columns, or a skipped group (rule 7 writes
        # a near-empty outliers_filtered.tsv for groups below its
        # min-samples threshold) -- fall through to None (merge_hits.
        # build_hit_table's existing "no cohort data" fallback) rather than
        # treating this as an error.

    omim_df = merge_hits.load_omim_df(args.omim) if args.omim else None

    hit_dfs = []
    for sample in args.samples:
        sample_variant_df = variant_df[variant_df['sample'].astype(str) == str(sample)]
        sample_ase_df = ase_df[ase_df['sample'].astype(str) == str(sample)]
        sample_junction_df = junction_df[junction_df['sample'].astype(str) == str(sample)]
        sample_cohort_junction_df = (
            cohort_junction_df[cohort_junction_df['sample'].astype(str) == str(sample)]
            if cohort_junction_df is not None else None
        )
        hit_df = merge_hits.build_hit_table(
            sample_variant_df, sample_ase_df, sample_junction_df, sample_cohort_junction_df,
            sample, omim_df,
        )
        hit_dfs.append(hit_df)
        print(f"  {sample}: {len(hit_df)} gene(s)")

    all_hits = pd.concat(hit_dfs, ignore_index=True) if hit_dfs else pd.DataFrame()

    outdir = os.path.dirname(args.outfile)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    all_hits.to_csv(args.outfile, sep='\t', index=False)
    print(f"\nSaved group-level hits ({len(all_hits)} total row(s)) to {args.outfile}")


if __name__ == "__main__":
    main()
