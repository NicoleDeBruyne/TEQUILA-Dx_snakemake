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
        help="Path to OMIM data. If omitted, phenotypes/inheritance_patterns/haploinsufficient are filled with '.'/False.")
    parser.add_argument("--gene-expression-matrix", required=False, default=None,
        help="Group-level targeted-panel CPTM matrix (rule _10D's <outprefix>_matrix.tsv, from "
             "quantify_gene_by_assignment.py) -- one row per gene, one column per sample in this "
             "(bed_id, sample_type) group. If omitted/missing, relative_gene_expression/"
             "cohort_relative_gene_expression/n_cohort are filled with '.' (merge_hits.build_hit_table's "
             "existing fallback convention for optional inputs).")
    parser.add_argument("--debug-sample", required=False, default=None,
        help="If set, print diagnostic detail (to stderr) for this one sample: the gene list surviving "
             "the per-sample variant_df slice, and the gene list in build_hit_table's output -- useful "
             "for tracking down a gene that unexpectedly disappears for one specific sample. No effect "
             "on output files. Off by default.")
    parser.add_argument("--debug-gene", required=False, default=None,
        help="Used together with --debug-sample: also print the matching variant_df row(s) for this "
             "gene (if any) within that sample's slice, and whether the gene made it into "
             "build_hit_table's output. Ignored if --debug-sample isn't also set.")
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
    present -- matching merge_hits.load_cohort_junction_df's expected schema.

    A row exists in outliers_filtered.tsv if ANY of the six metrics was a
    significant outlier for that junction, but the row still carries a raw
    delta_{metric}/modz_{metric} value for every metric regardless of
    whether that particular metric was the one that triggered -- e.g. a
    junction flagged only via junction_full_IR_ratio (event_type="full_IR")
    still has a real-looking delta_junction_PSI value sitting right next to
    it, even though junction_PSI was never actually significant here.
    event_type (identify_cohort_junction_outliers.py's sig_df construction)
    is comma-joined across every metric that DID trigger for that junction,
    so a metric's delta value is only kept if event_type contains one of
    that metric's own event strings; otherwise it's blanked to '' (empty
    string, not NaN/pd.NA -- build_phased_junction_df below joins these
    into a single semicolon/comma-separated string per gene, and pd.NA
    stringifies to the literal text '<NA>' rather than disappearing, which
    would leak into the final merged_all_hits.tsv instead of being caught
    by merge_hits.build_hit_table's trailing fillna('.'). An empty string
    join segment is already handled correctly by max_deltas' parse_vals,
    same as a genuinely blank/missing entry)."""
    out = pd.DataFrame()
    out['sample'] = df['sample']
    out['phasing'] = df['phasing']
    out['gene'] = df['gene']
    out['junction'] = df['junction']
    out['jxn_coverage'] = df['junction_coverage'] if 'junction_coverage' in df.columns else ''

    # Mirrors identify_cohort_junction_outliers.py's _METRIC_EVENTS exactly
    # -- which event-type strings each metric can produce.
    _METRIC_EVENTS = {
        'junction_PSI_approx':    {'alt_5ss_approx', 'alt_3ss_approx', 'exon_skipping_approx', 'exon_inclusion_approx'},
        'junction_PSI':           {'alt_5ss', 'alt_3ss', 'exon_skipping', 'exon_inclusion'},
        '5ss_IR_ratio':           {'5ss_IR'},
        '3ss_IR_ratio':           {'3ss_IR'},
        'junction_full_IR_ratio': {'full_IR'},
        'junction_IPA_ratio':     {'IPA'},
    }
    _METRIC_TO_OUTCOL = {
        'junction_PSI':           'delta_PSI',
        'junction_PSI_approx':    'delta_PSI_approx',
        '5ss_IR_ratio':           'delta_5ss_IR',
        '3ss_IR_ratio':           'delta_3ss_IR',
        'junction_full_IR_ratio': 'delta_full_IR',
        'junction_IPA_ratio':     'delta_IPA',
    }

    event_type = df['event_type'] if 'event_type' in df.columns else pd.Series([''] * len(df), index=df.index)
    fired_events = event_type.fillna('').apply(lambda s: set(str(s).split(',')) if s else set())

    for metric, outcol in _METRIC_TO_OUTCOL.items():
        if f'delta_{metric}' in df.columns:
            raw = df[f'delta_{metric}']
        elif f'modz_{metric}' in df.columns:
            raw = df[f'modz_{metric}']
        else:
            out[outcol] = ''
            continue
        is_outlier_for_metric = fired_events.apply(lambda evs, m=metric: bool(evs & _METRIC_EVENTS[m]))
        out[outcol] = raw.where(is_outlier_for_metric, other='')

    out['annotation'] = df['junction_type'] if 'junction_type' in df.columns else '.'
    out['event'] = df['event_type'] if 'event_type' in df.columns else '.'
    out['sample_count'] = df['n_sample_outlier_junction_PSI'] if 'n_sample_outlier_junction_PSI' in df.columns else ''
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

    gene_expression_df = None
    if args.gene_expression_matrix and os.path.isfile(args.gene_expression_matrix):
        gene_expression_df = merge_hits.load_gene_expression_df(args.gene_expression_matrix)
    # else: falls through to merge_hits.build_hit_table's existing "no
    # gene expression data" fallback (same pattern as omim_df/
    # cohort_junction_tsv being optional above).

    hit_dfs = []
    for sample in args.samples:
        sample_variant_df = variant_df[variant_df['sample'].astype(str) == str(sample)]
        sample_ase_df = ase_df[ase_df['sample'].astype(str) == str(sample)]
        sample_junction_df = junction_df[junction_df['sample'].astype(str) == str(sample)]
        sample_cohort_junction_df = (
            cohort_junction_df[cohort_junction_df['sample'].astype(str) == str(sample)]
            if cohort_junction_df is not None else None
        )

        if args.debug_sample and str(sample) == args.debug_sample:
            genes_in_slice = sorted(sample_variant_df['gene'].dropna().unique().tolist())
            print(f"  [DEBUG] {sample}: {len(sample_variant_df)} variant_df row(s) after per-sample slice, "
                  f"genes: {genes_in_slice}", file=sys.stderr)
            if args.debug_gene:
                match = sample_variant_df[sample_variant_df['gene'] == args.debug_gene]
                print(f"  [DEBUG] {sample}: {len(match)} row(s) for gene=={args.debug_gene!r} in the slice:",
                      file=sys.stderr)
                if not match.empty:
                    print(match.to_string(), file=sys.stderr)

        hit_df = merge_hits.build_hit_table(
            sample_variant_df, sample_ase_df, sample_junction_df, sample_cohort_junction_df,
            sample, omim_df, gene_expression_df,
        )

        if args.debug_sample and str(sample) == args.debug_sample:
            genes_in_output = sorted(hit_df['gene'].dropna().unique().tolist())
            print(f"  [DEBUG] {sample}: {len(hit_df)} row(s) in build_hit_table's output, "
                  f"genes: {genes_in_output}", file=sys.stderr)
            if args.debug_gene:
                present = args.debug_gene in genes_in_output
                print(f"  [DEBUG] {sample}: gene=={args.debug_gene!r} present in output: {present}", file=sys.stderr)

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