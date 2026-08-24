#! /usr/bin/env python3

# scripts/merge_hits.py
# Builds one sample's ranked candidate-hits table from that sample's
# variant/ASE/junction/cohort-junction rows. Exposes load_*()/build_hit_table()
# as a library (imported by scripts/merge_group_hits.py, which calls
# build_hit_table() once per sample in a group and concatenates the results)
# as well as a standalone single-sample CLI (below) for manual use/debugging.

import argparse
import os
import re
import pandas as pd
import numpy as np
import warnings

# Ignore FutureWarnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)

# Maps each raw delta source column to its output-column suffix. The
# GTEx-comparison junction_df only ever has 'delta_PSI' (one metric); the
# cohort-comparison cohort_junction_df has all six (see
# scripts/merge_group_hits.py's normalize_cohort_junction_df).
_DELTA_COL_TO_SUFFIX = {
    'delta_PSI':        'deltaPSI',
    'delta_PSI_approx': 'deltaPSIapprox',
    'delta_5ss_IR':     'delta5ssIR',
    'delta_3ss_IR':     'delta3ssIR',
    'delta_full_IR':    'deltaFullIR',
    'delta_IPA':        'deltaIPA',
}

# Which delta-metric suffixes exist per prefix -- '' (GTEx-comparison) only
# ever has deltaPSI; 'cohort_' has all six. max_deltas() (inside
# build_hit_table) takes the max magnitude across whichever of these are
# present for a given prefix: no per-metric weighting, just the single
# largest |delta| across all populated metrics.
_DELTA_SUFFIXES_BY_PREFIX = {
    '':        ('deltaPSI',),
    'cohort_': ('deltaPSI', 'deltaPSIapprox', 'delta5ssIR', 'delta3ssIR', 'deltaFullIR', 'deltaIPA'),
}


def parse_args():
    """ Parse command line arguments """

    parser = argparse.ArgumentParser(description='Merge hits from variant, ASE, and junction analyses for a single sample.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file')
    parser.add_argument('--sample-name', required=True, type=str, help='Name of the sample')
    parser.add_argument('--variant-hits', required=True, type=str, help='Path to variant hit results. At a minimum, file should contain columns: \
                                                                            sample, chrom, pos, ref, alt, GT, gnomAD_AF, CLNSIG, gene, CADD_PHRED, SpliceAI, num_callers, sample_count, \
                                                                            ANNOVAR_AAChange.refGene, ANNOVAR_GeneDetail.refGene')
    parser.add_argument('--ase-hits', required=True, type=str, help='Path to ASE hit results. At a minimum, file should contain columns: sample, gene, ratio')
    parser.add_argument('--junction-hits', required=True, type=str, help='Path to junction hit results. At a minimum, file should contain columns: \
                                                                            sample, gene, phasing, junction, delta_PSI, sample_count, annotation, event')
    parser.add_argument('--cohort-junction-hits', required=False, default=None, type=str, help='Path to cohort-comparison junction hit results (rules/7_cohort_junction_analysis.smk), \
                                                                            same schema as --junction-hits but the sample compared against the rest of its cohort rather than GTEx. \
                                                                            If omitted, the cohort_* junction columns are filled with "." rather than annotated.')
    parser.add_argument('--omim', required=False, default=None, type=str, help='Path to OMIM data. At a minimum, file should contain columns: approved_gene_symbol, phenotypes, inheritance_patterns. \
                                                                                  If omitted, the phenotypes/inheritance_patterns columns are filled with "." rather than annotated.')
    
    return parser.parse_args()


def load_variant_df(path):
    """Read + preprocess a variant-hits tsv (single-sample or group-level --
    same schema either way, just more rows for the latter)."""
    variant_df = pd.read_csv(path, sep='\t', usecols=[
        'sample', 'chrom', 'pos', 'ref', 'alt', 'GT', 'gnomAD_AF', 'CLNSIG', 'gene', 'CADD_PHRED', 'SpliceAI',
        'num_callers', 'sample_count', 'ANNOVAR_AAChange.refGene', 'ANNOVAR_GeneDetail.refGene',
    ]).drop_duplicates()
    variant_df = variant_df.rename(columns={'sample_count': 'variant_nsamples'})
    variant_df = variant_df[variant_df['gene'] != '.']
    # 'gene' lists every BED-panel gene a variant overlaps (comma-separated
    # if more than one) -- explode so a variant overlapping multiple genes
    # is grouped into each of those genes' hit rows individually.
    variant_df['gene'] = variant_df['gene'].str.split(',')
    variant_df = variant_df.explode('gene')
    variant_df['GT'] = variant_df['GT'].str.replace('|', '/', regex=False).str.replace('1/0', '0/1', regex=False)
    variant_df['variant_ID'] = variant_df.apply(lambda x: f"{x.chrom}-{x.pos}-{x.ref}-{x.alt}", axis=1).drop_duplicates()
    # ANNOVAR only ever populates one of these two per variant (AAChange.refGene
    # for exonic variants with a codon change to report, GeneDetail.refGene for
    # everything else it has position detail for -- splicing/UTR/intronic/etc.)
    # -- coalesce into a single human-readable consequence column rather than
    # carrying two mostly-empty columns through to the final hits table.
    aachange = variant_df['ANNOVAR_AAChange.refGene']
    genedetail = variant_df['ANNOVAR_GeneDetail.refGene']
    variant_df['variant_consequence'] = aachange.where(aachange.notna() & (aachange != '.'), genedetail)
    variant_df['variant_consequence'] = variant_df['variant_consequence'].fillna('.')
    return variant_df


def load_ase_df(path):
    """Read + preprocess an ASE-hits tsv (single-sample or group-level)."""
    ase_df = pd.read_csv(path, sep='\t', usecols=['sample', 'gene', 'ratio', 'sample_count']).drop_duplicates()
    return ase_df.rename(columns={'ratio': 'ASE_ratio', 'sample_count': 'ASE_nsamples'})


def load_junction_df(path):
    """Read a GTEx-comparison junction-hits tsv (single-sample or group-level) -- one metric (delta_PSI)."""
    return pd.read_csv(path, sep='\t', usecols=[
        'sample', 'gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event',
    ]).drop_duplicates()


def load_cohort_junction_df(path):
    """Read a cohort-comparison junction-hits tsv (already normalized to one
    delta column per metric by scripts/merge_group_hits.py's
    normalize_cohort_junction_df -- six metrics, unlike load_junction_df's one)."""
    return pd.read_csv(path, sep='\t', usecols=[
        'sample', 'gene', 'phasing', 'junction', 'jxn_coverage',
        'delta_PSI', 'delta_PSI_approx', 'delta_5ss_IR', 'delta_3ss_IR', 'delta_full_IR', 'delta_IPA',
        'sample_count', 'annotation', 'event',
    ]).drop_duplicates()


def load_omim_df(path):
    omim_df = pd.read_csv(path, sep='\t', usecols=['approved_gene_symbol', 'phenotypes', 'inheritance_patterns'])
    return omim_df.rename(columns={'approved_gene_symbol': 'gene'})


def build_phased_junction_df(df, prefix, delta_cols=('delta_PSI',)):
    """Aggregate a junction_df-shaped table (gene, phasing, junction,
    jxn_coverage, one-or-more delta columns, sample_count, annotation,
    event) into one gene-level row per phasing tier, with columns named
    '{prefix}bulk_jxns', '{prefix}bulk_jxn_coverage', one
    '{prefix}bulk_{suffix}' per entry in delta_cols, etc. Used for both
    the GTEx-comparison junction_df (prefix='', delta_cols=('delta_PSI',))
    and the cohort-comparison cohort_junction_df (prefix='cohort_',
    delta_cols=all six metrics)."""
    tiers = {}
    for phasing, sep in (('bulk', ';'), ('hap1', ','), ('hap2', ',')):
        agg_kwargs = {
            prefix + phasing + '_jxns':         ('junction', lambda x, sep=sep: sep.join(map(str, x))),
            prefix + phasing + '_jxn_coverage': ('jxn_coverage', lambda x, sep=sep: sep.join(map(str, x))),
        }
        for col in delta_cols:
            suffix = _DELTA_COL_TO_SUFFIX[col]
            agg_kwargs[prefix + phasing + '_' + suffix] = (col, lambda x, sep=sep: sep.join(map(str, x)))
        agg_kwargs[prefix + phasing + '_jxn_annotation'] = ('annotation', lambda x, sep=sep: sep.join(map(str, x)))
        agg_kwargs[prefix + phasing + '_jxn_event']      = ('event', lambda x, sep=sep: sep.join(map(str, x)))
        agg_kwargs[prefix + phasing + '_jxn_nsamples']   = ('sample_count', lambda x, sep=sep: sep.join(map(str, x)))
        tiers[phasing] = (
            df[df['phasing'] == phasing]
                .sort_values('junction')
                .groupby('gene')
                .agg(**agg_kwargs)
                .reset_index()
        )
    merged = pd.merge(tiers['bulk'], tiers['hap1'], on='gene', how='outer')
    merged = pd.merge(merged, tiers['hap2'], on='gene', how='outer')
    return merged.drop_duplicates()


def build_hit_table(variant_df, ase_df, junction_df, cohort_junction_df, sample_name, omim_df=None):
    """Build one sample's ranked candidate-hits table. All four input
    DataFrames are assumed already loaded (via load_*() above) and already
    filtered down to this sample only -- this function itself is agnostic
    to whether they came from single-sample files or were filtered out of
    group-level tables. cohort_junction_df/omim_df may be None (same
    fallback behavior as omitting --cohort-junction-hits/--omim on the CLI:
    cohort_* columns filled with '.', phenotypes/inheritance_patterns
    filled with '.')."""

    # Create modified DataFrames
    # mod_variant_df should contain genes, variant_pos, variant_GT, variant_CLNSIG, variant_nsamples
    mod_variant_df = (
        variant_df
            .sort_values(['gene', 'variant_ID'])
            .groupby(['gene', 'variant_ID'], sort=False)
            .agg(
                GT=('GT', lambda x: ','.join(x)),
                gnomAD_AF=('gnomAD_AF', lambda x: ','.join(x.dropna().astype(str))),
                CLNSIG=('CLNSIG', lambda x: ','.join(x.dropna().astype(str))),
                CADD_PHRED=('CADD_PHRED', lambda x: ','.join(x.dropna().astype(str))),
                SpliceAI=('SpliceAI', lambda x: ','.join(x.dropna().astype(str))),
                variant_consequence=('variant_consequence', lambda x: ','.join(x.dropna().astype(str))),
                num_callers=('num_callers', lambda x: ','.join(x.astype(str))),
                variant_nsamples=('variant_nsamples', lambda x: ','.join(x.astype(str))),
            )
            .reset_index()
            .groupby('gene', sort=False)
            .agg(
                variant_ID=('variant_ID', ';'.join),
                variant_GT=('GT', ';'.join),
                variant_gnomAD_AF=('gnomAD_AF', ';'.join),
                variant_CLNSIG=('CLNSIG', ';'.join),
                variant_CADD_PHRED=('CADD_PHRED', ';'.join),
                variant_SpliceAI=('SpliceAI', ';'.join),
                variant_consequence=('variant_consequence', ';'.join),
                variant_num_callers=('num_callers', ';'.join),
                variant_nsamples=('variant_nsamples', ';'.join),
            )
            .reset_index()
    )
    # mod_junction_df should contain genes, bulk_jxns, bulk_jxn_coverage, bulk_deltaPSI, bulk_jxn_nsamples, 
    # hap1_jxns, hap1_jxn_coverage, hap1_deltaPSI, hap1_jxn_nsamples, hap2_jxns, hap2_jxn_coverage, hap2_delta_PSI, hap2_jxn_nsamples
    mod_junction_df = build_phased_junction_df(junction_df, '')
    if cohort_junction_df is not None:
        mod_cohort_junction_df = build_phased_junction_df(
            cohort_junction_df, 'cohort_',
            delta_cols=('delta_PSI', 'delta_PSI_approx', 'delta_5ss_IR', 'delta_3ss_IR', 'delta_full_IR', 'delta_IPA'),
        )
    else:
        mod_cohort_junction_df = None

    # Merge hits
    hit_df = pd.merge(mod_variant_df, ase_df, on='gene', how='outer')
    hit_df = pd.merge(hit_df, mod_junction_df, on='gene', how='outer')
    if mod_cohort_junction_df is not None:
        hit_df = pd.merge(hit_df, mod_cohort_junction_df, on='gene', how='outer')
    else:
        # No cohort-comparison junction data provided -- keep the same
        # output schema, just unannotated, rather than dropping these
        # columns entirely (mirrors the phenotypes/inheritance_patterns
        # fallback below when omim_df is omitted).
        for phasing in ('bulk', 'hap1', 'hap2'):
            for suffix in ('jxns', 'jxn_coverage', 'deltaPSI', 'deltaPSIapprox', 'delta5ssIR',
                           'delta3ssIR', 'deltaFullIR', 'deltaIPA', 'jxn_annotation', 'jxn_event', 'jxn_nsamples'):
                hit_df['cohort_' + phasing + '_' + suffix] = '.'
    if omim_df is not None:
        hit_df = pd.merge(hit_df, omim_df, on='gene', how='left')
    else:
        # No OMIM data provided -- keep the same output schema, just
        # unannotated, rather than dropping these columns entirely.
        hit_df['phenotypes'] = '.'
        hit_df['inheritance_patterns'] = '.'
    hit_df = hit_df.drop_duplicates()

    # Create boolean columns to indicate whether there is a candidate variant or allele-specific expression
    hit_df["variant"] = hit_df["variant_ID"].notna()
    hit_df["pathogenic_variant"] = (
        hit_df["variant_CLNSIG"]
            .astype(str)
            .str.contains(r"Pathogenic|Likely_pathogenic", regex=True)
    )
    hit_df["ASE"] = pd.to_numeric(hit_df["ASE_ratio"], errors="coerce").notna()

    # Add a column to indicate whether a gene has strong, moderate, weak, splicing dysregulation

    def max_deltas(row, prefix=''):
        """ Extract the max |delta| across all populated delta metrics for a gene, per phasing tier """
        def parse_vals(s):
            if not pd.notna(s):
                return []
            vals = []
            for v in re.split('[,;]', s):
                v = v.strip()
                if not v:
                    continue
                try:
                    vals.append(float(v))
                except ValueError:
                    # Non-numeric sentinel from rules/7_cohort_junction_analysis.smk's
                    # identify_cohort_junction_outliers.py -- "low_n" (too few
                    # cohort samples with good coverage), "error" (beta_binomial
                    # fit failed), or "no_variance" (modified_zscore: zero
                    # variance in the reference distribution). These mean "not
                    # statistically testable", not a numeric value -- skip
                    # rather than crash, same as a blank/missing entry.
                    continue
            return vals

        def max_for_phasing(phasing):
            all_vals = []
            for suffix in _DELTA_SUFFIXES_BY_PREFIX[prefix]:
                all_vals.extend(parse_vals(row.get(prefix + phasing + '_' + suffix)))
            return max(abs(v) for v in all_vals) if all_vals else np.nan

        max_bulk = max_for_phasing('bulk')
        max_hap1 = max_for_phasing('hap1')
        max_hap2 = max_for_phasing('hap2')

        return max_bulk, max_hap1, max_hap2

    def inspect_row(row, prefix=''):
        """ Determine if a gene has splicing dysregulation """
        max_bulk, max_hap1, max_hap2 = max_deltas(row, prefix)
        dominant = any(x in row['inheritance_patterns'] for x in ['AD', 'XLD']) if pd.notna(row.get('inheritance_patterns')) else False
        if (
            (max_bulk >= 0.5)
            or (max_hap1 >= 0.5 and max_hap2 >= 0.5)
            or (dominant and (max_bulk >= 0.2 or max_hap1 >= 0.5 or max_hap2 >= 0.5))
        ):
            return "Strong"
        elif (max_bulk >= 0.2 or max_hap1 >= 0.5 or max_hap2 >= 0.5):
            return "Moderate"
        elif (max_bulk > 0 or max_hap1 > 0 or max_hap2 > 0):
            return "Weak"
        else:
            return "None"
    hit_df['outlier_junction'] = hit_df.apply(inspect_row, axis=1)
    if mod_cohort_junction_df is not None:
        # NOTE: cohort_bulk_deltaPSI (etc.) is delta_junction_PSI for groups
        # analyzed with the beta_binomial method (a PSI-scale value in
        # [-1, 1], same units the 0.2/0.5 thresholds below assume) but
        # modz_junction_PSI for groups analyzed with modified_zscore (an
        # unbounded z-score) -- rules/7_cohort_junction_analysis.smk picks
        # the method per group based on its sample count. Strong/Moderate/
        # Weak is therefore not on a consistent scale across groups that
        # used different methods; treat cross-group comparisons of this
        # column with that in mind.
        hit_df['cohort_outlier_junction'] = hit_df.apply(lambda row: inspect_row(row, prefix='cohort_'), axis=1)
    else:
        hit_df['cohort_outlier_junction'] = '.'

    # Fill missing values
    # Cast to object before fillna: on pandas >=3.0, filling a still-numeric
    # column (e.g. ASE_ratio, which is never string-joined the way the
    # variant/junction columns are) with the string "." raises
    # LossySetitemError instead of silently upcasting like older pandas did.
    hit_df = hit_df.astype(object)
    hit_df.fillna(".", inplace=True)

    # Rank hits by a points-based score:
    #   - variant:           2 points if pathogenic_variant, else 1 point if
    #                        variant, else 0 (these are NOT additive --
    #                        pathogenic_variant implies variant, but the
    #                        category caps at 2, it isn't 1 + 2 = 3)
    #   - ASE:               2 points if True, else 0
    #   - outlier_junction:  2 points if Strong/Moderate, 1 if Weak, else 0
    #   - cohort_outlier_junction: same 2/1/0 scale as outlier_junction
    def _variant_points(row):
        if row['pathogenic_variant']:
            return 2
        if row['variant']:
            return 1
        return 0

    def _junction_points(tier):
        if tier in ('Strong', 'Moderate'):
            return 2
        if tier == 'Weak':
            return 1
        return 0

    hit_df['score'] = (
        hit_df.apply(_variant_points, axis=1)
        + hit_df['ASE'].apply(lambda x: 2 if x is True else 0)
        + hit_df['outlier_junction'].apply(_junction_points)
        + hit_df['cohort_outlier_junction'].apply(_junction_points)
    )

    # Tiebreak chain (applied only when two genes have the same score),
    # in priority order:
    #   1. pathogenic_variant (True first)
    #   2. ASE (True first)
    #   3. Strong outlier_junction OR Strong cohort_outlier_junction (True first)
    #   4. max magnitude of bulk deltaPSI (outlier_junction) or bulk
    #      cohort deltaPSI/deltaPSIapprox/delta5ssIR/delta3ssIR/deltaFullIR/
    #      deltaIPA (cohort_outlier_junction) -- whichever is larger,
    #      descending
    #   5. max CADD score across the gene's variant(s), descending
    #   6. gene name, alphabetically
    def _cadd_max(row):
        vals = []
        for v in re.split('[,;]', str(row.get('variant_CADD_PHRED'))):
            v = v.strip()
            if not v:
                continue
            try:
                vals.append(float(v))
            except ValueError:
                continue
        return max(vals) if vals else -1

    hit_df['_tb_strong'] = (
        (hit_df['outlier_junction'] == 'Strong') | (hit_df['cohort_outlier_junction'] == 'Strong')
    )
    def _max_bulk_delta(row):
        a = max_deltas(row, '')[0]
        b = max_deltas(row, 'cohort_')[0]
        vals = [v for v in (a, b) if pd.notna(v)]
        return max(vals) if vals else -1

    hit_df['_tb_max_bulk_delta'] = hit_df.apply(_max_bulk_delta, axis=1)
    hit_df['_tb_max_cadd'] = hit_df.apply(_cadd_max, axis=1)

    hit_df.sort_values(
        by=['score', 'pathogenic_variant', 'ASE', '_tb_strong', '_tb_max_bulk_delta', '_tb_max_cadd', 'gene'],
        ascending=[False, False, False, False, False, False, True],
        na_position='last',
        inplace=True
    )
    hit_df.drop(columns=['_tb_strong', '_tb_max_bulk_delta', '_tb_max_cadd'], inplace=True)
    hit_df['ranking'] = np.arange(1, len(hit_df) + 1)

    # Reorder
    hit_df['sample'] = sample_name
    hit_df = hit_df[[
        'sample', 'gene', 'phenotypes', 'inheritance_patterns', 'ranking', 'score', 'variant', 'pathogenic_variant', 'ASE', 'outlier_junction', 'cohort_outlier_junction',
        'variant_ID', 'variant_GT', 'variant_gnomAD_AF',  'variant_CLNSIG', 'variant_CADD_PHRED', 'variant_SpliceAI', 'variant_consequence', 'variant_num_callers', 'variant_nsamples',
        'ASE_ratio', 'ASE_nsamples', 
        'bulk_jxns', 'bulk_jxn_coverage', 'bulk_deltaPSI', 'bulk_jxn_annotation', 'bulk_jxn_event', 'bulk_jxn_nsamples',
        'hap1_jxns', 'hap1_jxn_coverage', 'hap1_deltaPSI', 'hap1_jxn_annotation', 'hap1_jxn_event', 'hap1_jxn_nsamples',
        'hap2_jxns', 'hap2_jxn_coverage', 'hap2_deltaPSI', 'hap2_jxn_annotation', 'hap2_jxn_event', 'hap2_jxn_nsamples',
        'cohort_bulk_jxns', 'cohort_bulk_jxn_coverage',
        'cohort_bulk_deltaPSI', 'cohort_bulk_deltaPSIapprox', 'cohort_bulk_delta5ssIR', 'cohort_bulk_delta3ssIR', 'cohort_bulk_deltaFullIR', 'cohort_bulk_deltaIPA',
        'cohort_bulk_jxn_annotation', 'cohort_bulk_jxn_event', 'cohort_bulk_jxn_nsamples',
        'cohort_hap1_jxns', 'cohort_hap1_jxn_coverage',
        'cohort_hap1_deltaPSI', 'cohort_hap1_deltaPSIapprox', 'cohort_hap1_delta5ssIR', 'cohort_hap1_delta3ssIR', 'cohort_hap1_deltaFullIR', 'cohort_hap1_deltaIPA',
        'cohort_hap1_jxn_annotation', 'cohort_hap1_jxn_event', 'cohort_hap1_jxn_nsamples',
        'cohort_hap2_jxns', 'cohort_hap2_jxn_coverage',
        'cohort_hap2_deltaPSI', 'cohort_hap2_deltaPSIapprox', 'cohort_hap2_delta5ssIR', 'cohort_hap2_delta3ssIR', 'cohort_hap2_deltaFullIR', 'cohort_hap2_deltaIPA',
        'cohort_hap2_jxn_annotation', 'cohort_hap2_jxn_event', 'cohort_hap2_jxn_nsamples'
    ]]
    return hit_df


def main():
    """ Main function -- standalone single-sample CLI. """

    # Parse command line arguments
    args = parse_args()
    for attr, value in vars(args).items():
        if value == "None":
            setattr(args, attr, None)
    print(f"\nMerging variant, ASE, and outlier junction hits for sample: {args.sample_name}")

    variant_df = load_variant_df(args.variant_hits)
    ase_df = load_ase_df(args.ase_hits)
    junction_df = load_junction_df(args.junction_hits)
    cohort_junction_df = load_cohort_junction_df(args.cohort_junction_hits) if args.cohort_junction_hits else None
    omim_df = load_omim_df(args.omim) if args.omim else None

    hit_df = build_hit_table(variant_df, ase_df, junction_df, cohort_junction_df, args.sample_name, omim_df)

    # Make output directory if it doesn't exist
    if args.outfile:
        outdir = os.path.dirname(args.outfile)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir)

    # Save hits
    hit_df.to_csv(args.outfile, sep='\t', index=False)
    print(f"Saved merged hits to {args.outfile}")

if __name__ == "__main__":
    main()
