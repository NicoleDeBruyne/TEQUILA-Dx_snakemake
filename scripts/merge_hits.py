#! /usr/bin/env python3

import argparse
import warnings
import subprocess
import os
import re
import pandas as pd
import numpy as np
import warnings

# Ignore FutureWarnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)

def parse_args():
    """ Parse command line arguments """

    parser = argparse.ArgumentParser(description='Merge hits from variant, ASE, and junction analyses for a single sample.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file')
    parser.add_argument('--sample-name', required=True, type=str, help='Name of the sample')
    parser.add_argument('--variant-hits', required=True, type=str, help='Path to variant hit results. At a minimum, file should contain columns: \
                                                                            chrom, pos, ref, alt, GT, gnomAD_AF, CLNSIG, gene, CADD_PHRED, SpliceAI, num_callers, sample_count, \
                                                                            ANNOVAR_AAChange.refGene, ANNOVAR_GeneDetail.refGene')
    parser.add_argument('--ase-hits', required=True, type=str, help='Path to ASE hit results. At a minimum, file should contain columns: gene, ratio')
    parser.add_argument('--junction-hits', required=True, type=str, help='Path to junction hit results. At a minimum, file should contain columns: \
                                                                            gene, phasing, junction, delta_PSI, sample_count, annotation, event')
    parser.add_argument('--cohort-junction-hits', required=False, default=None, type=str, help='Path to cohort-comparison junction hit results (rules/7_cohort_junction_analysis.smk), \
                                                                            same schema as --junction-hits but the sample compared against the rest of its cohort rather than GTEx. \
                                                                            If omitted, the cohort_* junction columns are filled with "." rather than annotated.')
    parser.add_argument('--omim', required=False, default=None, type=str, help='Path to OMIM data. At a minimum, file should contain columns: approved_gene_symbol, phenotypes, inheritance_patterns. \
                                                                                  If omitted, the phenotypes/inheritance_patterns columns are filled with "." rather than annotated.')
    
    return parser.parse_args()

def main():
    """ Main function """

    # Parse command line arguments
    args = parse_args()
    for attr, value in vars(args).items():
        if value == "None":
            setattr(args, attr, None)
    print(f"\nMerging variant, ASE, and outlier junction hits for sample: {args.sample_name}")

    # Read information from input files
    variant_df = pd.read_csv(args.variant_hits, sep='\t', usecols=['chrom', 'pos', 'ref', 'alt', 'GT', 'gnomAD_AF', 'CLNSIG', 'gene', 'CADD_PHRED', 'SpliceAI', 'num_callers', 'sample_count',
                                                                     'ANNOVAR_AAChange.refGene', 'ANNOVAR_GeneDetail.refGene']).drop_duplicates()
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
    ase_df = pd.read_csv(args.ase_hits, sep='\t', usecols=['gene', 'ratio', 'sample_count']).drop_duplicates()
    ase_df = ase_df.rename(columns={'ratio': 'ASE_ratio', 'sample_count': 'ASE_nsamples'})
    junction_df = pd.read_csv(args.junction_hits, sep='\t', usecols=['gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event']).drop_duplicates()
    if args.cohort_junction_hits:
        # Normalized to the same schema as --junction-hits by
        # split_group_hits_by_sample.py (which coalesces the beta_binomial-
        # vs modified_zscore-method column differences from rule 7's raw
        # output into these same names), so it can be aggregated with the
        # exact same logic below.
        cohort_junction_df = pd.read_csv(args.cohort_junction_hits, sep='\t', usecols=['gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event']).drop_duplicates()
    else:
        cohort_junction_df = None
    if args.omim:
        omim_df = pd.read_csv(args.omim, sep='\t', usecols=['approved_gene_symbol', 'phenotypes', 'inheritance_patterns'])
        omim_df = omim_df.rename(columns={'approved_gene_symbol': 'gene'})
    else:
        omim_df = None

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
    def build_phased_junction_df(df, prefix):
        """Aggregate a junction_df-shaped table (gene, phasing, junction,
        jxn_coverage, delta_PSI, sample_count, annotation, event) into one
        gene-level row per phasing tier, with columns named
        '{prefix}bulk_jxns', '{prefix}bulk_jxn_coverage', etc. Used for both
        the GTEx-comparison junction_df (prefix='') and the cohort-comparison
        cohort_junction_df (prefix='cohort_')."""
        tiers = {}
        for phasing, sep in (('bulk', ';'), ('hap1', ','), ('hap2', ',')):
            tiers[phasing] = (
                df[df['phasing'] == phasing]
                    .sort_values('junction')
                    .groupby('gene')
                    .agg(**{
                        prefix + phasing + '_jxns':          ('junction', lambda x, sep=sep: sep.join(map(str, x))),
                        prefix + phasing + '_jxn_coverage':  ('jxn_coverage', lambda x, sep=sep: sep.join(map(str, x))),
                        prefix + phasing + '_deltaPSI':      ('delta_PSI', lambda x, sep=sep: sep.join(map(str, x))),
                        prefix + phasing + '_jxn_annotation': ('annotation', lambda x, sep=sep: sep.join(map(str, x))),
                        prefix + phasing + '_jxn_event':     ('event', lambda x, sep=sep: sep.join(map(str, x))),
                        prefix + phasing + '_jxn_nsamples':  ('sample_count', lambda x, sep=sep: sep.join(map(str, x))),
                    })
                    .reset_index()
            )
        merged = pd.merge(tiers['bulk'], tiers['hap1'], on='gene', how='outer')
        merged = pd.merge(merged, tiers['hap2'], on='gene', how='outer')
        return merged.drop_duplicates()

    mod_junction_df = build_phased_junction_df(junction_df, '')
    if cohort_junction_df is not None:
        mod_cohort_junction_df = build_phased_junction_df(cohort_junction_df, 'cohort_')
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
        # fallback below when --omim is omitted).
        for phasing in ('bulk', 'hap1', 'hap2'):
            for suffix in ('jxns', 'jxn_coverage', 'deltaPSI', 'jxn_annotation', 'jxn_event', 'jxn_nsamples'):
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
        """ Extract the max deltaPSI values for a gene """
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
        bulk_vals = parse_vals(row.get(prefix + 'bulk_deltaPSI'))
        hap1_vals = parse_vals(row.get(prefix + 'hap1_deltaPSI'))
        hap2_vals = parse_vals(row.get(prefix + 'hap2_deltaPSI'))
        max_bulk = max(abs(v) for v in bulk_vals) if bulk_vals else np.nan
        max_hap1 = max(abs(v) for v in hap1_vals) if hap1_vals else np.nan
        max_hap2 = max(abs(v) for v in hap2_vals) if hap2_vals else np.nan

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

    # Rank hits
    hit_df.sort_values(
        by=['ASE', 'outlier_junction', 'variant', 'bulk_jxn_coverage', 'hap1_jxn_coverage', 'hap2_jxn_coverage'],
        key=lambda col: (
            col.map({True: 0, False: 1, '.': 2}) if col.name == 'ASE' else
            col.map({'Strong': 0, 'Moderate': 1, 'Weak': 2, '.': 3}) if col.name == 'outlier_junction' else
            col.map({True: 0, False: 1, '.': 2}) if col.name == 'variant' else
            col.apply(lambda x: max([float(v) for v in re.split('[,;]', str(x)) if v.strip()]) if pd.notna(x) and x != '.' else 0)
        ),
        ascending=[True, True, True, False, False, False],
        inplace=True
    )
    hit_df['ranking'] = np.arange(1, len(hit_df) + 1)

    # Reorder
    hit_df['sample']=args.sample_name
    hit_df=hit_df[[
        'sample', 'gene', 'phenotypes', 'inheritance_patterns', 'ranking', 'variant', 'pathogenic_variant', 'ASE', 'outlier_junction', 'cohort_outlier_junction',
        'variant_ID', 'variant_GT', 'variant_gnomAD_AF',  'variant_CLNSIG', 'variant_CADD_PHRED', 'variant_SpliceAI', 'variant_consequence', 'variant_num_callers', 'variant_nsamples',
        'ASE_ratio', 'ASE_nsamples', 
        'bulk_jxns', 'bulk_jxn_coverage', 'bulk_deltaPSI', 'bulk_jxn_annotation', 'bulk_jxn_event', 'bulk_jxn_nsamples',
        'hap1_jxns', 'hap1_jxn_coverage', 'hap1_deltaPSI', 'hap1_jxn_annotation', 'hap1_jxn_event', 'hap1_jxn_nsamples',
        'hap2_jxns', 'hap2_jxn_coverage', 'hap2_deltaPSI', 'hap2_jxn_annotation', 'hap2_jxn_event', 'hap2_jxn_nsamples',
        'cohort_bulk_jxns', 'cohort_bulk_jxn_coverage', 'cohort_bulk_deltaPSI', 'cohort_bulk_jxn_annotation', 'cohort_bulk_jxn_event', 'cohort_bulk_jxn_nsamples',
        'cohort_hap1_jxns', 'cohort_hap1_jxn_coverage', 'cohort_hap1_deltaPSI', 'cohort_hap1_jxn_annotation', 'cohort_hap1_jxn_event', 'cohort_hap1_jxn_nsamples',
        'cohort_hap2_jxns', 'cohort_hap2_jxn_coverage', 'cohort_hap2_deltaPSI', 'cohort_hap2_jxn_annotation', 'cohort_hap2_jxn_event', 'cohort_hap2_jxn_nsamples'
    ]]

    # Make output directory if it doesn't exist
    if args.outfile:
        outdir = os.path.dirname(args.outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    # Save hits
    hit_df.to_csv(args.outfile, sep='\t', index=False)
    print(f"Saved merged hits to {args.outfile}")

if __name__ == "__main__":
    main()