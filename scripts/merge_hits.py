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
                                                                            chrom, pos, ref, alt, GT, gnomAD_AF, CLNSIG, gene, CADD_PHRED, SpliceAI, num_callers, sample_count')
    parser.add_argument('--ase-hits', required=True, type=str, help='Path to ASE hit results. At a minimum, file should contain columns: gene, ratio')
    parser.add_argument('--junction-hits', required=True, type=str, help='Path to junction hit results. At a minimum, file should contain columns: \
                                                                            gene, phasing, junction, delta_PSI, sample_count')
    parser.add_argument('--omim', required=True, type=str, help='Path to OMIM data. At a minimum, file should contain columns: approved_gene_symbol, phenotypes, inheritance_patterns')
    
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
    variant_df = pd.read_csv(args.variant_hits, sep='\t', usecols=['chrom', 'pos', 'ref', 'alt', 'GT', 'gnomAD_AF', 'CLNSIG', 'gene', 'CADD_PHRED', 'SpliceAI', 'num_callers', 'sample_count']).drop_duplicates()
    variant_df = variant_df.rename(columns={'sample_count': 'variant_nsamples'})
    variant_df = variant_df[variant_df['gene'] != '.']
    # 'gene' lists every BED-panel gene a variant overlaps (comma-separated
    # if more than one) -- explode so a variant overlapping multiple genes
    # is grouped into each of those genes' hit rows individually.
    variant_df['gene'] = variant_df['gene'].str.split(',')
    variant_df = variant_df.explode('gene')
    variant_df['GT'] = variant_df['GT'].str.replace('|', '/', regex=False).str.replace('1/0', '0/1', regex=False)
    variant_df['variant_ID'] = variant_df.apply(lambda x: f"{x.chrom}-{x.pos}-{x.ref}-{x.alt}", axis=1).drop_duplicates()
    ase_df = pd.read_csv(args.ase_hits, sep='\t', usecols=['gene', 'ratio', 'sample_count']).drop_duplicates()
    ase_df = ase_df.rename(columns={'ratio': 'ASE_ratio', 'sample_count': 'ASE_nsamples'})
    junction_df = pd.read_csv(args.junction_hits, sep='\t', usecols=['gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count']).drop_duplicates()
    omim_df = pd.read_csv(args.omim, sep='\t', usecols=['approved_gene_symbol', 'phenotypes', 'inheritance_patterns'])
    omim_df = omim_df.rename(columns={'approved_gene_symbol': 'gene'})

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
                variant_num_callers=('num_callers', ';'.join),
                variant_nsamples=('variant_nsamples', ';'.join),
            )
            .reset_index()
    )
    # mod_junction_df should contain genes, bulk_jxns, bulk_jxn_coverage, bulk_deltaPSI, bulk_jxn_nsamples, 
    # hap1_jxns, hap1_jxn_coverage, hap1_deltaPSI, hap1_jxn_nsamples, hap2_jxns, hap2_jxn_coverage, hap2_delta_PSI, hap2_jxn_nsamples
    bulk_junction_df = (
        junction_df[junction_df['phasing'] == 'bulk']
            .sort_values('junction')
            .groupby('gene')
            .agg(
                bulk_jxns=('junction', lambda x: ';'.join(map(str, x))),
                bulk_jxn_coverage=('jxn_coverage', lambda x: ';'.join(map(str, x))),
                bulk_deltaPSI=('delta_PSI', lambda x: ';'.join(map(str, x))),
                bulk_jxn_nsamples=('sample_count', lambda x: ';'.join(map(str, x))),
            )
            .reset_index()
    )
    hap1_junction_df = (
        junction_df[junction_df['phasing'] == 'hap1']
            .sort_values('junction')
            .groupby('gene')
            .agg(
                hap1_jxns=('junction', lambda x: ','.join(map(str, x))),
                hap1_jxn_coverage=('jxn_coverage', lambda x: ','.join(map(str, x))),
                hap1_deltaPSI=('delta_PSI', lambda x: ','.join(map(str, x))),
                hap1_jxn_nsamples=('sample_count', lambda x: ','.join(map(str, x)))
            )
            .reset_index()
    )
    hap2_junction_df = (
        junction_df[junction_df['phasing'] == 'hap2']
            .sort_values('junction')
            .groupby('gene')
            .agg(
                hap2_jxns=('junction', lambda x: ','.join(map(str, x))),
                hap2_jxn_coverage=('jxn_coverage', lambda x: ','.join(map(str, x))),
                hap2_deltaPSI=('delta_PSI', lambda x: ','.join(map(str, x))),
                hap2_jxn_nsamples=('sample_count', lambda x: ','.join(map(str, x)))
            )
            .reset_index()
    )
    mod_junction_df = pd.merge(bulk_junction_df, hap1_junction_df, on='gene', how='outer')
    mod_junction_df = pd.merge(mod_junction_df, hap2_junction_df, on='gene', how='outer')
    mod_junction_df = mod_junction_df.drop_duplicates()

    # Merge hits
    hit_df = pd.merge(mod_variant_df, ase_df, on='gene', how='outer')
    hit_df = pd.merge(hit_df, mod_junction_df, on='gene', how='outer')
    hit_df = pd.merge(hit_df, omim_df, on='gene', how='left')
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

    def max_deltas(row):
        """ Extract the max deltaPSI values for a gene """
        def parse_vals(s):
            if pd.notna(s):
                return [float(v) for v in re.split('[,;]', s) if v.strip()]
            return []
        bulk_vals = parse_vals(row.get('bulk_deltaPSI'))
        hap1_vals = parse_vals(row.get('hap1_deltaPSI'))
        hap2_vals = parse_vals(row.get('hap2_deltaPSI'))
        max_bulk = max(abs(v) for v in bulk_vals) if bulk_vals else np.nan
        max_hap1 = max(abs(v) for v in hap1_vals) if hap1_vals else np.nan
        max_hap2 = max(abs(v) for v in hap2_vals) if hap2_vals else np.nan

        return max_bulk, max_hap1, max_hap2
    def inspect_row(row):
        """ Determine if a gene has splicing dysregulation """
        max_bulk, max_hap1, max_hap2 = max_deltas(row)
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

    # Fill missing values
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
        'sample', 'gene', 'phenotypes', 'inheritance_patterns', 'ranking', 'variant', 'pathogenic_variant', 'ASE', 'outlier_junction', 
        'variant_ID', 'variant_GT', 'variant_gnomAD_AF',  'variant_CLNSIG', 'variant_CADD_PHRED', 'variant_SpliceAI', 'variant_num_callers', 'variant_nsamples',
        'ASE_ratio', 'ASE_nsamples', 
        'bulk_jxns', 'bulk_jxn_coverage', 'bulk_deltaPSI', 'bulk_jxn_nsamples', 'hap1_jxns', 'hap1_jxn_coverage', 'hap1_deltaPSI', 'hap1_jxn_nsamples', 'hap2_jxns', 'hap2_jxn_coverage', 'hap2_deltaPSI', 'hap2_jxn_nsamples'
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
