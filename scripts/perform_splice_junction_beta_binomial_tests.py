#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2024.01.17
# Adapted from Robert Wang (Xing Lab)
# Optimized: 2025

import os, argparse, warnings
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import betabinom, beta
from pandas.errors import PerformanceWarning
from math import ceil
import concurrent.futures
from statsmodels.stats.multitest import multipletests
import time

warnings.filterwarnings('ignore', category=PerformanceWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Identifies splice junctions mapping to a user-defined gene region that '
                    'show unusually high or low usage frequencies within a sample of interest '
                    'relative to tissue-matched GTEx controls')
    parser.add_argument('--jxn-info-file', required=True,
        help='Path to TSV file with junction information for the region of interest.')
    parser.add_argument('--gtexfile', required=True,
        help='Path to read count matrix for splice junctions discovered across tissue-matched GTEx controls.')
    parser.add_argument('--outfile', required=True, type=str,
        help='Output file to write merged junctions.')
    parser.add_argument('--sample-coverage-threshold', type=int, default=20)
    parser.add_argument('--gtex-coverage-threshold', type=int, default=20)
    parser.add_argument('--PSI-rescale-factor', type=float, default=1e-3)
    parser.add_argument('--gtex-n-threshold', type=int, default=100)
    parser.add_argument('--phasing-threshold', default=0.8)
    parser.add_argument('--junction-to-gene-coverage-ratio', default=0, type=float)
    parser.add_argument('--annotation-file', type=str)
    parser.add_argument('--threads', type=int, default=1)
    return parser.parse_args()


def parse_gtf_splice_junctions(gtf_file):
    """Parse a GTF file and return a dict of splice junctions (chr_ss1_ss2) to annotation type."""

    transcripts = defaultdict(list)
    transcript_type = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, start, end, attrs = fields[0], int(fields[3]), int(fields[4]), fields[8]
            tid = ttype = None
            canonical_flag = False
            for attr in attrs.split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                parts = attr.split(" ", 1)
                if len(parts) != 2:
                    continue
                key, val = parts
                val = val.strip('"')
                if key == "transcript_id":
                    tid = val
                elif key == "transcript_type":
                    ttype = val
                elif key == "tag" and val == "Ensembl_canonical":
                    canonical_flag = True
            if tid:
                transcripts[tid].append((chrom, start, end))
                if ttype:
                    transcript_type[tid] = (
                        f"canonical:{ttype}" if canonical_flag else f"annotated:{ttype}"
                    )

    print(f"Parsed {len(transcripts)} transcripts from GTF.")
    print(f"Transcript types found: {set(transcript_type.values())}")

    priority_order = [
        "canonical:protein_coding",
        "canonical:protein_coding_CDS_not_defined",
        "canonical:lncRNA",
        "canonical:nonsense_mediated_decay",
        "canonical:non_stop_decay",
        "canonical:retained_intron",
        "canonical:other"
        "annotated:protein_coding",
        "annotated:protein_coding_CDS_not_defined",
        "annotated:lncRNA",
        "annotated:nonsense_mediated_decay",
        "annotated:non_stop_decay",
        "annotated:retained_intron",
        "annotated:other"
    ]
    priority_rank = {typ: i for i, typ in enumerate(priority_order)}

    # Normalise unknown transcript types to "other"
    for tid in list(transcript_type):
        if transcript_type[tid] not in priority_rank:
            transcript_type[tid] = "other"

    junction_types = defaultdict(set)
    for tid, exons in transcripts.items():
        if len(exons) > 1:
            exons_sorted = sorted(exons, key=lambda x: x[1])
            ttype = transcript_type.get(tid, "other")
            for (chrom, _, end1), (_, start2, _) in zip(exons_sorted, exons_sorted[1:]):
                junction_types[f"{chrom}_{end1+1}_{start2-1}"].add(ttype)

    junction_annotation = {
        jxn: min(types, key=lambda x: priority_rank.get(x, len(priority_order)))
        for jxn, types in junction_types.items()
    }
    return junction_annotation


def calculate_coverage(df, col_prefix=''):
    """Calculate coverage for each splice site and junction."""
    try:
        jxn_col = col_prefix + 'jxn_alignment_count'
        df[jxn_col] = pd.to_numeric(df[jxn_col], errors='coerce').fillna(0).astype(int)
        ss1_sums = df.groupby('ss1')[jxn_col].sum().astype(int)
        ss2_sums = df.groupby('ss2')[jxn_col].sum().astype(int)
        df[col_prefix + 'ss1_coverage'] = df['ss1'].map(ss1_sums).fillna(0).astype(int)
        df[col_prefix + 'ss2_coverage'] = df['ss2'].map(ss2_sums).fillna(0).astype(int)
        df[col_prefix + 'jxn_coverage'] = (
            df[col_prefix + 'ss1_coverage'] + df[col_prefix + 'ss2_coverage'] - df[jxn_col]
        ).astype(int)
    except Exception as e:
        print(f"Error calculating coverage: {e}")


def calculate_PSI(df, PSI_rescale_factor, col_prefix=''):
    """Calculate PSI values."""
    try:
        jxn_col = col_prefix + 'jxn_alignment_count'
        cov_col = col_prefix + 'jxn_coverage'
        df[jxn_col] = pd.to_numeric(df[jxn_col], errors='coerce').fillna(0).astype(int)
        df[cov_col] = pd.to_numeric(df[cov_col], errors='coerce').fillna(0).astype(int)
        # Vectorised PSI — avoids row-wise apply
        cov_num = df[cov_col].astype(float)
        jxn_num = df[jxn_col].astype(float)
        psi = np.where(cov_num == 0, np.nan, jxn_num / cov_num)
        df[col_prefix + 'sample_PSI'] = np.where(np.isnan(psi), "n/a", psi)
        df[col_prefix + 'rescaled_sample_PSI'] = np.where(
            np.isnan(psi), "n/a",
            psi * (1 - 2 * PSI_rescale_factor) + PSI_rescale_factor
        )
    except Exception as e:
        print(f"Error calculating PSI: {e}")


def fit_beta_dist(x, tol, n_threshold):
    """Fit a beta distribution on values in x."""
    x = pd.to_numeric(x, errors='coerce')
    x = x[~np.isnan(x)]

    if len(x) < n_threshold:
        return (len(x), "low_n", "low_n", "low_n")
    if x.var() < tol:
        return (len(x), x.mean() / tol, (1 - x.mean()) / tol, x.mean())
    try:
        alpha_value, beta_value = beta.fit(x, floc=0, fscale=1)[0:2]
        expected_PSI = alpha_value / (alpha_value + beta_value)
        return (len(x), alpha_value, beta_value, expected_PSI)
    except Exception:
        return (len(x), "error", "error", "error")


def beta_binomial_test(x, n, alpha_value, beta_value):
    """Compute the probability of observing a value as extreme as x from a beta distribution."""
    if any(np.isnan(val) for val in [x, n, alpha_value, beta_value]):
        return "n/a"
    x = round(x)
    n = round(n)
    lte_x = betabinom.cdf(x, n, alpha_value, beta_value)
    gte_x = betabinom.cdf(n - x, n, beta_value, alpha_value)
    p_value = np.clip(2 * min(lte_x, gte_x), 0, 1)
    if np.isnan(p_value):
        return "error"
    return p_value


def process_region(jxn_info_df_filtered, gtex_df_filtered, region, sample_coverage_threshold,
                   gtex_coverage_threshold, PSI_rescale_factor, gtex_n_threshold,
                   phasing_threshold, annotated_junctions, report_outdir):
    """Process region of interest."""

    start_time = time.time()
    report = os.path.join(report_outdir, f"{region.replace(':', '_').replace('-', '_')}_report.tsv")

    with open(report, 'w') as report_file:
        report_file.write(f"Processing region {region}...\n\n")

        ############################## STEP 1: INITIALIZE DATAFRAMES ##############################

        bulk_df = jxn_info_df_filtered[jxn_info_df_filtered['phasing'] == 'bulk'].copy()
        phasing_values = jxn_info_df_filtered['phasing'].unique()
        haplotype_specific = 'hap1' in phasing_values and 'hap2' in phasing_values
        if haplotype_specific:
            hap1_df = jxn_info_df_filtered[jxn_info_df_filtered['phasing'] == 'hap1'].copy()
            hap2_df = jxn_info_df_filtered[jxn_info_df_filtered['phasing'] == 'hap2'].copy()
        del jxn_info_df_filtered

        gtex_samples = gtex_df_filtered.columns[1:]
        gtex_df_filtered.columns = [gtex_df_filtered.columns[0]] + \
                                    [col + '_jxn_alignment_count' for col in gtex_samples]

        if len(bulk_df['gene'].unique()) > 1:
            report_file.write(f"Error: Multiple genes for region {region}. Exiting...\n")
            return
        if len(bulk_df['gene_alignment_count'].unique()) > 1:
            report_file.write(f"Error: Multiple gene alignment counts for region {region}. Exiting...\n")
            return
        gene = bulk_df['gene'].unique()[0]
        bulk_alignment_count = bulk_df['gene_alignment_count'].unique()[0]

        report_file.write(f"Gene: {gene}\n")
        report_file.write(f"\nAlignment count in sample of interest: {bulk_alignment_count}\n")
        if haplotype_specific:
            hap1_alignment_count = hap1_df['gene_alignment_count'].unique()[0]
            hap2_alignment_count = hap2_df['gene_alignment_count'].unique()[0]
            report_file.write(f"  Haplotype 1 alignment count: {hap1_alignment_count}\n")
            report_file.write(f"  Haplotype 2 alignment count: {hap2_alignment_count}\n")
        report_file.write(f"Number of junctions in sample of interest: {len(bulk_df)}\n\n")
        report_file.write(f"Number of junctions in GTEx: {len(gtex_df_filtered)}\n\n")
        if len(gtex_df_filtered) == 0:
            report_file.write(f"\nNo junctions found over region {region} in GTEx samples. Exiting...\n")
            return

        ############################## STEP 2: ADD MISSING JUNCTIONS ##############################

        # Vectorised splice-site extraction — faster than list comprehension on the index
        def _add_ss(df):
            parts = df.index.str.split('_')
            df['ss1'] = parts.map(lambda p: p[0] + '_' + p[1])
            df['ss2'] = parts.map(lambda p: p[0] + '_' + p[2])

        _add_ss(bulk_df)
        _add_ss(gtex_df_filtered)
        if haplotype_specific:
            _add_ss(hap1_df)
            _add_ss(hap2_df)

        bulk_df_full = pd.DataFrame()
        for sample in bulk_df['sample'].unique():
            bulk_df_sample = bulk_df[bulk_df['sample'] == sample]
            bulk_df_sample = pd.concat([
                bulk_df_sample,
                gtex_df_filtered[~gtex_df_filtered.index.isin(bulk_df_sample.index)][['ss1', 'ss2']]
                    .assign(sample=sample, phasing='bulk', region=region, gene=gene,
                            gene_alignment_count=bulk_alignment_count, jxn_alignment_count=0)
            ])
            bulk_df_full = pd.concat([bulk_df_full, bulk_df_sample])

        gtex_df_full = pd.concat([
            gtex_df_filtered,
            bulk_df[~bulk_df.index.isin(gtex_df_filtered.index)][['ss1', 'ss2']]
                .assign(**{s + '_jxn_alignment_count': 0 for s in gtex_samples})
        ])
        del bulk_df, bulk_df_sample

        if haplotype_specific:
            hap1_df_full = pd.DataFrame()
            for sample in hap1_df['sample'].unique():
                hap1_df_sample = hap1_df[hap1_df['sample'] == sample]
                hap1_df_sample = pd.concat([
                    hap1_df_sample,
                    gtex_df_full[~gtex_df_full.index.isin(hap1_df_sample.index)][['ss1', 'ss2']]
                        .assign(sample=sample, phasing='hap1', region=region, gene=gene,
                                gene_alignment_count=hap1_alignment_count, jxn_alignment_count=0)
                ])
                hap1_df_full = pd.concat([hap1_df_full, hap1_df_sample])

            hap2_df_full = pd.DataFrame()
            for sample in hap2_df['sample'].unique():
                hap2_df_sample = hap2_df[hap2_df['sample'] == sample]
                hap2_df_sample = pd.concat([
                    hap2_df_sample,
                    gtex_df_full[~gtex_df_full.index.isin(hap2_df_sample.index)][['ss1', 'ss2']]
                        .assign(sample=sample, phasing='hap2', region=region, gene=gene,
                                gene_alignment_count=hap2_alignment_count, jxn_alignment_count=0)
                ])
                hap2_df_full = pd.concat([hap2_df_full, hap2_df_sample])
            del hap1_df, hap2_df, hap1_df_sample, hap2_df_sample

        if bulk_df_full.index.symmetric_difference(gtex_df_full.index).any():
            report_file.write(f"Error: unexpected mismatch between junctions in sample of interest and GTEx samples. Exiting...\n")
            return
        report_file.write(f"There are {len(bulk_df_full)} total junctions to analyze.\n\n")

        ############################## STEP 3: CALCULATE COVERAGE ##############################

        report_file.write(f"Calculating coverage...\n")
        calculate_coverage(bulk_df_full)
        if haplotype_specific:
            calculate_coverage(hap1_df_full)
            calculate_coverage(hap2_df_full)
        for sample in gtex_samples:
            calculate_coverage(gtex_df_full, sample + '_')

        ############################## STEP 4: CALCULATE PSI VALUES ##############################

        report_file.write(f"Calculating PSI values...\n")
        calculate_PSI(bulk_df_full, PSI_rescale_factor)
        if haplotype_specific:
            calculate_PSI(hap1_df_full, PSI_rescale_factor)
            calculate_PSI(hap2_df_full, PSI_rescale_factor)
        for sample in gtex_samples:
            calculate_PSI(gtex_df_full, PSI_rescale_factor, sample + '_')
            low_coverage_mask = pd.to_numeric(gtex_df_full[sample + '_jxn_coverage'], errors='coerce') < gtex_coverage_threshold
            gtex_df_full.loc[low_coverage_mask, sample + '_rescaled_sample_PSI'] = np.nan

        gtex_df_rescaled_PSI = gtex_df_full[[f'{s}_rescaled_sample_PSI' for s in gtex_samples]].copy()
        del gtex_df_full

        ############################## STEP 5: FIT BETA DISTRIBUTION ON GTEx DATA ##############################

        report_file.write(f"Fitting beta distributions on PSI values from GTEx samples...\n")
        gtex_df_rescaled_PSI[['num_gtex_samples_with_good_coverage', 'alpha', 'beta', 'expected_PSI']] = \
            gtex_df_rescaled_PSI.apply(
                lambda row: pd.Series(fit_beta_dist(np.array(row), tol=PSI_rescale_factor, n_threshold=gtex_n_threshold)),
                axis=1
            )
        final_df = bulk_df_full.merge(
            gtex_df_rescaled_PSI[['num_gtex_samples_with_good_coverage', 'alpha', 'beta', 'expected_PSI']],
            left_index=True, right_index=True
        )

        ############################## STEP 6: RUN BETA-BINOMIAL TESTS ##############################

        if haplotype_specific:
            hap1_coverage = hap1_df_full['jxn_coverage'].reindex(final_df.index, fill_value=0)
            hap2_coverage = hap2_df_full['jxn_coverage'].reindex(final_df.index, fill_value=0)
            coverage_mask = (hap1_coverage + hap2_coverage) > (phasing_threshold * final_df['jxn_coverage'])
            hap1_df_filtered = hap1_df_full[coverage_mask.reindex(hap1_df_full.index, fill_value=False)].merge(
                final_df[['num_gtex_samples_with_good_coverage', 'alpha', 'beta', 'expected_PSI']],
                left_index=True, right_index=True)
            hap2_df_filtered = hap2_df_full[coverage_mask.reindex(hap2_df_full.index, fill_value=False)].merge(
                final_df[['num_gtex_samples_with_good_coverage', 'alpha', 'beta', 'expected_PSI']],
                left_index=True, right_index=True)
            final_df = pd.concat([final_df, hap1_df_filtered, hap2_df_filtered])
            num_hap = len(final_df[final_df['phasing'].isin(['hap1', 'hap2'])])
            report_file.write(f"Running beta-binomial tests for {len(final_df)-num_hap} bulk "
                               f"and {num_hap} haplotype-specific junctions.\n")
        else:
            report_file.write(f"Running beta-binomial tests for {len(final_df)} junctions.\n")

        # Vectorised numeric coercion before row-wise apply
        jxn_num = pd.to_numeric(final_df['jxn_alignment_count'], errors='coerce')
        cov_num = pd.to_numeric(final_df['jxn_coverage'], errors='coerce')
        alpha_num = pd.to_numeric(final_df['alpha'], errors='coerce')
        beta_num = pd.to_numeric(final_df['beta'], errors='coerce')
        final_df['p_value'] = [
            beta_binomial_test(ceil(x), ceil(n), a, b)
            for x, n, a, b in zip(jxn_num, cov_num, alpha_num, beta_num)
        ]

        psi_num = pd.to_numeric(final_df['rescaled_sample_PSI'], errors='coerce')
        exp_num = pd.to_numeric(final_df['expected_PSI'], errors='coerce')
        final_df['delta_PSI'] = np.where(
            np.isnan(psi_num) | np.isnan(exp_num), "n/a",
            psi_num - exp_num
        )

        ############################## STEP 6: ADD FLAGS ##############################

        final_df['flag'] = ""
        final_df['flag'] = np.where(
            pd.to_numeric(final_df['jxn_alignment_count'], errors='coerce') == 0,
            final_df['flag'] + ';not_detected_in_sample', final_df['flag'])
        final_df['flag'] = np.where(
            (0 < pd.to_numeric(final_df['jxn_alignment_count'], errors='coerce')) &
            (pd.to_numeric(final_df['jxn_alignment_count'], errors='coerce') < 5),
            final_df['flag'] + ';low_jxn_alignment_count', final_df['flag'])
        final_df['flag'] = np.where(
            (0 < pd.to_numeric(final_df['jxn_alignment_count'], errors='coerce')) &
            (0 < pd.to_numeric(final_df['jxn_coverage'], errors='coerce')) &
            (pd.to_numeric(final_df['jxn_coverage'], errors='coerce') < sample_coverage_threshold),
            final_df['flag'] + ';low_coverage', final_df['flag'])
        final_df['flag'] = np.where(
            ~final_df.index.isin(gtex_df_filtered.index),
            final_df['flag'] + ';not_detected_in_gtex', final_df['flag'])

        final_df['flag'] = final_df['flag'].str[1:]  # strip leading semicolon
        final_df['flag'] = final_df['flag'].apply(lambda x: x if x else "no_flag")

        ############################## STEP 7: IDENTIFY ANNOTATED JUNCTIONS ##############################

        if annotated_junctions:
            final_df["annotation"] = final_df.index.map(annotated_junctions).fillna("unannotated")
            n_annotated = (final_df["annotation"] != "unannotated").sum()
            n_novel = (final_df["annotation"] == "unannotated").sum()
            report_file.write(f"Identified {n_annotated} annotated junctions and {n_novel} novel junctions in region {region}.\n")
        else:
            final_df["annotation"] = "n/a"

        ############################## STEP 8: SAVE RESULTS ##############################

        final_df.reset_index(inplace=True)
        final_df.rename(columns={'index': 'junction'}, inplace=True)
        final_df.fillna("n/a", inplace=True)
        final_df = final_df[['sample', 'phasing', 'region', 'gene', 'gene_alignment_count', 'junction',
                              'jxn_alignment_count', 'ss1_coverage', 'ss2_coverage', 'jxn_coverage',
                              'sample_PSI', 'rescaled_sample_PSI', 'num_gtex_samples_with_good_coverage',
                              'alpha', 'beta', 'expected_PSI', 'delta_PSI', 'p_value', 'flag', 'annotation']]

        report_file.write(f"Region {region} processed successfully in {time.time() - start_time:.2f} seconds.\n\n")

    return final_df


def main():
    """Main script."""

    print(f"\n\n\n******************************************************************************************")
    print(f"Detecting splice junction outliers...")
    print(f"******************************************************************************************\n")

    args = parse_args()

    if not os.path.exists(args.jxn_info_file):
        print(f"\nERROR: {args.jxn_info_file} not found.")
        return
    if not os.path.exists(args.gtexfile):
        print(f"\nERROR: {args.gtexfile} not found.")
        return

    outdir = os.path.dirname(args.outfile)
    os.makedirs(outdir, exist_ok=True)
    report_outdir = os.path.join(outdir, 'reports')
    os.makedirs(report_outdir, exist_ok=True)

    jxn_info_df = pd.read_csv(args.jxn_info_file, sep='\t', keep_default_na=False, header=0,
                               dtype={'sample': str, 'phasing': str, 'region': str, 'gene': str,
                                      'region_alignment_count': int, 'junction': str,
                                      'jxn_alignment_count': int})

    print(f"\nReading GTEx file {args.gtexfile}...")
    gtex_df = pd.read_csv(args.gtexfile, sep='\t', index_col=0)

    if args.annotation_file:
        if not os.path.exists(args.annotation_file):
            print(f"\nERROR: {args.annotation_file} not found.")
            return
        print(f"\nReading GTF file {args.annotation_file}...")
        annotated_junctions = parse_gtf_splice_junctions(args.annotation_file)
        print(f"Identified {len(annotated_junctions)} annotated junctions in the GTF file.")
    else:
        annotated_junctions = {}

    regions = jxn_info_df['region'].unique()
    print(f"\nBegin processing {len(regions)} regions using {args.threads} threads. This may take a while...")

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for region in regions:
            jxn_info_df_filtered = jxn_info_df[jxn_info_df['region'] == region].copy()
            jxn_info_df_filtered.set_index('junction', inplace=True)

            # Pre-parse region components for the GTEx filter — avoids per-row lambda parsing
            reg_chrom, reg_coords = region.split(':')
            reg_start, reg_end = map(int, reg_coords.split('-'))
            # GTEx junction IDs are formatted "chr:start-end:strand" (e.g. "chr1:11212-12009:+"),
            # unlike every other junction ID in this pipeline (get_splice_junction_counts_by_region.py,
            # make_junction_count_matrix.py, merge_and_filter_junction_results.py), which all use
            # "chr_start_end". Strand is ignored -- start is always < end regardless of strand.
            idx_chrom_coord = gtex_df.index.str.split(':')
            gtex_chrom = idx_chrom_coord.map(lambda p: p[0])
            idx_coords = idx_chrom_coord.map(lambda p: p[1]).str.split('-')
            gtex_start = idx_coords.map(lambda p: int(p[0]))
            gtex_end = idx_coords.map(lambda p: int(p[1]))
            gtex_mask = (
                (gtex_chrom == reg_chrom) &
                (gtex_start >= reg_start) &
                (gtex_end <= reg_end)
            )
            gtex_df_filtered = gtex_df[gtex_mask].copy()
            # Normalize to "chr_start_end" so downstream index comparisons against sample
            # junctions (_add_ss, .isin() in process_region) match correctly instead of
            # silently never matching two different string formats for the same junction.
            gtex_df_filtered.index = (
                gtex_chrom[gtex_mask].astype(str) + '_' +
                gtex_start[gtex_mask].astype(str) + '_' +
                gtex_end[gtex_mask].astype(str)
            )

            futures.append(executor.submit(
                process_region, jxn_info_df_filtered, gtex_df_filtered, region,
                args.sample_coverage_threshold, args.gtex_coverage_threshold,
                args.PSI_rescale_factor, args.gtex_n_threshold, args.phasing_threshold,
                annotated_junctions, report_outdir))

    results = []
    for future in concurrent.futures.as_completed(futures):
        try:
            res = future.result()
            if res is not None:
                results.append(res)
        except Exception as e:
            print(f"Error encountered: {e}")

    df = pd.concat(results, ignore_index=True)

    mask = (
        (pd.to_numeric(df['num_gtex_samples_with_good_coverage'], errors='coerce') >= args.gtex_n_threshold) &
        (pd.to_numeric(df['jxn_coverage'], errors='coerce') >= args.sample_coverage_threshold) &
        (pd.to_numeric(df['jxn_coverage'], errors='coerce') /
         pd.to_numeric(df['gene_alignment_count'], errors='coerce') >= args.junction_to_gene_coverage_ratio)
    )

    print('Correcting for multiple testing...')
    col_idx = df.columns.get_loc("p_value")
    df.insert(col_idx + 1, "padj", "n/a")
    df['padj'] = df['padj'].astype(object)
    numeric_p_values = pd.to_numeric(df.loc[mask, 'p_value'], errors='coerce')
    padj_values = multipletests(numeric_p_values.dropna(), method='fdr_by')[1]
    df.loc[mask & numeric_p_values.notna(), 'padj'] = padj_values

    df = df.sort_values(by=['sample', 'junction', 'phasing'])
    df.to_csv(args.outfile, sep='\t', index=False)
    print(f'Merged junctions written to {args.outfile}')


if __name__ == '__main__':
    main()
