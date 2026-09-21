

import os, argparse, warnings, traceback
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import betabinom
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
                    'relative to tissue-matched GTEx controls, using precomputed beta fits '
                    '(see fit_gtex_beta_distributions.py / fit_novel_junction_beta_distributions.py).')
    parser.add_argument('--jxn-info-file', required=True,
        help='Path to TSV file with junction information for the region of interest.')
    parser.add_argument('--gtex-beta-fits', required=True,
        help='Path to fit_gtex_beta_distributions.py\'s (_5B1) output for this tissue.')
    parser.add_argument('--novel-beta-fits', required=True,
        help='Path to fit_novel_junction_beta_distributions.py\'s (_5B2) output for this sample+tissue.')
    parser.add_argument('--outfile', required=True, type=str,
        help='Output file to write merged junctions.')
    parser.add_argument('--sample-coverage-threshold', type=int, default=20)
    parser.add_argument('--gtex-coverage-threshold', type=int, default=20)
    parser.add_argument('--PSI-rescale-factor', type=float, default=1e-3)
    parser.add_argument('--gtex-n-threshold', type=int, default=100)
    parser.add_argument('--phasing-threshold', type=float, default=0.8)
    parser.add_argument('--junction-to-gene-coverage-ratio', default=0, type=float)
    parser.add_argument('--annotation-file', type=str)
    parser.add_argument('--threads', type=int, default=1)
    return parser.parse_args()


def parse_gtf_splice_junctions(gtf_file):

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
    try:
        jxn_col = col_prefix + 'jxn_alignment_count'
        cov_col = col_prefix + 'jxn_coverage'
        df[jxn_col] = pd.to_numeric(df[jxn_col], errors='coerce').fillna(0).astype(int)
        df[cov_col] = pd.to_numeric(df[cov_col], errors='coerce').fillna(0).astype(int)
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
        traceback.print_exc()


def beta_binomial_test_vectorized(x, n, alpha_value, beta_value):
    x = np.asarray(x, dtype=float)
    n = np.asarray(n, dtype=float)
    alpha_value = np.asarray(alpha_value, dtype=float)
    beta_value = np.asarray(beta_value, dtype=float)

    valid = ~(np.isnan(x) | np.isnan(n) | np.isnan(alpha_value) | np.isnan(beta_value))

    x_r = np.round(x)
    n_r = np.round(n)

    result = np.full(x.shape, "n/a", dtype=object)

    if valid.any():
        lte_x = betabinom.cdf(x_r[valid], n_r[valid], alpha_value[valid], beta_value[valid])
        gte_x = betabinom.cdf(n_r[valid] - x_r[valid], n_r[valid], beta_value[valid], alpha_value[valid])
        p_value = np.clip(2 * np.minimum(lte_x, gte_x), 0, 1)

        valid_result = np.where(np.isnan(p_value), "error", p_value)
        result[valid] = valid_result

    return result


def process_region(jxn_info_df_filtered, gtex_fits_filtered, region, sample_coverage_threshold,
                   PSI_rescale_factor, phasing_threshold, annotated_junctions, report_outdir):

    start_time = time.time()
    report = os.path.join(report_outdir, f"{region.replace(':', '_').replace('-', '_')}_report.tsv")

    with open(report, 'w') as report_file:
        report_file.write(f"Processing region {region}...\n\n")


        bulk_df = jxn_info_df_filtered[jxn_info_df_filtered['phasing'] == 'bulk'].copy()
        phasing_values = jxn_info_df_filtered['phasing'].unique()
        haplotype_specific = 'hap1' in phasing_values and 'hap2' in phasing_values
        if haplotype_specific:
            hap1_df = jxn_info_df_filtered[jxn_info_df_filtered['phasing'] == 'hap1'].copy()
            hap2_df = jxn_info_df_filtered[jxn_info_df_filtered['phasing'] == 'hap2'].copy()
        del jxn_info_df_filtered

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
        report_file.write(f"Number of GTEx-covered junctions in region: {len(gtex_fits_filtered)}\n\n")
        if len(gtex_fits_filtered) == 0:
            report_file.write(f"\nNo junctions found over region {region} in the GTEx beta-fit tables. Exiting...\n")
            return


        def _add_ss(df):
            parts = df.index.str.split('_')
            df['ss1'] = parts.map(lambda p: p[0] + '_' + p[1])
            df['ss2'] = parts.map(lambda p: p[0] + '_' + p[2])

        bulk_df = bulk_df.set_index('junction') if 'junction' in bulk_df.columns else bulk_df
        _add_ss(bulk_df)
        if haplotype_specific:
            hap1_df = hap1_df.set_index('junction') if 'junction' in hap1_df.columns else hap1_df
            hap2_df = hap2_df.set_index('junction') if 'junction' in hap2_df.columns else hap2_df
            _add_ss(hap1_df)
            _add_ss(hap2_df)

        bulk_df_full = pd.DataFrame()
        for sample in bulk_df['sample'].unique():
            bulk_df_sample = bulk_df[bulk_df['sample'] == sample]
            missing = gtex_fits_filtered[~gtex_fits_filtered.index.isin(bulk_df_sample.index)]
            padding = pd.DataFrame(index=missing.index)
            parts = padding.index.str.split('_')
            padding['ss1'] = parts.map(lambda p: p[0] + '_' + p[1])
            padding['ss2'] = parts.map(lambda p: p[0] + '_' + p[2])
            padding['sample'] = sample
            padding['phasing'] = 'bulk'
            padding['region'] = region
            padding['gene'] = gene
            padding['gene_alignment_count'] = bulk_alignment_count
            padding['jxn_alignment_count'] = 0
            bulk_df_sample = pd.concat([bulk_df_sample, padding])
            bulk_df_full = pd.concat([bulk_df_full, bulk_df_sample])
        del bulk_df, bulk_df_sample

        gtex_fits_full = pd.concat([
            gtex_fits_filtered,
            pd.DataFrame(
                index=bulk_df_full.index[~bulk_df_full.index.isin(gtex_fits_filtered.index)].unique()
            )
        ])
        gtex_fits_full = gtex_fits_full.reindex(bulk_df_full.index.unique())

        if haplotype_specific:
            hap1_df_full = pd.DataFrame()
            for sample in hap1_df['sample'].unique():
                hap1_df_sample = hap1_df[hap1_df['sample'] == sample]
                missing = gtex_fits_full[~gtex_fits_full.index.isin(hap1_df_sample.index)]
                padding = pd.DataFrame(index=missing.index)
                parts = padding.index.str.split('_')
                padding['ss1'] = parts.map(lambda p: p[0] + '_' + p[1])
                padding['ss2'] = parts.map(lambda p: p[0] + '_' + p[2])
                padding['sample'] = sample
                padding['phasing'] = 'hap1'
                padding['region'] = region
                padding['gene'] = gene
                padding['gene_alignment_count'] = hap1_alignment_count
                padding['jxn_alignment_count'] = 0
                hap1_df_sample = pd.concat([hap1_df_sample, padding])
                hap1_df_full = pd.concat([hap1_df_full, hap1_df_sample])

            hap2_df_full = pd.DataFrame()
            for sample in hap2_df['sample'].unique():
                hap2_df_sample = hap2_df[hap2_df['sample'] == sample]
                missing = gtex_fits_full[~gtex_fits_full.index.isin(hap2_df_sample.index)]
                padding = pd.DataFrame(index=missing.index)
                parts = padding.index.str.split('_')
                padding['ss1'] = parts.map(lambda p: p[0] + '_' + p[1])
                padding['ss2'] = parts.map(lambda p: p[0] + '_' + p[2])
                padding['sample'] = sample
                padding['phasing'] = 'hap2'
                padding['region'] = region
                padding['gene'] = gene
                padding['gene_alignment_count'] = hap2_alignment_count
                padding['jxn_alignment_count'] = 0
                hap2_df_sample = pd.concat([hap2_df_sample, padding])
                hap2_df_full = pd.concat([hap2_df_full, hap2_df_sample])
            del hap1_df, hap2_df, hap1_df_sample, hap2_df_sample

        report_file.write(f"There are {len(bulk_df_full)} total junctions to analyze.\n\n")


        report_file.write(f"Calculating coverage...\n")
        calculate_coverage(bulk_df_full)
        if haplotype_specific:
            calculate_coverage(hap1_df_full)
            calculate_coverage(hap2_df_full)


        report_file.write(f"Calculating PSI values...\n")
        calculate_PSI(bulk_df_full, PSI_rescale_factor)
        if haplotype_specific:
            calculate_PSI(hap1_df_full, PSI_rescale_factor)
            calculate_PSI(hap2_df_full, PSI_rescale_factor)


        report_file.write(f"Merging precomputed GTEx beta fits...\n")
        fit_cols = ['num_gtex_samples_with_good_coverage', 'alpha', 'beta', 'expected_PSI',
                    'p1_PSI', 'p99_PSI', 'in_gtex_matrix']
        final_df = bulk_df_full.merge(
            gtex_fits_full[fit_cols], left_index=True, right_index=True, how='left'
        )
        final_df['in_gtex_matrix'] = final_df['in_gtex_matrix'].fillna(False)


        if haplotype_specific:
            hap1_coverage = pd.to_numeric(hap1_df_full['jxn_coverage'], errors='coerce').reindex(final_df.index, fill_value=0)
            hap2_coverage = pd.to_numeric(hap2_df_full['jxn_coverage'], errors='coerce').reindex(final_df.index, fill_value=0)
            bulk_coverage_num = pd.to_numeric(final_df['jxn_coverage'], errors='coerce')
            coverage_mask = (hap1_coverage + hap2_coverage) > (phasing_threshold * bulk_coverage_num)
            hap1_df_filtered = hap1_df_full[coverage_mask.reindex(hap1_df_full.index, fill_value=False)].merge(
                gtex_fits_full[fit_cols], left_index=True, right_index=True, how='left')
            hap2_df_filtered = hap2_df_full[coverage_mask.reindex(hap2_df_full.index, fill_value=False)].merge(
                gtex_fits_full[fit_cols], left_index=True, right_index=True, how='left')
            hap1_df_filtered['in_gtex_matrix'] = hap1_df_filtered['in_gtex_matrix'].fillna(False)
            hap2_df_filtered['in_gtex_matrix'] = hap2_df_filtered['in_gtex_matrix'].fillna(False)
            final_df = pd.concat([final_df, hap1_df_filtered, hap2_df_filtered])
            num_hap = len(final_df[final_df['phasing'].isin(['hap1', 'hap2'])])
            report_file.write(f"Running beta-binomial tests for {len(final_df)-num_hap} bulk "
                               f"and {num_hap} haplotype-specific junctions.\n")
        else:
            report_file.write(f"Running beta-binomial tests for {len(final_df)} junctions.\n")

        jxn_num = pd.to_numeric(final_df['jxn_alignment_count'], errors='coerce')
        cov_num = pd.to_numeric(final_df['jxn_coverage'], errors='coerce')
        alpha_num = pd.to_numeric(final_df['alpha'], errors='coerce')
        beta_num = pd.to_numeric(final_df['beta'], errors='coerce')
        final_df['p_value'] = beta_binomial_test_vectorized(
            np.ceil(jxn_num), np.ceil(cov_num), alpha_num, beta_num
        )

        psi_num = pd.to_numeric(final_df['rescaled_sample_PSI'], errors='coerce')
        p1_num  = pd.to_numeric(final_df['p1_PSI'], errors='coerce')
        p99_num = pd.to_numeric(final_df['p99_PSI'], errors='coerce')
        delta_num = np.where(
            psi_num > p99_num, psi_num - p99_num,
            np.where(psi_num < p1_num, psi_num - p1_num, 0.0)
        )
        delta_PSI = pd.Series(delta_num, index=final_df.index).astype(object)
        delta_PSI[np.isnan(psi_num) | np.isnan(p1_num) | np.isnan(p99_num)] = "n/a"
        delta_PSI[final_df['expected_PSI'] == 'low_n'] = 'low_n'
        delta_PSI[final_df['expected_PSI'] == 'error'] = 'error'
        final_df['delta_PSI'] = delta_PSI


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
            ~final_df['in_gtex_matrix'].astype(bool),
            final_df['flag'] + ';not_detected_in_gtex', final_df['flag'])

        final_df['flag'] = final_df['flag'].str[1:]
        final_df['flag'] = final_df['flag'].apply(lambda x: x if x else "no_flag")


        if annotated_junctions:
            final_df["annotation"] = final_df.index.map(annotated_junctions).fillna("unannotated")
            n_annotated = (final_df["annotation"] != "unannotated").sum()
            n_novel = (final_df["annotation"] == "unannotated").sum()
            report_file.write(f"Identified {n_annotated} annotated junctions and {n_novel} novel junctions in region {region}.\n")
        else:
            final_df["annotation"] = "n/a"


        final_df.reset_index(inplace=True)
        final_df.rename(columns={'index': 'junction'}, inplace=True)
        final_df.fillna("n/a", inplace=True)
        final_df = final_df[['sample', 'phasing', 'region', 'gene', 'gene_alignment_count', 'junction',
                              'jxn_alignment_count', 'ss1_coverage', 'ss2_coverage', 'jxn_coverage',
                              'sample_PSI', 'rescaled_sample_PSI', 'num_gtex_samples_with_good_coverage',
                              'alpha', 'beta', 'expected_PSI', 'p1_PSI', 'p99_PSI', 'delta_PSI', 'p_value', 'flag', 'annotation']]

        report_file.write(f"Region {region} processed successfully in {time.time() - start_time:.2f} seconds.\n\n")

    return final_df


def main():

    print(f"\n\n\n******************************************************************************************")
    print(f"Testing splice junction usage against precomputed GTEx beta fits...")
    print(f"******************************************************************************************\n")

    args = parse_args()

    if not os.path.exists(args.jxn_info_file):
        print(f"\nERROR: {args.jxn_info_file} not found.")
        return
    if not os.path.exists(args.gtex_beta_fits):
        print(f"\nERROR: {args.gtex_beta_fits} not found.")
        return
    if not os.path.exists(args.novel_beta_fits):
        print(f"\nERROR: {args.novel_beta_fits} not found.")
        return

    outdir = os.path.dirname(args.outfile)
    os.makedirs(outdir, exist_ok=True)
    report_outdir = os.path.join(outdir, 'reports')
    os.makedirs(report_outdir, exist_ok=True)

    jxn_info_df = pd.read_csv(args.jxn_info_file, sep='\t', keep_default_na=False, header=0,
                               dtype={'sample': str, 'phasing': str, 'region': str, 'gene': str,
                                      'gene_alignment_count': int, 'junction': str,
                                      'jxn_alignment_count': int})

    print(f"\nReading precomputed GTEx beta fits {args.gtex_beta_fits}...")
    gtex_fits = pd.read_csv(args.gtex_beta_fits, sep='\t', comment='#', index_col=0)

    print(f"Reading precomputed novel-junction beta fits {args.novel_beta_fits}...")
    novel_fits = pd.read_csv(args.novel_beta_fits, sep='\t', index_col=0)

    all_fits = pd.concat([gtex_fits, novel_fits])
    all_fits['in_gtex_matrix'] = all_fits['in_gtex_matrix'].astype(bool)
    all_fits = all_fits[~all_fits.index.duplicated(keep='first')]

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

    idx_chrom_coord = all_fits.index.str.split('_')
    fits_chrom = idx_chrom_coord.map(lambda p: p[0])
    fits_start = idx_chrom_coord.map(lambda p: int(p[1]))
    fits_end = idx_chrom_coord.map(lambda p: int(p[2]))

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for region in regions:
            jxn_info_df_filtered = jxn_info_df[jxn_info_df['region'] == region].copy()

            reg_chrom, reg_coords = region.split(':')
            reg_start, reg_end = map(int, reg_coords.split('-'))
            fits_mask = (
                (fits_chrom == reg_chrom) &
                (fits_start >= reg_start) &
                (fits_end <= reg_end)
            )
            gtex_fits_filtered = all_fits[fits_mask].copy()

            futures.append(executor.submit(
                process_region, jxn_info_df_filtered, gtex_fits_filtered, region,
                args.sample_coverage_threshold, args.PSI_rescale_factor, args.phasing_threshold,
                annotated_junctions, report_outdir))

    results = []
    for future in concurrent.futures.as_completed(futures):
        try:
            res = future.result()
            if res is not None:
                results.append(res)
        except Exception as e:
            print(f"Error encountered: {e}")
            traceback.print_exc()

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
    padj_values = multipletests(numeric_p_values.dropna(), method='fdr_bh')[1]
    df.loc[mask & numeric_p_values.notna(), 'padj'] = padj_values

    df = df.sort_values(by=['sample', 'junction', 'phasing'])
    df.to_csv(args.outfile, sep='\t', index=False)
    print(f'Merged junctions written to {args.outfile}')


if __name__ == '__main__':
    main()
