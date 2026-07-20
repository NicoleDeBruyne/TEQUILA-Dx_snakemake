#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Optimized: 2025

# Extracts splice junctions from a BAM file mapping to specified regions of interest.

import os, argparse, warnings, pysam
import pandas as pd
from pandas.errors import PerformanceWarning
import concurrent.futures

warnings.filterwarnings('ignore', category=PerformanceWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# Pre-compiled CIGAR pattern used by every worker — compiled once at import time
import re
_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Identifies splice junctions (format: chrom_intronStart_intronEnd) mapping to a user-defined gene region')
    parser.add_argument('--mapping-file', required=True,
        help='Path to a TSV file with the gene regions of interest to analyze. Expected columns: '
             'sample_id, region (chr:start-end), gene, bamfile, hap1file (if available), hap2file (if available)')
    parser.add_argument('--outfile', required=True, help="Path to the output file")
    parser.add_argument('--threads', type=int, default=1,
        help='Number of threads to use for parallel processing. Default: 1')
    return parser.parse_args()

def get_splice_junctions_from_CIGAR(cigar, start, region_chrom, region_start, region_end):
    """Extract splice junctions from a single alignment.

    Accepts pre-parsed region fields to avoid redundant string splitting per alignment.
    Returns a list of junction strings: chrom_intronStart_intronEnd.
    """
    pos = start
    sj_list = []

    for length_str, operation in _CIGAR_RE.findall(cigar):
        length = int(length_str)

        if operation == 'P':
            warnings.warn(f"Padding operation 'P' encountered in CIGAR string: {cigar}")
        elif operation in 'SHI':
            pass  # no position change
        elif operation in 'MXD=':
            pos += length
        elif operation == 'N':
            junction_start = pos + 1
            junction_end = pos + length
            if junction_start < region_end and junction_end > region_start:
                sj_list.append(f"{region_chrom}_{junction_start}_{junction_end}")
            pos += length

    return sj_list

def extract_splice_junctions_from_BAM(bamfile, region):
    """Extract coordinates of all splice junctions from a BAM file for a region.

    Returns (sj_dict, alignment_count).
    """
    # Parse region once
    region_chrom, region_boundaries = region.split(':')
    region_start, region_end = map(int, region_boundaries.split('-'))

    sj_dict = {}
    alignment_count = 0

    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        for alignment in bam.fetch(region=region):
            if alignment.is_secondary:
                continue
            alignment_count += 1
            cigar = alignment.cigarstring
            if cigar is None:
                continue
            for junction in get_splice_junctions_from_CIGAR(
                    cigar, alignment.reference_start, region_chrom, region_start, region_end):
                if junction in sj_dict:
                    sj_dict[junction] += 1
                else:
                    sj_dict[junction] = 1

    return sj_dict, alignment_count

def _build_df(sj_dict, sample_id, phasing, region, gene, gene_alignment_count):
    """Helper: build a junction DataFrame from a dict."""
    if not sj_dict:
        return pd.DataFrame(columns=['sample', 'phasing', 'region', 'gene',
                                     'gene_alignment_count', 'junction', 'jxn_alignment_count'])
    df = pd.DataFrame(
        {'junction': list(sj_dict.keys()), 'jxn_alignment_count': list(sj_dict.values())}
    )
    df['sample'] = sample_id
    df['phasing'] = phasing
    df['region'] = region
    df['gene'] = gene
    df['gene_alignment_count'] = gene_alignment_count
    return df[['sample', 'phasing', 'region', 'gene', 'gene_alignment_count', 'junction', 'jxn_alignment_count']]

def process_region(sample_id, gene, region, bamfile, report_outdir, hap1file=None, hap2file=None):
    """Process one region of interest and return a combined junction DataFrame."""

    report = os.path.join(report_outdir,
        f"{sample_id}_{gene}_{region.replace(':', '_').replace('-', '_')}_report.tsv")

    with open(report, 'w') as report_file:
        report_file.write(f"Obtaining splice junction information for {sample_id} over {gene} ({region})\n")

        # Bulk BAM
        report_file.write("Extracting splice junctions from the sample of interest...\n")
        sj_dict, gene_alignment_count = extract_splice_junctions_from_BAM(bamfile, region)
        if not sj_dict:
            report_file.write("    No junctions found in sample of interest. Exiting...\n")
            return None

        sj_df = _build_df(sj_dict, sample_id, 'bulk', region, gene, gene_alignment_count)
        report_file.write(
            f"    Processed {gene_alignment_count} primary/supplementary alignments "
            f"and found {len(sj_dict)} junctions.\n")

        bulk_junctions = set(sj_dict.keys())
        del sj_dict

        if hap1file is None or hap2file is None:
            return sj_df

        # Haplotype-specific BAMs
        report_file.write("Processing haplotype-specific data...\n")

        hap1_dict, hap1_count = extract_splice_junctions_from_BAM(hap1file, region)
        report_file.write(
            f"    Haplotype 1: Processed {hap1_count} primary/supplementary alignments "
            f"and found {len(hap1_dict)} junctions.\n")
        hap1_df = _build_df(hap1_dict, sample_id, 'hap1', region, gene, hap1_count)
        # Fill in zero-count junctions seen in bulk but not in hap1
        missing_hap1 = bulk_junctions - hap1_dict.keys()
        if missing_hap1:
            hap1_df = pd.concat([hap1_df, pd.DataFrame({
                'sample': sample_id, 'phasing': 'hap1', 'region': region,
                'gene': gene, 'gene_alignment_count': hap1_count,
                'junction': list(missing_hap1), 'jxn_alignment_count': 0
            })], ignore_index=True)
        del hap1_dict

        hap2_dict, hap2_count = extract_splice_junctions_from_BAM(hap2file, region)
        report_file.write(
            f"    Haplotype 2: Processed {hap2_count} primary/supplementary alignments "
            f"and found {len(hap2_dict)} junctions.\n")
        hap2_df = _build_df(hap2_dict, sample_id, 'hap2', region, gene, hap2_count)
        missing_hap2 = bulk_junctions - hap2_dict.keys()
        if missing_hap2:
            hap2_df = pd.concat([hap2_df, pd.DataFrame({
                'sample': sample_id, 'phasing': 'hap2', 'region': region,
                'gene': gene, 'gene_alignment_count': hap2_count,
                'junction': list(missing_hap2), 'jxn_alignment_count': 0
            })], ignore_index=True)
        del hap2_dict

    return pd.concat([sj_df, hap1_df, hap2_df], ignore_index=True)

def main():
    """Main script."""

    print("\n\n\n******************************************************************************************")
    print("Getting splice junction counts...")
    print("******************************************************************************************\n")

    args = parse_args()

    outdir = os.path.dirname(args.outfile)
    os.makedirs(outdir, exist_ok=True)
    report_outdir = os.path.join(outdir, 'reports')
    os.makedirs(report_outdir, exist_ok=True)

    print(f"\nBegin processing regions of interest from mapping file {args.mapping_file}...")
    mapping_df = pd.read_csv(args.mapping_file, sep='\t')

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(
                process_region,
                row.iloc[0],           # sample_id
                row.iloc[2],           # gene
                row.iloc[1],           # region
                row.iloc[3],           # bamfile
                report_outdir,
                row.iloc[4] if pd.notna(row.iloc[4]) else None,
                row.iloc[5] if pd.notna(row.iloc[5]) else None,
            )
            for _, row in mapping_df.iterrows()
        ]

        results = []
        for future in concurrent.futures.as_completed(futures):
            try:
                res = future.result()
                if res is not None:
                    results.append(res)
            except Exception as e:
                print(f"Error encountered: {e}")

    combined_df = pd.concat(results, ignore_index=True)
    combined_df.to_csv(args.outfile, sep='\t', index=False)
    print(f"\nSplice junction information saved to {args.outfile}\n")

if __name__ == '__main__':
    main()
