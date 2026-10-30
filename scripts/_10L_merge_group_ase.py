# 10L - merge ASE outliers across the group and drop genes with ASE in too many samples

import argparse
import pandas as pd
import os
from palette import SELECT_ORDER, use_fonts
from matplotlib import rcParams
from plot_utils import plot_outlier_counts
use_fonts()
rcParams['pdf.fonttype'] = 42

def parse_args():
    
    parser = argparse.ArgumentParser(description="Remove genes with outlier ASE from long-read RNA-seq data that are present in a large number of samples (i.e. from the same TEQUILA panel and/or sequencing batch).")
    parser.add_argument("--infiles", nargs="+", help="Input TSV files containing the outlier genes to be filtered (output from run_ase_analysis.py).", required=True)
    parser.add_argument("--outprefix", help="Prefix for output files.", required=True)
    parser.add_argument("--title", default="Number of Genes with Allele-specific Expression by Sample", help="Title for the plot.")
    parser.add_argument('--min-haplotype-ratio', default=0, type=float, help='Exclude genes with haplotype ratio < min-haplotype-ratio, for example, if you believe these are false positives.')
    parser.add_argument('--minor-haplotype-frequency-threshold', default=0.4, type=float, help='A gene is considered an outlier when its minor-haplotype fraction (ratio) is below this value. (default: 0.4)')
    parser.add_argument('--padj-threshold', default=0.05, type=float, help='Maximum adjusted p-value to be considered an outlier. (default: 0.05)')
    parser.add_argument("--sample-number-threshold", type=int, help="The maximum number of samples that a gene can be an outlier in to be retained.")
    parser.add_argument("--plot", action="store_true", help="If set, generates a plot of the number of genes with outlier ASE according to the haplotype-ratio and padj thresholds before and after sample number filtering.")
    
    return parser.parse_args()


def main():

    args = parse_args()
    
    if len(args.infiles) == 0:
        raise ValueError("No input files provided.")
    for input_file in args.infiles:
        if not os.path.isfile(input_file):
            raise FileNotFoundError(f"Input file {input_file} not found.")

    print(f"Reading in {len(args.infiles)} input files...")
    merged_df = pd.DataFrame()
    samples = set()
    for input_file in args.infiles:
        df = pd.read_csv(input_file, sep="\t")
        if df.empty:
            print(f"WARNING: Input file {input_file} is empty.")
            continue
        samples.update(df['sample'].unique())
        df = df[(df['ratio'] >= args.min_haplotype_ratio) & (df['ratio'] < args.minor_haplotype_frequency_threshold) & (df['padj'] <= args.padj_threshold)]
        df.insert(0, "input_file", input_file)
        merged_df = pd.concat([merged_df, df], ignore_index=True)

    if not 'gene' in merged_df.columns:
        raise ValueError("Input files must contain 'gene' column.")

    # recurrence filter
    merged_df['sample_count'] = merged_df.groupby('gene')['sample'].transform('nunique')

    if args.sample_number_threshold:
        filtered_df = merged_df[merged_df['sample_count'] <= args.sample_number_threshold]
    else:
        filtered_df = merged_df.copy()

    os.makedirs(os.path.dirname(args.outprefix), exist_ok=True)
    filtered_df.to_csv(args.outprefix + ".tsv", sep="\t", index=False)
    print(f"Merged ASE results and saved to {args.outprefix + '.tsv'}.")

    if args.plot:
        legends = [
            (SELECT_ORDER[0], f'Genes with minor haplotype fraction < {args.minor_haplotype_frequency_threshold} and padj ≤ {args.padj_threshold}'), \
            (SELECT_ORDER[1], f'Genes with allelic imbalance in ≤ {args.sample_number_threshold} samples')
        ]
        plot_outlier_counts(samples, [merged_df, filtered_df], legends, args.title, f"{args.outprefix}_counts.pdf")

if __name__ == "__main__":
    main()