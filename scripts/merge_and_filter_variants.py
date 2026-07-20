#!/usr/bin/env python

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2024.10.23

import argparse
import pandas as pd
import os
from math import ceil
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def parse_args():
    """ Parse command line arguments """
    
    parser = argparse.ArgumentParser(description="Remove variants from long-read RNA-seq data that are present in a large number of samples (i.e. from the same TEQUILA panel and/or sequencing batch).")
    parser.add_argument("--infiles", nargs="+", help="Input TSV file(s) containing the variants to be filtered (output from compile_variants.py).", required=True)
    parser.add_argument("--outprefix", help="Prefix for output files.", required=True)
    parser.add_argument("--title", default="Variant Counts", help="Title for the plot.")
    parser.add_argument("--include-multiallelic", action="store_true", help="Include multiallelic variants in the analysis. By default, they are excluded.")
    parser.add_argument("--filters", type=str, nargs="*", help="List of variant filters to apply. Example: --filters 'PASS' 'LowQual'")
    parser.add_argument("--min-DP", type=int, help="Minimum depth of coverage for variants")
    parser.add_argument("--max-DP", type=int, help="Maximum depth of coverage for variants")
    parser.add_argument("--min-DP-SNV", type=int, help="Minimum depth of coverage for SNVs. Overrides --min-DP")
    parser.add_argument("--max-DP-SNV", type=int, help="Maximum depth of coverage for SNVs. Overrides --max-DP")
    parser.add_argument("--min-DP-indel", type=int, help="Minimum depth of coverage for indels. Overrides --min-DP")
    parser.add_argument("--max-DP-indel", type=int, help="Maximum depth of coverage for indels. Overrides --max-DP")
    parser.add_argument("--num-callers-threshold", type=int, help="Minimum number of variant callers that must call a variant for it to be retained.")
    parser.add_argument("--num-callers-threshold-SNV", type=int, help="Minimum number of variant callers that must call a SNV for it to be retained. Overrides --num-callers-threshold")
    parser.add_argument("--num-callers-threshold-indel", type=int, help="Minimum number of variant callers that must call an indel for it to be retained. Overrides --num-callers-threshold")
    parser.add_argument("--sample-number-threshold", type=int, help="The maximum number of samples that a variant can be present in to be retained.")
    parser.add_argument("--plot", action="store_true", help="If set, generates a plot of the number of variants.")
    parser.add_argument("--plot-variant-type", type=str, nargs="+", choices=["SNV", "indel", "all"], default=["all"], \
                                            help="Type of variants to include in the plot (SNV, indel, or all). Can specify multiple types to output multiple plots, \
                                                                        e.g. --plot-variant-type 'SNV' 'indel'. By default, all types are included.")
    parser.add_argument("--plot-before-sample-number-filtering", action="store_true", help="If set, generates a plot of the number of variants before sample number filtering.")

    return parser.parse_args()

def plot_variant_counts(df, plot_variant_type, outfile, title):
    """ Generate a plot of the number of variants """

    # Make text editable in PDF/PS output
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # Safely copy the DataFrame to avoid modifying the original
    df = df.copy()

    # Create variant_id column
    df['variant_id'] = df['chrom'].astype(str) + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']

    # Create variant_caller column
    df['variant_caller'] = df['name'].apply(lambda x: 'longcallR' if x.endswith('_longcallR') else 
                                                     ('NanoTS' if x.endswith('_nanoTS') else 
                                                     ('Clair3-RNA' if x.endswith('_clair3_rna') else 
                                                      ('Deepvariant' if x.endswith('_deepvariant') else '')))
    )

    # Drop rows with empty variant_caller
    if '' in df['variant_caller'].values:
        print('WARNING: Empty variant_caller values found. Dropping these rows:')
        dropped_rows = df[df['variant_caller'] == '']
        print(dropped_rows)
        df = df[df['variant_caller'] != '']

    # Create sample_id column
    df['sample_id'] = df['name'].apply(lambda x: x.replace('_longcallR', '').replace('_nanoTS', '').replace('_clair3_rna', '').replace('_deepvariant', ''))

    # Get the full list of all sample IDs present in the data (before any variant-type subsetting)
    all_sample_ids = sorted(df['sample_id'].unique())

    # Define the color map
    color_map = {
        # 1 caller (red)
        'longcallR': '#c23637',
        'NanoTS': '#d95d5b',
        'Clair3-RNA': '#ea9a9c',
        'Deepvariant': '#f8c9c6',

        # 2 callers (blue)
        'longcallR,NanoTS': '#002b58',
        'longcallR,Clair3-RNA': '#0068a9',
        'longcallR,Deepvariant': '#4c8fca',
        'NanoTS,Clair3-RNA': '#91c4e9',
        'NanoTS,Deepvariant': '#a9cce6',
        'Clair3-RNA,Deepvariant': '#c1e3fa',

        # 3 callers (green)
        'longcallR,NanoTS,Clair3-RNA': '#3d892e',
        'longcallR,NanoTS,Deepvariant': '#57aa3e',
        'longcallR,Clair3-RNA,Deepvariant': '#95c36e',
        'NanoTS,Clair3-RNA,Deepvariant': '#d3e4be',

        # 4 callers (orange)
        'longcallR,NanoTS,Clair3-RNA,Deepvariant': '#f48f3e'
    }

    # Make a stacked barplot with subplots for each plot_variant_type
    num_rows = len(plot_variant_type)
    fig, axes = plt.subplots(nrows=num_rows, ncols=1, figsize=(12, 6), sharex=True)

    # Ensure axes is always iterable (handles num_rows == 1)
    if num_rows == 1:
        axes = [axes]

    # Count the counts for each combination of variant callers,
    # reindexing to include all samples even those with 0 variants.
    def get_plot_data(df, color_map, all_sample_ids):
        # Group by variant_id and sample_id to find which variant callers called each variant
        caller_order = {'longcallR': 0, 'NanoTS': 1, 'Clair3-RNA': 2, 'Deepvariant': 3}
        caller_combinations = df.groupby(['sample_id', 'variant_id'])['variant_caller']\
            .apply(lambda x: ','.join(sorted(x.unique(), key=lambda y: caller_order[y])))\
            .reset_index()
        # Count the number of variants for each combination of callers
        combination_counts = caller_combinations.groupby(['sample_id', 'variant_caller']).size().unstack(fill_value=0)
        # Reindex to include all samples (fills missing samples with 0)
        combination_counts = combination_counts.reindex(all_sample_ids, fill_value=0)
        # Reorder the columns
        col_order = [col for col in [
            # 4 callers
            'longcallR,NanoTS,Clair3-RNA,Deepvariant',

            # 3 callers
            'longcallR,NanoTS,Clair3-RNA',
            'longcallR,NanoTS,Deepvariant',
            'longcallR,Clair3-RNA,Deepvariant',
            'NanoTS,Clair3-RNA,Deepvariant',
            'Clair3-RNA,Deepvariant,longcallR',

            # 2 callers
            'longcallR,NanoTS',
            'longcallR,Clair3-RNA',
            'longcallR,Deepvariant',
            'NanoTS,Clair3-RNA',
            'NanoTS,Deepvariant',
            'Clair3-RNA,Deepvariant',

            # 1 caller
            'longcallR',
            'NanoTS'
            'Clair3-RNA',
            'Deepvariant',
        ] if col in combination_counts.columns]
        combination_counts = combination_counts[col_order]
        # Select colors
        colors = [color_map.get(col, '#000000') for col in combination_counts.columns]
        return combination_counts, colors
    
    # Plot counts
    i = 0
    if "all" in plot_variant_type:
        combination_counts, colors = get_plot_data(df, color_map, all_sample_ids)
        if not combination_counts.empty:
            combination_counts.plot(kind='bar', stacked=True, color=colors, ax=axes[i])
        else:
            print("WARNING: No variants found for 'all' category; skipping this subplot.")
        axes[i].set_ylabel('Number of Variants')
        i += 1
    if "SNV" in plot_variant_type:
        SNV_df = df[df['ref'].str.len() == 1]
        SNV_df = SNV_df[SNV_df['alt'].str.len() == 1]
        SNV_combination_counts, colors = get_plot_data(SNV_df, color_map, all_sample_ids)
        if not SNV_combination_counts.empty:
            SNV_combination_counts.plot(kind='bar', stacked=True, color=colors, ax=axes[i])
        else:
            print("WARNING: No SNVs found; skipping this subplot.")
        axes[i].set_ylabel('Number of SNVs')
        i += 1
    if "indel" in plot_variant_type:
        indel_df = df[(df['ref'].str.len() > 1) | (df['alt'].str.len() > 1)]
        indel_combination_counts, colors = get_plot_data(indel_df, color_map, all_sample_ids)
        if not indel_combination_counts.empty:
            indel_combination_counts.plot(kind='bar', stacked=True, color=colors, ax=axes[i])
        else:
            print("WARNING: No indels found; skipping this subplot.")
        axes[i].set_ylabel('Number of Indels')
    ticklabels = axes[i].xaxis.get_ticklabels()
    if ticklabels:
        axes[i].tick_params(axis='x', labelsize=ticklabels[0].get_size() * 0.5)
    axes[i].set_xlabel('Sample ID')
    for ax in axes:
        ax.legend().set_visible(False)

    # Add legend
    handles = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
    axes[-1].legend(handles=handles, title="Variant Callers", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Finish plotting
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()
    print(f"Plot saved to {outfile}")

def main():
    """ Main function """

    # Parse command line arguments
    args = parse_args()
    
    # Check that at least one input file was provided and that all input files are valid
    if len(args.infiles) == 0:
        raise ValueError("No input files provided.")
    for input_file in args.infiles:
        if not os.path.isfile(input_file):
            raise FileNotFoundError(f"Input file {input_file} not found.")

    # Read in the input TSV files
    merged_df = pd.DataFrame()
    for input_file in args.infiles:
        df = pd.read_csv(input_file, sep="\t")
        df.insert(0, "input_file", input_file)  # Add a column to indicate the input file
        df = df.loc[~df.duplicated()]
        if df.empty:
            print(f"WARNING: Input file {input_file} is empty.")
            continue
        merged_df = pd.concat([merged_df, df], ignore_index=True)

    # Remove multiallelic variants
    if not args.include_multiallelic:
        merged_df = merged_df[~merged_df['alt'].str.contains(',')]
        print(f"Removed multiallelic variants.")

    # Filter variants based on provided filters
    if args.filters:
        merged_df = merged_df[merged_df['filter'].isin(args.filters)]
        print(f"Filtered for variants with filters: {', '.join(args.filters)}")

    # Split into SNVs and indels
    SNV_df = merged_df[(merged_df['ref'].str.len() == 1) & (merged_df['alt'].str.len() == 1)]
    indel_df = merged_df[(merged_df['ref'].str.len() > 1) | (merged_df['alt'].str.len() > 1)]

    # Filter variants based on depth of coverage
    if args.min_DP or args.min_DP_SNV or args.max_DP or args.max_DP_SNV:
        initial_SNV_count = SNV_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        min_DP_SNV = args.min_DP_SNV if args.min_DP_SNV else args.min_DP if args.min_DP else None
        max_DP_SNV = args.max_DP_SNV if args.max_DP_SNV else args.max_DP if args.max_DP else None
        SNV_df = SNV_df[SNV_df['DP'] >= min_DP_SNV] if min_DP_SNV else SNV_df
        SNV_df = SNV_df[SNV_df['DP'] <= max_DP_SNV] if max_DP_SNV else SNV_df
        filtered_SNV_count = SNV_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        print(f"Filtered SNV variants for DP >= {min_DP_SNV} (from {initial_SNV_count} to {filtered_SNV_count})")
    if args.min_DP or args.min_DP_indel or args.max_DP or args.max_DP_indel:
        initial_indel_count = indel_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        min_DP_indel = args.min_DP_indel if args.min_DP_indel else args.min_DP if args.min_DP else None
        max_DP_indel = args.max_DP_indel if args.max_DP_indel else args.max_DP if args.max_DP else None
        indel_df = indel_df[indel_df['DP'] >= min_DP_indel] if min_DP_indel else indel_df
        indel_df = indel_df[indel_df['DP'] <= max_DP_indel] if max_DP_indel else indel_df
        filtered_indel_count = indel_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        print(f"Filtered indel variants for DP >= {min_DP_indel} (from {initial_indel_count} to {filtered_indel_count})")

    # Filter variants based on num_callers
    SNV_df['num_callers'] = SNV_df.groupby(['input_file', 'chrom', 'pos', 'ref', 'alt'])['name'].transform('nunique')
    indel_df['num_callers'] = indel_df.groupby(['input_file', 'chrom', 'pos', 'ref', 'alt'])['name'].transform('nunique')
    if args.num_callers_threshold or args.num_callers_threshold_SNV:
        initial_SNV_count = SNV_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        num_callers_threshold_SNV = args.num_callers_threshold_SNV if args.num_callers_threshold_SNV else args.num_callers_threshold
        SNV_df = SNV_df[SNV_df['num_callers'] >= num_callers_threshold_SNV]
        filtered_SNV_count = SNV_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        print(f"Filtered SNV variants for num_callers >= {num_callers_threshold_SNV} (from {initial_SNV_count} to {filtered_SNV_count})")
    if args.num_callers_threshold or args.num_callers_threshold_indel:
        initial_indel_count = indel_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        num_callers_threshold_indel = args.num_callers_threshold_indel if args.num_callers_threshold_indel else args.num_callers_threshold
        indel_df = indel_df[indel_df['num_callers'] >= num_callers_threshold_indel]
        filtered_indel_count = indel_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        print(f"Filtered indel variants for num_callers >= {num_callers_threshold_indel} (from {initial_indel_count} to {filtered_indel_count})")

    # Merge SNV and indel DataFrames
    merged_df = pd.concat([SNV_df, indel_df], ignore_index=True).sort_values(by=['input_file', 'chrom', 'pos'])

    # Filter variants based on sample_number_threshold
    if args.sample_number_threshold:
        merged_df['sample_count'] = merged_df.groupby(['chrom', 'pos', 'ref', 'alt'])['input_file'].transform('nunique')
        initial_count = merged_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        filtered_df = merged_df[merged_df['sample_count'] <= args.sample_number_threshold]
        filtered_count = filtered_df[['input_file', 'chrom', 'pos', 'ref', 'alt']].drop_duplicates().shape[0]
        print(f"Filtered variants for sample_count <= {args.sample_number_threshold} (from {initial_count} to {filtered_count})")
    else:
        filtered_df = merged_df.copy()
        print("No sample number filtering applied.")

    # Write the filtered DataFrame to file
    os.makedirs(os.path.dirname(args.outprefix), exist_ok=True)
    filtered_df.to_csv(args.outprefix + '.tsv', sep="\t", index=False)
    print(f"Saved to {args.outprefix + '.tsv'}.")

    # Generate a plot of the number of variants
    if args.plot:
        plot_variant_counts(filtered_df, args.plot_variant_type, args.outprefix + ".pdf", args.title)
        if args.plot_before_sample_number_filtering:
            plot_variant_counts(merged_df, args.plot_variant_type, args.outprefix + "_before_sample_number_filtering.pdf", args.title + " - Before Filtering by Sample Number")

if __name__ == "__main__":
    main()