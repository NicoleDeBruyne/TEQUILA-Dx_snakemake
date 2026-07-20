#!/usr/bin/env python

# Author: Nicole DeBruyne (Lin Lab)
# Original Date: 2024.10.23

import argparse
import pandas as pd
import numpy as np
from itertools import combinations
from collections import defaultdict
from scipy.stats import beta, betabinom
import os
from math import ceil
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

def parse_args():
    """ Parse command line arguments """
    parser = argparse.ArgumentParser(description="Remove outlier junctions from long-read RNA-seq data.")
    parser.add_argument("--infiles", nargs="+", required=True, help="Input TSV files (merged junctions).")
    parser.add_argument("--outprefix", required=True, help="Prefix for output files.")
    parser.add_argument('--jxn-coverage-threshold', default=20, type=int, help='Minimum junction coverage to be considered an outlier.')
    parser.add_argument('--padj-threshold', default=0.05, type=float, help='Maximum adjusted p-value to be considered an outlier.')
    parser.add_argument('--delta-PSI-threshold', default=0.1, type=float, help='Minimum abs(delta_PSI) to be considered an outlier.')
    parser.add_argument("--delta-PSI-direction", choices=["positive", "negative"], help="If set, only keep junctions with delta PSI in specified direction.")
    parser.add_argument("--event-types", nargs="+", choices=["exon_skipping", "exon_inclusion", "alt_ss1", "alt_ss2", "single"], help="Event types to include when filtering junctions.")
    parser.add_argument("--include-upreg-annotation-types", nargs="+", help="Annotation types to include when filtering upregulated junctions. Excludes all others.")
    parser.add_argument("--include-downreg-annotation-types", nargs="+", help="Annotation types to include when filtering downregulated junctions. Excludes all others.")
    parser.add_argument("--exclude-upreg-annotation-types", nargs="+", help="Annotation types to exclude when filtering upregulated junctions. Includes all others.")
    parser.add_argument("--exclude-downreg-annotation-types", nargs="+", help="Annotation types to exclude when filtering downregulated junctions. Includes all others.")
    parser.add_argument("--sample-number-threshold", type=int, help="Max number of samples that a junction can be an outlier in to be retained")
    parser.add_argument("--filter-by-cohort-IQR", action="store_true", help="Require junctions to also be outliers relative to the cohort based on IQR.")
    parser.add_argument("--fit-beta-dist-on-cohort", action="store_true", help="Fit beta distributions on cohort PSI values.")
    parser.add_argument("--compare-PSI-to-cohort-median", action="store_true", help="Compare each sample's PSI to cohort median. Ignored if --fit-beta-dist-on-cohort is set.")
    parser.add_argument("--delta-PSI-vs-cohort-threshold", type=float, help="Min abs(delta_PSI_vs_cohort) to filter when using cohort statistics (beta distribution or median).")
    parser.add_argument("--plot", action="store_true", help="If set, generates plots.")
    parser.add_argument("--title", default="Number of Genes with Outlier Junction Hits by Sample", help="Title for the plot.")
    return parser.parse_args()

##################################################
# Helper functions
##################################################

def define_events(df):
    df[['chr', 'ss1', 'ss2']] = df['junction'].str.split('_', expand=True)
    df['ss1'] = df['ss1'].astype(int)
    df['ss2'] = df['ss2'].astype(int)

    junction_events = defaultdict(list)

    for (sample, phasing, chrom), group in df.groupby(['sample', 'phasing', 'chr']):
        group = group.sort_values(['ss1', 'ss2']).reset_index(drop=True)
        junctions = group.to_dict('records')
        event_dict = defaultdict(list)

        for comb in combinations(junctions, 3):
            j0, j1, j2 = comb
            if j0['ss1'] == j1['ss1'] and j1['ss2'] == j2['ss2']:
                if j0['delta_PSI'] < 0 and j1['delta_PSI'] > 0 and j2['delta_PSI'] < 0:
                    for j in [j0, j1, j2]:
                        event_dict[j['junction']].append('exon_skipping')
                elif j0['delta_PSI'] > 0 and j1['delta_PSI'] < 0 and j2['delta_PSI'] > 0:
                    for j in [j0, j1, j2]:
                        event_dict[j['junction']].append('exon_inclusion')

        for comb in combinations(junctions, 2):
            j0, j1 = comb
            if j0['ss1'] == j1['ss1']:
                for j in [j0, j1]:
                    event_dict[j['junction']].append('alt_ss2')
            if j0['ss2'] == j1['ss2']:
                for j in [j0, j1]:
                    event_dict[j['junction']].append('alt_ss1')

        for j in junctions:
            junc = j['junction']
            events = event_dict[junc]
            p1 = [e for e in events if e in ['exon_skipping', 'exon_inclusion']]
            p2 = [e for e in events if e in ['alt_ss1', 'alt_ss2']]
            final_events = p1 if p1 else p2 if p2 else ['single']
            junction_events[(j['sample'], j['phasing'], junc)] = ';'.join(final_events)

    df['event'] = df.set_index(['sample', 'phasing', 'junction']).index.map(junction_events.get).fillna('single')
    df = df.drop(columns=['chr', 'ss1', 'ss2'])
    return df

def fit_beta_dist(x, tol, n_threshold):
    x = pd.to_numeric(x, errors='coerce')
    x = x[~np.isnan(x)]
    if len(x) < n_threshold:
        return (len(x), "low_n", "low_n", "low_n")
    if x.var() < tol:
        return (len(x), x.mean() / tol, (1 - x.mean()) / tol, x.mean())
    try:
        alpha_value, beta_value = beta.fit(x, floc = 0, fscale = 1)[0:2]
        expected_PSI = alpha_value / (alpha_value + beta_value)
        return (len(x), alpha_value, beta_value, expected_PSI)
    except:
        return (len(x), "error", "error", "error")

def beta_binomial_test(x, n, alpha_value, beta_value):
    if any([np.isnan(val) for val in [x, n, alpha_value, beta_value]]):
        return "n/a"
    x = round(x)
    n = round(n)
    lte_x = betabinom.cdf(x, n, alpha_value, beta_value)
    gte_x = betabinom.cdf(n - x, n, beta_value, alpha_value)
    p_value = np.clip(2 * min(lte_x, gte_x), 0, 1)
    if np.isnan(p_value):
        return "error"
    return p_value

def compute_cohort_statistics(df, args):
    print(f"    Computing cohort statistics for {len(df['junction'].unique())} junctions...")
    jxns_interest = df['junction'].unique()
    dfs_bulk = []

    for input_file in args.infiles:
        df_bulk = pd.read_csv(input_file, sep="\t")
        if df_bulk.empty:
            continue
        df_bulk[['chr', 'ss1', 'ss2']] = df_bulk['junction'].str.split('_', n=2, expand=True)
        missing_jxns = set(jxns_interest) - set(df_bulk['junction'])
        if missing_jxns:
            df_missing = pd.DataFrame({'junction': list(missing_jxns)})
            df_missing['phasing'] = 'bulk'
            df_missing[['chr', 'ss1', 'ss2']] = df_missing['junction'].str.split('_', n=2, expand=True)
            ss1_sums = df_bulk.groupby('ss1')['ss1_coverage'].sum()
            ss2_sums = df_bulk.groupby('ss2')['ss2_coverage'].sum()
            df_missing['ss1_coverage'] = df_missing['ss1'].map(ss1_sums).fillna(0)
            df_missing['ss2_coverage'] = df_missing['ss2'].map(ss2_sums).fillna(0)
            df_missing['jxn_coverage'] = df_missing['ss1_coverage'] + df_missing['ss2_coverage']
            df_missing['sample_PSI'] = 0
            df_missing['rescaled_sample_PSI'] = 0.001
            df_bulk = pd.concat([df_bulk, df_missing], ignore_index=True)
        df_bulk = df_bulk[(df_bulk['junction'].isin(jxns_interest)) &
                          (df_bulk['phasing'] == 'bulk') &
                          (df_bulk['jxn_coverage'] >= 20)]
        df_bulk.insert(0, "input_file", input_file)
        df_bulk = df_bulk[['input_file', 'junction', 'sample_PSI', 'rescaled_sample_PSI']]
        df_bulk["rescaled_sample_PSI"] = pd.to_numeric(df_bulk["rescaled_sample_PSI"], errors="coerce")
        dfs_bulk.append(df_bulk)

    if len(dfs_bulk) == 0:
        print("    No bulk samples with sufficient coverage found; skipping cohort statistics computation.")
        return df

    all_bulk_df = pd.concat(dfs_bulk, ignore_index=True)
    all_bulk_df = all_bulk_df.pivot_table(index='junction', columns='input_file', values='rescaled_sample_PSI', aggfunc='first')
    
    if args.filter_by_cohort_IQR:
        print("    Computing cohort IQR statistics...")
        all_bulk_numeric = all_bulk_df.apply(pd.to_numeric, errors="coerce")
        q1 = all_bulk_numeric.quantile(0.25, axis=1)
        q3 = all_bulk_numeric.quantile(0.75, axis=1)
        median = all_bulk_numeric.median(axis=1)
        fence_df = pd.DataFrame({"cohort_q1": q1, "cohort_median": median, "cohort_q3": q3})
        df = df.merge(fence_df, left_on="junction", right_index=True, how="left")

    elif args.fit_beta_dist_on_cohort:
        print("    Fitting beta distributions on cohort bulk PSI values...")
        all_bulk_df[['cohort_n', 'cohort_alpha', 'cohort_beta', 'cohort_expected_PSI']] = all_bulk_df.apply(
            lambda row: pd.Series(fit_beta_dist(np.array(row), tol=1e-3, n_threshold=30)), axis=1
        )
        df = df.merge(all_bulk_df[['cohort_n', 'cohort_alpha', 'cohort_beta', 'cohort_expected_PSI']],
                      left_on='junction', right_index=True, how='left')
        df['delta_PSI_vs_cohort'] = np.where(
            (np.isnan(pd.to_numeric(df['rescaled_sample_PSI'], errors='coerce')) |
             np.isnan(pd.to_numeric(df['cohort_expected_PSI'], errors='coerce'))),
            "n/a",
            pd.to_numeric(df['cohort_expected_PSI'], errors='coerce') - pd.to_numeric(df['rescaled_sample_PSI'], errors='coerce')
        )

        print("    Calculating cohort-level p-values...")
        df['cohort_p_value'] = df.apply(lambda row: beta_binomial_test(
            ceil(pd.to_numeric(row['jxn_alignment_count'], errors='coerce')),
            ceil(pd.to_numeric(row['jxn_coverage'], errors='coerce')),
            pd.to_numeric(row['cohort_alpha'], errors='coerce'),
            pd.to_numeric(row['cohort_beta'], errors='coerce')
        ), axis=1)

    elif args.compare_PSI_to_cohort_median:
        print("    Calculating cohort median PSI values...")
        all_bulk_df['cohort_median_PSI'] = all_bulk_df.median(axis=1, skipna=True)
        df = df.merge(all_bulk_df[['cohort_median_PSI']], left_on='junction', right_index=True, how='left')
        df['delta_PSI_vs_cohort'] = np.where(
            (np.isnan(pd.to_numeric(df['sample_PSI'], errors='coerce')) |
             np.isnan(pd.to_numeric(df['cohort_median_PSI'], errors='coerce'))),
            "n/a",
            pd.to_numeric(df['sample_PSI'], errors='coerce') - pd.to_numeric(df['cohort_median_PSI'], errors='coerce')
        )

    return df

##################################################
# Plotting functions
##################################################

def plot_outlier_counts(samples, dfs, legends, suptitle, outfile):
    """Generate bar plots and box plots."""

    # Count unique genes per sample for each dataset
    counts_list = []
    for df in dfs:
        counts_list.append(df.groupby('sample')['gene'].nunique().sort_index())

    # Ensure all samples are represented in each series in the same order
    sample_order = sorted(samples)
    for i, counts in enumerate(counts_list):
        counts_list[i] = counts.reindex(sample_order, fill_value=0)

    # ---------------- Bar plots ----------------
    barplot_file = os.path.splitext(outfile)[0] + "_barplot.pdf"
    fig, axes = plt.subplots(len(dfs), 1, figsize=(max(16, len(sample_order)*0.2), max(6, len(dfs)*4)), sharex=True)
    if len(dfs) == 1:
        axes = [axes]  # ensure iterable

    for ax, counts, (color, label) in zip(axes, counts_list, legends):
        ax.bar(counts.index, counts.values, color=color)
        ax.set_title(f'Average={counts.mean():.2f}, Median={counts.median()}, Range={counts.min()}–{counts.max()}')
        ax.set_ylabel('Count')
        ax.tick_params(axis='x', rotation=90)
        ax.legend(handles=[mpatches.Patch(color=color, label=label)],
                  loc='center left', bbox_to_anchor=(1.0, 0.5))

    axes[-1].set_xlabel('Sample')
    fig.suptitle(suptitle, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(barplot_file)
    plt.close()
    print(f"Bar plot saved to {barplot_file}")

    # ---------------- Box plot ----------------
    boxplot_file = os.path.splitext(outfile)[0] + "_boxplot.pdf"
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Data for boxplot: one box per dataset
    data = [counts.values for counts in counts_list]
    colors = [c for c, _ in legends]

    # Draw empty boxes with black edges and black median
    box = ax.boxplot(
        data,
        patch_artist=False,
        medianprops=dict(color='black'),
        boxprops=dict(color='black'),
        whiskerprops=dict(color='black'),
        capprops=dict(color='black'),
        showfliers=False
    )

    # Plot points with jitter and dataset colors
    jitter_strength = 0.08
    for i, counts in enumerate(data):
        x = np.full(len(counts), i + 1) + np.random.uniform(-jitter_strength, jitter_strength, len(counts))
        ax.scatter(x, counts, color=colors[i], alpha=0.7, s=20)

    # Clean x-axis (no labels, just spacing)
    ax.set_xticks(range(1, len(dfs) + 1))
    ax.set_xticklabels([])

    ax.set_ylabel('Number of genes')
    ax.set_title(suptitle)

    # Legend with stats included, with stats on new line
    legend_handles = []
    for counts, color, (_, label) in zip(counts_list, colors, legends):
        stats_text = f"{label}\nAvg={counts.mean():.2f}, Med={counts.median()}, Range={counts.min()}–{counts.max()}"
        legend_handles.append(mpatches.Patch(color=color, label=stats_text))
    ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1.05, 0.5))

    plt.tight_layout()
    plt.savefig(boxplot_file, bbox_inches='tight')
    plt.close()
    print(f"Box plot saved to {boxplot_file}")

def plot_outlier_types(df, suptitle, outfile):
    """ Make a pie chart of outlier types: phasing type, annotation type, event type """

    annotation_colors = {
        'unannotated': '#d95d5b',
        'protein_coding': '#4c8fca',
        'protein_coding_CDS_not_defined': '#57aa3e',
        'nonsense_mediated_decay': '#f48f3e',
        'retained_intron': '#b271ab',
        'lncRNA': '#42b4b5',
        'other': '#a0a0a0'
    }
    event_alphas = {
        'exon_skipping/inclusion': 1.0,
        'alt_ss': 0.6,
        'single': 0.3,
        'other': 0.3
    }

    phasing_groups = {
        'bulk': df[df['phasing'] == 'bulk'],
        'haplotype-specific': df[df['phasing'].isin(['hap1','hap2'])],
        'other': df[~df['phasing'].isin(['bulk','hap1','hap2'])]
    }
    annotation_groups = {
        'unannotated': df[df['annotation'] == 'unannotated'],
        'protein_coding': df[df['annotation'] == 'annotated:protein_coding'],
        'protein_coding_CDS_not_defined': df[df['annotation'] == 'annotated:protein_coding_CDS_not_defined'],
        'nonsense_mediated_decay': df[df['annotation'] == 'annotated:nonsense_mediated_decay'],
        'retained_intron': df[df['annotation'] == 'annotated:retained_intron'],
        'lncRNA': df[df['annotation'] == 'annotated:lncRNA'],
        'other': df[
            (df['annotation'] == 'other') |
            (~df['annotation'].isin(['unannotated','annotated:protein_coding','annotated:protein_coding_CDS_not_defined',
                'annotated:nonsense_mediated_decay','annotated:retained_intron','lncRNA']))
        ]
    }
    event_groups = {
        'exon_skipping/inclusion': df[df['event'].str.contains('exon_skipping|exon_inclusion')],
        'alt_ss': df[df['event'].str.contains('alt_ss1|alt_ss2')],
        'single': df[df['event'] == 'single'],
        'other': df[
            (df['event'] != 'single') &
            (~df['event'].str.contains('exon_skipping|exon_inclusion|alt_ss1|alt_ss2'))
        ]
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, (phasing_label, ph_df) in zip(axes, phasing_groups.items()):
        sizes, colors = [], []
        for ann_label, ann_df in annotation_groups.items():
            sub_df = ph_df[ph_df.index.isin(ann_df.index)]
            for evt_label, evt_df in event_groups.items():
                count = len(sub_df[sub_df.index.isin(evt_df.index)])
                if count > 0:
                    sizes.append(count)
                    base_color = annotation_colors[ann_label]
                    rgba = list(plt.matplotlib.colors.to_rgba(base_color))
                    rgba[-1] = event_alphas[evt_label]
                    colors.append(rgba)
        if sizes:
            ax.pie(sizes, labels=None, colors=colors, startangle=90, wedgeprops=dict(edgecolor='white'))
        ax.set(aspect="equal")
        ax.set_title(f"{phasing_label.capitalize()} (n={sum(sizes)})")

    used_annotations = set()
    used_events = set()
    for ax, (phasing_label, ph_df) in zip(axes, phasing_groups.items()):
        for ann_label, ann_df in annotation_groups.items():
            sub_df = ph_df[ph_df.index.isin(ann_df.index)]
            for evt_label, evt_df in event_groups.items():
                count = len(sub_df[sub_df.index.isin(evt_df.index)])
                if count > 0:
                    used_annotations.add(ann_label)
                    used_events.add(evt_label)
    ann_patches = [mpatches.Patch(color=annotation_colors[a], label=a) for a in annotation_colors if a in used_annotations]
    evt_patches = [mpatches.Patch(facecolor='gray', alpha=event_alphas[e], label=e) for e in event_alphas if e in used_events]
    leg1 = fig.legend(handles=ann_patches, title="Annotation", loc="center right", bbox_to_anchor=(1.15, 0.7))
    fig.add_artist(leg1)
    fig.legend(handles=evt_patches, title="Event type", loc="center right", bbox_to_anchor=(1.15, 0.3))

    plt.suptitle(suptitle)
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight")
    plt.close()
    print(f"Pie charts saved to {outfile}")

def main():
    args = parse_args()

    if args.include_upreg_annotation_types and args.exclude_upreg_annotation_types:
        raise ValueError("--include-upreg-annotation-types and --exclude-upreg-annotation-types are mutually exclusive")
    if args.include_downreg_annotation_types and args.exclude_downreg_annotation_types:
        raise ValueError( "--include-downreg-annotation-types and --exclude-downreg-annotation-types are mutually exclusive")

    # Print summary
    print(f"Preparing to filter for junctions with coverage >= {args.jxn_coverage_threshold}, padj <= {args.padj_threshold}, and abs(delta_PSI) >= {args.delta_PSI_threshold}.")
    if args.delta_PSI_direction:
        print(f"Will keep an additional file with only outlier junctions with delta_PSI direction: {args.delta_PSI_direction}")
    if args.event_types:
        print(f"Will keep an additional file with only outlier junctions with event types: {', '.join(args.event_types)}")
    if args.include_upreg_annotation_types or args.exclude_upreg_annotation_types:
        if args.delta_PSI_direction and args.delta_PSI_direciton=="negative":
            print(f"User specified negative delta_PSI_direction, so upregulated junctions cannot be filtered by annotation type.")
        elif args.include_upreg_annotation_types:
            print(f"Will keep an additional file with upregulated outlier junctions only with annotation types: {', '.join(args.include_upreg_annotation_types)}")
        else:
            print(f"Will keep an additional file excluding upregulated outlier junctions with annotation types: {', '.join(args.exclude_upreg_annotation_types)}")
    if args.include_downreg_annotation_types or args.exclude_downreg_annotation_types:
        if args.delta_PSI_direction and args.delta_PSI_direciton=="positive":
            print(f"User specified positive delta_PSI_direction, so downregulation junctions cannot be filtered by annotation type.")
        elif args.include_upreg_annotation_types or args.exclude_upreg_annotation_types:
            if args.include_downreg_annotation_types:
                print(f"    and will keep only downregulated outlier junctions with annotation types: {', '.join(args.include_downreg_annotation_types)}")
            else:
                print(f"    and excluding downregulated outlier junctions with annotation types: {', '.join(args.exclude_downreg_annotation_types)}")
        else:
            if args.include_downreg_annotation_types:
                print(f"Will keep an additional file with only downregulated outlier junctions with annotation types: {', '.join(args.include_downreg_annotation_types)}")
            else:
                print(f"Will keep an additional file excluding downregulated outlier junctions with annotation types: {', '.join(args.exclude_downreg_annotation_types)}")
    if (args.filter_by_cohort_IQR or ((args.fit_beta_dist_on_cohort or args.compare_PSI_to_cohort_median) and args.delta_PSI_vs_cohort_threshold)):
        if args.filter_by_cohort_IQR:
            print(f"Will keep an additional file with only junctions which are IQR outliers.")
        elif args.fit_beta_dist_on_cohort:
            print(f"Will fit beta distributions on cohort PSI values and keep an additional file with only junctions where abs(delta_PSI_vs_cohort) >= {args.delta_PSI_vs_cohort_threshold}.")
        elif args.compare_PSI_to_cohort_median:
            print(f"Will compore PSI value to the cohort median and keep an additional file with only junctions where abs(delta_PSI_vs_cohort) >= {args.delta_PSI_vs_cohort_threshold}.")
    if args.sample_number_threshold:
        print(f"Will keep an additional file with only junctions which are outliers in <= {args.sample_number_threshold} samples.")

    # Input checks
    if len(args.infiles) == 0:
        raise ValueError("No input files provided.")
    for input_file in args.infiles:
        if not os.path.isfile(input_file):
            raise FileNotFoundError(f"Input file {input_file} not found.")

    # Stage 1: read and apply jxn_coverage, padj, deltaPSI
    print(f"\nReading in {len(args.infiles)} input files...")
    dfs = []
    samples = set()
    for input_file in args.infiles:
        df = pd.read_csv(input_file, sep="\t", usecols=['sample', 'phasing', 'gene', 'junction', 'jxn_alignment_count', 'jxn_coverage', 
                                                            'sample_PSI', 'rescaled_sample_PSI', 'padj', 'delta_PSI', 'annotation'])
        if df.empty:
            print(f"WARNING: Input file {input_file} is empty.")
            continue
        samples.add(df['sample'].iloc[0])
        mask = (df['jxn_coverage'] >= args.jxn_coverage_threshold) & \
            (df['padj'] <= args.padj_threshold) & \
            (df['delta_PSI'].abs() >= args.delta_PSI_threshold)
        df = df.loc[mask]
        if df.empty:
            continue
        df.insert(0, "input_file", input_file)
        dfs.append(df)
    if len(dfs) == 0:
        raise RuntimeError("No data after primary filtering.")
    outlier_df = pd.concat(dfs, ignore_index=True)
    jxns = outlier_df['junction'].unique()
    print(f"Filtered for junctions with junction coverage >= {args.jxn_coverage_threshold}, padj <= {args.padj_threshold}, and abs(delta PSI) >= {args.delta_PSI_threshold}:\n"
        f"    {len(outlier_df)} total hits ({len(outlier_df['junction'].unique())} unique junctions) in {len(outlier_df['sample'].unique())} samples."
    )

    # Define alternative splicing events
    outlier_df = define_events(outlier_df.copy())

    # Add sample_count column
    outlier_df['sample_count'] = outlier_df.groupby('junction')['sample'].transform('nunique')

    # Compute cohort statistics if requested
    if len(samples) >= 8 and (args.filter_by_cohort_IQR or ((args.fit_beta_dist_on_cohort or args.compare_PSI_to_cohort_median) and args.delta_PSI_vs_cohort_threshold)):
        print(f"\nComputing cohort-level statistics...")
        outlier_df = compute_cohort_statistics(outlier_df, args)
    elif args.filter_by_cohort_IQR or ((args.fit_beta_dist_on_cohort or args.compare_PSI_to_cohort_median) and args.delta_PSI_vs_cohort_threshold):
        if len(samples) < 8:
            print("WARNING: Fewer than 8 samples were provided. Skipping cohort statistics computation.")
            args.filter_by_cohort_IQR = False
            args.fit_beta_dist_on_cohort = False
            args.compare_PSI_to_cohort_median = False
            args.delta_PSI_vs_cohort_threshold = None

    # Save file
    os.makedirs(os.path.dirname(args.outprefix), exist_ok=True)
    outfile_stage1 = f"{args.outprefix}_{args.jxn_coverage_threshold}jxncov_{args.padj_threshold}padj_{args.delta_PSI_threshold}deltaPSI.tsv"
    outlier_df.to_csv(outfile_stage1, sep="\t", index=False)

    # Stage 2: delta_PSI_direction
    if args.delta_PSI_direction:
        df_stage2 = outlier_df.copy()
        if args.delta_PSI_direction == "positive":
            df_stage2 = df_stage2[df_stage2['delta_PSI'] > 0]
        else:
            df_stage2 = df_stage2[df_stage2['delta_PSI'] < 0]
        outfile_stage2 = outfile_stage1.replace(".tsv", f"_{args.delta_PSI_direction}.tsv")
        df_stage2.to_csv(outfile_stage2, sep="\t", index=False)
        print(f"After delta-PSI-direction filtering: {len(df_stage2)} rows ({len(df_stage2['junction'].unique())} unique junctions).")

    # Stage 3: Filter by event_types
    if args.event_types:
        df_stage3 = df_stage2.copy() if "df_stage2" in locals() else outlier_df.copy()
        df_stage3 = df_stage3[df_stage3['event'].apply(lambda x: any(evt in x.split(';') for evt in args.event_types))]
        outfile_stage3 = (outfile_stage2 if "outfile_stage2" in locals() else
                          outfile_stage1).replace(".tsv", f"_event.tsv")
        df_stage3.to_csv(outfile_stage3, sep="\t", index=False)
        print(f"After event-type filtering: {len(df_stage3)} rows ({len(df_stage3['junction'].unique())} unique junctions).")

    # Stage 4: Filter by annotation_types
    if args.include_upreg_annotation_types or args.include_downreg_annotation_types:
        df_stage4 = df_stage3.copy() if "df_stage3" in locals() else (
            df_stage2.copy() if "df_stage2" in locals() else outlier_df.copy()
        )
        mask = pd.Series(False, index=df_stage4.index)
        if args.include_upreg_annotation_types and args.include_downreg_annotation_types:
            mask = (
                ((df_stage4["delta_PSI"] > 0) & df_stage4["annotation"].isin(args.include_upreg_annotation_types)) |
                ((df_stage4["delta_PSI"] < 0) & df_stage4["annotation"].isin(args.include_downreg_annotation_types))
            )
        elif args.include_upreg_annotation_types:
            mask = (
                ((df_stage4["delta_PSI"] > 0) & df_stage4["annotation"].isin(args.include_upreg_annotation_types)) |
                (df_stage4["delta_PSI"] < 0)
            )
        elif args.include_downreg_annotation_types:
            mask = (
                ((df_stage4["delta_PSI"] < 0) & df_stage4["annotation"].isin(args.include_downreg_annotation_types)) |
                (df_stage4["delta_PSI"] > 0)
            )
        df_stage4 = df_stage4[mask]
        outfile_stage4 = (outfile_stage3 if "outfile_stage3" in locals() else
                          outfile_stage2 if "outfile_stage2" in locals() else
                          outfile_stage1).replace(".tsv", "_annotation.tsv")
        df_stage4.to_csv(outfile_stage4, sep="\t", index=False)
        print(f"After annotation-type filtering: {len(df_stage4)} rows ({len(df_stage4['junction'].unique())} unique junctions).")

    # Stage 5: Filter by delta_PSI_vs_cohort threshold
    if args.filter_by_cohort_IQR or ((args.fit_beta_dist_on_cohort or args.compare_PSI_to_cohort_median) and args.delta_PSI_vs_cohort_threshold):
        df_stage5 = df_stage4.copy() if "df_stage4" in locals() else (
            df_stage3.copy() if "df_stage3" in locals() else (
                df_stage2.copy() if "df_stage2" in locals() else outlier_df.copy()
            )
        )
        
        # Apply threshold filtering
        if args.filter_by_cohort_IQR and 'cohort_q1' in df_stage5.columns and 'cohort_q3' in df_stage5.columns:
            iqr = df_stage5['cohort_q3'] - df_stage5['cohort_q1']
            lower_fence = df_stage5['cohort_q1'] - 1.5 * iqr
            upper_fence = df_stage5['cohort_q3'] + 1.5 * iqr
            sample_psi = pd.to_numeric(df_stage5["sample_PSI"], errors="coerce")
            is_outlier = (sample_psi < lower_fence) | (sample_psi > upper_fence)
            df_stage5 = df_stage5[is_outlier]
        elif args.delta_PSI_vs_cohort_threshold and 'delta_PSI_vs_cohort' in df_stage5.columns:
            delta_psi_cohort = pd.to_numeric(df_stage5['delta_PSI_vs_cohort'], errors='coerce')
            df_stage5 = df_stage5[delta_psi_cohort.abs() >= args.delta_PSI_vs_cohort_threshold]

        outfile_stage5 = (outfile_stage4 if "outfile_stage4" in locals() else
                          outfile_stage3 if "outfile_stage3" in locals() else
                          outfile_stage2 if "outfile_stage2" in locals() else
                          outfile_stage1).replace(".tsv", "_cohortIQR.tsv" if args.filter_by_cohort_IQR else f"_{args.delta_PSI_vs_cohort_threshold}deltaPSIvsCohort.tsv")
        df_stage5.to_csv(outfile_stage5, sep="\t", index=False)
        print(f"    After cohort-level filtering: {len(df_stage5)} rows ({len(df_stage5['junction'].unique())} unique junctions).")

    # Stage 6: Filter by sample_number_threshold
    if args.sample_number_threshold:
        df_stage6 = df_stage5.copy() if "df_stage5" in locals() else (
            df_stage4.copy() if "df_stage4" in locals() else (
                df_stage3.copy() if "df_stage3" in locals() else (
                    df_stage2.copy() if "df_stage2" in locals() else outlier_df.copy()
                )
            )
        )
        df_stage6 = df_stage6[df_stage6['sample_count'] <= args.sample_number_threshold]
        outfile_stage6 = (outfile_stage5 if "outfile_stage5" in locals() else
                        outfile_stage4 if "outfile_stage4" in locals() else
                        outfile_stage3 if "outfile_stage3" in locals() else
                        outfile_stage2 if "outfile_stage2" in locals() else
                        outfile_stage1).replace(".tsv", f"_{args.sample_number_threshold}samples.tsv")
        df_stage6.to_csv(outfile_stage6, sep="\t", index=False)
        print(f"After sample-number filtering: {len(df_stage6)} rows ({len(df_stage6['junction'].unique())} unique junctions).")

    # Plotting if requested
    if args.plot:
        dfs_for_plot = []
        legends = []
        if 'outlier_df' in locals() and not outlier_df.empty:
            dfs_for_plot.append(outlier_df)
            legends.append(('#d95d5b', f'coverage ≥ {args.jxn_coverage_threshold}, abs(ΔPSI) ≥ {args.delta_PSI_threshold}, and padj ≤ {args.padj_threshold}'))
        if 'df_stage2' in locals() and not df_stage2.empty:
            dfs_for_plot.append(df_stage2)
            legends.append(('#4c8fca', f"+ upregulated junctions (delta PSI > 0)" if args.delta_PSI_direction == "positive" else f"+ downregulated junctions (delta PSI < 0)"))
        if 'df_stage3' in locals() and not df_stage3.empty:
            dfs_for_plot.append(df_stage3)
            legends.append(('#57aa3e', f"+ event types: \n        " + ',\n        '.join(args.event_types)))
        if 'df_stage4' in locals() and not df_stage4.empty:
            dfs_for_plot.append(df_stage4)
            legend_parts = []
            if args.include_upreg_annotation_types:
                legend_parts.append("+ upregulated junctions with annotation:\n        " + ",\n        ".join(args.include_upreg_annotation_types))
            elif args.exclude_upreg_annotation_types:
                legend_parts.append("- upregulated junctions with annotation:\n        " + ",\n        ".join(args.exclude_upreg_annotation_types))
            else:
                legend_parts.append("+ all upregulated junctions")
            if args.include_downreg_annotation_types:
                legend_parts.append("+ downregulated junctions with annotation:\n        " + ",\n        ".join(args.include_downreg_annotation_types))
            elif args.exclude_downreg_annotation_types:
                legend_parts.append("- downregulated junctions with annotation:\n        " + ",\n        ".join(args.exclude_downreg_annotation_types))
            else:
                legend_parts.append("+ all downregulated junctions")
            legends.append(('#f48f3e', "\n".join(legend_parts)))
        if 'df_stage5' in locals() and not df_stage5.empty:
            dfs_for_plot.append(df_stage5)
            legends.append(('#b271ab', '+ IQR outlier' if args.filter_by_cohort_IQR else f'+ abs(ΔPSI_vs_cohort) ≥ {args.delta_PSI_vs_cohort_threshold}'))
        if 'df_stage6' in locals() and not df_stage6.empty:
            dfs_for_plot.append(df_stage6)
            legends.append(('#42b4b5', f'+ outlier in ≤ {args.sample_number_threshold} samples'))
        plot_outlier_counts(samples, dfs_for_plot, legends, args.title, f"{args.outprefix}_counts.pdf")
        plot_outlier_types(outlier_df,
                        f"Outlier Junction (coverage≥{args.jxn_coverage_threshold}, padj≤{args.padj_threshold}, abs(ΔPSI)≥{args.delta_PSI_threshold}) Types",
                        f"{args.outprefix}_types.pdf")

if __name__ == "__main__":
    main()
