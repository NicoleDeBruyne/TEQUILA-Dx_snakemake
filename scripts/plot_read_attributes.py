#!/usr/bin/env python3
"""
scripts/plot_read_attributes.py
Plots per-read length and quality distributions across the cohort, split by
on-target / off-target / unmapped status. Invoked by rules/9_plot_cohort_info.smk
(_9B).

Adapted from a FASTQ-based script pulled from Github 2026.04.08. This version
reads directly from BAM files instead of requiring a paired FASTQ (samples in
this pipeline's run config are provided as BAM only -- see README.md).

Read length and mean Phred quality are taken from each *primary* alignment's
SEQ/QUAL fields (secondary/supplementary records are skipped, since these can
be hard-clipped and would understate the true read length). For aligners that
don't hard-clip primary alignments (true of minimap2/nanopore long-read BAMs,
which this pipeline targets), this reproduces the FASTQ-derived length/quality
closely. If a BAM was written without quality strings (QUAL == "*"), Phred is
reported as 0 for that read; check the input BAM if quality plots look off.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import pysam
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42


def parse_args():
    """ Parse command-line arguments for read length and quality plotting """

    parser = argparse.ArgumentParser(description="Make a plot of read lengths and qualities from BAM files")
    parser.add_argument(
        "--mapping-file",
        required=True,
        help=("TSV file without a header and with the columns:"
              " name, bam, bed. The bed column is optional and"
              " defines the target regions for that row"))
    parser.add_argument(
        "--outprefix",
        required=True,
        help=("Prefix for output files _read_attributes.tsv,"
              " _summary.tsv, _read_lengths_boxplot.pdf, _read_qualities_boxplot.pdf"))
    parser.add_argument('--title')
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads")
    return parser.parse_args()


def get_on_target_readIDs(bam, bed):
    """Get on-target read IDs from a BAM file based on regions defined in a BED file."""

    on_target_readIDs = set()
    with pysam.AlignmentFile(bam, "rb", threads=2) as bamfile:
        with open(bed) as bedfile:
            for line in bedfile:
                if line.strip() and not line.startswith('#'):
                    chrom, start, end = line.split()[:3]
                    for read in bamfile.fetch(chrom, int(start), int(end)):
                        if not (read.is_unmapped or read.is_secondary or read.is_supplementary):
                            on_target_readIDs.add(read.query_name)
    return on_target_readIDs


def process_sample(name, bam, bed):
    """Process one sample directly from its BAM; returns list of dicts"""

    on_target_readIDs = get_on_target_readIDs(bam, bed) if bed else set()

    records = []
    with pysam.AlignmentFile(bam, "rb", threads=2) as bamfile:
        for read in bamfile.fetch(until_eof=True):
            # Skip secondary/supplementary records: they can be hard-clipped,
            # which would understate the read's true length, and would
            # double-count the same underlying read alongside its primary.
            if read.is_secondary or read.is_supplementary:
                continue

            seq_len = read.query_length or 0
            quals = read.query_qualities
            # Plain sum()/len() instead of np.mean(): per-call numpy overhead
            # dominates for arrays this small, and this runs once per read --
            # potentially millions of times per BAM.
            mean_phred = (sum(quals) / len(quals)) if quals is not None and len(quals) > 0 else 0.0

            if read.is_unmapped:
                target_type = "unmapped"
            elif bed:
                target_type = "on_target" if read.query_name in on_target_readIDs else "off_target"
            else:
                target_type = "mapped"

            records.append({
                "Sample": name,
                "ReadLength": seq_len,
                "Phred": mean_phred,
                "TargetType": target_type
            })

    return records


def main():
    """ Main function """

    args = parse_args()

    df_mapping = pd.read_csv(args.mapping_file, sep="\t", header=None)
    if df_mapping.shape[1] == 2:
        df_mapping.columns = ["name", "bam"]
        df_mapping["bed"] = None
    else:
        df_mapping.columns = ["name", "bam", "bed"]
    sample_order = df_mapping["name"].tolist()
    samples = df_mapping.to_dict(orient="records")

    print(f"Processing {len(df_mapping)} samples using {args.threads} threads...")

    dfs = []

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_sample, s["name"], s["bam"], s["bed"]): s["name"] for s in samples}
        for future in as_completed(futures):
            name = futures[future]
            try:
                result = future.result()
                dfs.append(pd.DataFrame(result))
                print(f"Finished {name}")
            except Exception as e:
                print(f"\u274c Error in {name}: {e}")

    print("Merging attribute info...")
    df_all = pd.concat(dfs, ignore_index=True)

    df_all = df_all.sort_values('Sample', key=lambda x: x.map({s: i for i, s in enumerate(sample_order)}))
    df_all.to_csv(f"{args.outprefix}_read_attributes.tsv", sep="\t", index=False)
    print(f"\n\u2705 Attributes saved to {args.outprefix}_read_attributes.tsv")

    print("\nGenerating summary statistics...")

    summary_overall = (
        df_all
        .groupby("Sample")
        .agg(
            MeanReadLength=("ReadLength", "mean"),
            MedianReadLength=("ReadLength", "median"),
            MeanPhred=("Phred", "mean"),
            MedianPhred=("Phred", "median"),
            TotalReads=("ReadLength", "count")
        )
        .reset_index()
    )
    summary_overall["TargetType"] = "ALL"

    summary_by_type = (
        df_all
        .groupby(["Sample", "TargetType"])
        .agg(
            MeanReadLength=("ReadLength", "mean"),
            MedianReadLength=("ReadLength", "median"),
            MeanPhred=("Phred", "mean"),
            MedianPhred=("Phred", "median"),
            TotalReads=("ReadLength", "count")
        )
        .reset_index()
    )

    summary = pd.concat([summary_overall, summary_by_type], ignore_index=True)

    default_target_order = ["on_target", "off_target", "mapped", "unmapped"]
    target_order = default_target_order + [
        t for t in summary['TargetType'].drop_duplicates()
        if t not in default_target_order
    ]
    sample_order = [s for s in sample_order if s in summary['Sample'].unique()]
    target_order = [t for t in target_order if t != "ALL" and t in summary['TargetType'].unique()]
    summary['Sample_order'] = summary['Sample'].map({s: i for i, s in enumerate(sample_order)})
    summary['Target_order'] = summary['TargetType'].map({t: i for i, t in enumerate(target_order)})
    summary_sorted = summary.sort_values(['Sample_order', 'Target_order']).drop(columns=['Sample_order', 'Target_order'])

    summary_file = f"{args.outprefix}_summary.tsv"
    summary_sorted.to_csv(summary_file, sep="\t", index=False)

    print(f"\u2705 Summary statistics saved to {summary_file}")

    # NOTE: deliberately not groupby(...).apply(...) here -- on newer pandas
    # (>=2.2), apply() drops the grouping columns from a function's result
    # whenever that result already contains them, turning "Sample"/"TargetType"
    # into index levels only. Every later re-groupby on df_ds by those same
    # column names would then raise KeyError. A manual per-group loop + concat
    # keeps them as ordinary columns regardless of pandas version.
    df_ds = pd.concat(
        [g.sample(min(len(g), 100_000), random_state=0) for _, g in df_all.groupby(["Sample", "TargetType"])],
        ignore_index=True,
    )

    palette_dict = dict(zip(target_order, sns.color_palette("pastel")))

    print("\nPlotting read length boxplot...")
    fig, ax = plt.subplots(figsize=(max(6, df_all['Sample'].nunique() * 0.5), 6))
    sns.boxplot(data=df_ds, x="Sample", y="ReadLength", hue="TargetType", order=sample_order,
                hue_order=target_order, palette=palette_dict, showfliers=False, width=0.7, ax=ax
                )
    ymax = ax.get_ylim()[1]
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, ymax * 1.05)
    ax.margins(y=0.05)
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(f"{args.outprefix}_read_lengths_boxplot.pdf")
    print(f"\u2705 Read length box plot saved to {args.outprefix}_read_lengths_boxplot.pdf")

    print("\nPlotting read quality boxplot...")
    fig, ax = plt.subplots(figsize=(max(6, df_all['Sample'].nunique() * 0.5), 6))
    sns.boxplot(data=df_ds, x="Sample", y="Phred", hue="TargetType", order=sample_order,
                hue_order=target_order, palette=palette_dict, showfliers=False, width=0.7, ax=ax
                )
    ymax = ax.get_ylim()[1]
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, ymax * 1.05)
    ax.margins(y=0.05)
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(f"{args.outprefix}_read_qualities_boxplot.pdf")
    print(f"\u2705 Read quality box plot saved to {args.outprefix}_read_qualities_boxplot.pdf")

    print("\nPlotting read length violin plot...")

    def filter_iqr(group):
        x = group["ReadLength"]
        q1 = x.quantile(0.25)
        q3 = x.quantile(0.75)
        iqr = q3 - q1
        upper_iqr = q3 + 3 * iqr
        upper_p95 = x.quantile(0.95)
        upper = max(upper_iqr, upper_p95)
        return group[x <= upper]

    # Same manual-loop reasoning as df_ds above.
    df_ds_filtered = pd.concat(
        [filter_iqr(g) for _, g in df_ds.groupby(["Sample", "TargetType"])],
        ignore_index=True,
    )
    fig, ax = plt.subplots(figsize=(max(6, df_all['Sample'].nunique() * 0.5), 6))
    sns.violinplot(data=df_ds_filtered, x="Sample", y="ReadLength", hue="TargetType", order=sample_order,
                    hue_order=target_order, palette=palette_dict, width=0.8, inner="quartile",
                    dodge=True, density_norm="width", ax=ax
                    )
    handles, labels = ax.get_legend_handles_labels()
    n = len(target_order)
    ax.legend(handles[:n], labels[:n])
    ymax = ax.get_ylim()[1]
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, ymax * 1.05)
    ax.margins(y=0.05)
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(f"{args.outprefix}_read_lengths_violin.pdf")
    print(f"\u2705 Read length violin plot saved to {args.outprefix}_read_lengths_violin.pdf")

    print("\nPlotting read quality violin plot...")
    fig, ax = plt.subplots(figsize=(max(6, df_all['Sample'].nunique() * 0.5), 6))
    sns.violinplot(data=df_ds, x="Sample", y="Phred", hue="TargetType", order=sample_order,
                    hue_order=target_order, palette=palette_dict, width=0.8, inner="quartile",
                    dodge=True, density_norm="width", ax=ax
                    )
    handles, labels = ax.get_legend_handles_labels()
    n = len(target_order)
    ax.legend(handles[:n], labels[:n])
    ymax = ax.get_ylim()[1]
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, ymax * 1.05)
    ax.margins(y=0.05)
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(f"{args.outprefix}_read_qualities_violin.pdf")
    print(f"\u2705 Read quality violin plot saved to {args.outprefix}_read_qualities_violin.pdf")


if __name__ == "__main__":
    main()