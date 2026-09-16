#!/usr/bin/env python3
"""
scripts/plot_on_target_rates.py
Cohort-level merge step: combines every sample's own {sample}_on_target.tsv
(written per-sample by scripts/get_on_target_rate.py via rules/6_sample_qc.smk's
_6A) into one cohort table, derives mapping_rate/on_target_rate from the raw
total/mapped/target counts, and plots them across the cohort. No BAM access
here -- adding/removing a sample from a cohort only reruns this cheap merge,
not the per-sample BAM scan. Invoked by rules/9_merge_results.smk's _9G.

Writes {outprefix}_on_target_rates.tsv, {outprefix}_mapping_rates.pdf, and
{outprefix}_on_target_rates.pdf (+ an alias-labeled copy of the latter).
"""

import argparse

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

from sample_alias import add_alias_map_arg, parse_alias_map, resolve


def parse_args():
    parser = argparse.ArgumentParser(description="Merge per-sample on-target TSVs and plot mapping/on-target rates across a cohort")
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Every sample's {sample}_on_target.tsv, in the desired plot order")
    parser.add_argument("--groups", nargs="*", default=[],
        help="Optional per-sample group/sample_type label (same order as --infiles) for bar coloring")
    parser.add_argument("--outprefix", required=True,
        help="Prefix for output files: <outprefix>_on_target_rates.tsv, _mapping_rates.pdf, _on_target_rates.pdf")
    parser.add_argument('--title')
    add_alias_map_arg(parser)
    return parser.parse_args()


def fmt(n):
    return f'{n/1e6:.1f}M' if n >= 1e6 else f'{n//1e3}k' if n >= 1e3 else str(n)


def rate_label(r):
    return '100%' if r == 100 else '>99%' if r >= 99 else f"{r:.0f}%" if r >= 10 else f"{r:.1f}%" if r >= 0.1 else '<0.1%'


def plot_mapping_figure(df, outprefix, width, title):
    fig, ax = plt.subplots(2, figsize=(width, 10))

    ax[0].bar(range(len(df)), df['mapping_rate'], color="#FFD676")  # yellow
    ax[0].bar(range(len(df)), 100 - df['mapping_rate'], bottom=df['mapping_rate'], color='#C4C4C4')  # grey
    for i, r in df.iterrows():
        ax[0].text(i, r['mapping_rate'] + 1, rate_label(r['mapping_rate']), ha='center', fontsize=8)
    ax[0].set_ylim(0, 105)
    ax[0].set_yticks([0, 20, 40, 60, 80, 100])
    ax[0].set_ylabel("Mapping Rate (%)")
    ax[0].set_xticks(range(len(df)))
    ax[0].set_xticklabels([])

    ax[1].bar(range(len(df)), df['mapped'] / 1e6, color="#FFD676")  # yellow
    ax[1].bar(range(len(df)), df['unmapped'] / 1e6, bottom=df['mapped'] / 1e6, color='#C4C4C4')  # grey
    offset = df['total'].max() * 0.02
    for i, r in df.iterrows():
        ax[1].text(i, (r['total'] + offset) / 1e6, f"{fmt(r['total'])}", ha='center', fontsize=8)
    ax[1].set_ylim(0, df['total'].max() * 1.05 / 1e6)
    ax[1].set_ylabel("Read Count (M)")
    ax[1].set_xticks(range(len(df)))
    ax[1].set_xticklabels(df['sample'], rotation=45, ha='right')

    if title:
        fig.suptitle(title)
    plt.tight_layout()
    out = f"{outprefix}_mapping_rates.pdf"
    plt.savefig(out, bbox_inches='tight')
    print(f"Saved figure: {out}")


def plot_ontarget_figure(df, bar_colors, outprefix, suffix, width, title):
    fig, ax = plt.subplots(2, figsize=(width, 10))

    ax[0].bar(range(len(df)), df['on_target_rate'], color=bar_colors)
    ax[0].bar(range(len(df)), 100 - df['on_target_rate'], bottom=df['on_target_rate'], color='#C4C4C4')  # grey
    for i, r in df.iterrows():
        ax[0].text(i, r['on_target_rate'] + 1, rate_label(r['on_target_rate']), ha='center', fontsize=8)
    ax[0].set_ylim(0, 105)
    ax[0].set_yticks([0, 20, 40, 60, 80, 100])
    ax[0].set_ylabel("On-target Rate (%)")
    ax[0].set_xticks(range(len(df)))
    ax[0].set_xticklabels([])

    off = df['mapped'] - df['target']
    ax[1].bar(range(len(df)), df['target'] / 1e6, color=bar_colors)
    ax[1].bar(range(len(df)), off / 1e6, bottom=df['target'] / 1e6, color='#C4C4C4')  # grey
    offset = df['mapped'].max() * 0.02
    for i, r in df.iterrows():
        ax[1].text(i, (r['mapped'] + offset) / 1e6, f"{fmt(r['mapped'])}", ha='center', fontsize=8)
    ax[1].set_ylim(0, df['mapped'].max() * 1.05 / 1e6)
    ax[1].set_ylabel("Read Count (M)")
    ax[1].set_xticks(range(len(df)))
    ax[1].set_xticklabels(df['sample'], rotation=45, ha='right')

    if title:
        fig.suptitle(title)
    plt.tight_layout()
    out = f"{outprefix}_on_target_rates{suffix}.pdf"
    plt.savefig(out, bbox_inches='tight')
    print(f"Saved figure: {out}")


def main():
    args = parse_args()

    per_sample_dfs = [pd.read_csv(f, sep="\t") for f in args.infiles]
    df = pd.concat(per_sample_dfs, ignore_index=True)

    groups = list(args.groups) if args.groups else [None] * len(args.infiles)
    if len(groups) != len(args.infiles):
        raise ValueError("--groups must have one entry per --infiles entry, or be omitted entirely")
    group_by_sample = dict(zip([d['sample'].iloc[0] for d in per_sample_dfs], groups))
    df['group'] = df['sample'].map(group_by_sample)

    df['unmapped'] = df['total'] - df['mapped']
    df['mapping_rate'] = df['mapped'] / df['total'] * 100
    df['on_target_rate'] = df['target'] / df['mapped'] * 100

    out_tsv = f"{args.outprefix}_on_target_rates.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved TSV: {out_tsv}")

    unique_groups = list(df['group'].dropna().unique())
    colors = [
        '#6997B9',  # blue
        '#BB6A68',  # red
        '#70A677',  # green
        '#D48653',  # orange
        '#A783A3',  # purple
    ]
    color_dict = {g: colors[i % len(colors)] for i, g in enumerate(unique_groups)} if unique_groups else {}
    bar_colors = df['group'].map(color_dict).fillna('#6997B9')  # blue

    # Linear scaling: reserve a fixed amount of horizontal space per sample so
    # labels don't overlap, with a floor for small cohorts. No upper cap --
    # large cohorts (e.g. AGS390 has ~200 samples) need a genuinely wide PDF.
    width = max(10, 0.22 * len(df))

    plot_mapping_figure(df, args.outprefix, width, args.title)
    plot_ontarget_figure(df, bar_colors, args.outprefix, "", width, args.title)

    # Alias-labeled copy: always produced (mirrors the real-ID figure
    # verbatim when --alias-map is empty), so the rule's declared output
    # exists regardless of whether this bed panel actually has any aliases.
    alias_map = parse_alias_map(args.alias_map)
    alias_df = df.copy()
    alias_df['sample'] = alias_df['sample'].apply(lambda s: resolve(s, alias_map))
    plot_ontarget_figure(alias_df, bar_colors, args.outprefix, "_alias", width, args.title)


if __name__ == "__main__":
    main()
