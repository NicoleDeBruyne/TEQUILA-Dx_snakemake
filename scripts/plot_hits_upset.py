#!/usr/bin/env python3

# scripts/plot_hits_upset.py
# UpSet-style breakdown of merged_all_hits.tsv's three hit categories --
# variant, ASE, junction (a hit in EITHER the GTEx comparison or the
# cohort comparison counts as one 'junction' hit) -- across a BED panel's
# whole cohort (all sample_types pooled, same scope as merged_all_hits.tsv
# itself).
#
# Unlike a standard UpSet plot (one bar per category combination, height =
# total element count), this shows the *distribution across samples* for
# each combination: a swarm plot of "how many genes did this sample have in
# this exact combination of categories", one non-overlapping dot per sample
# (colored by sample_type). Each combination gets its own group of
# side-by-side sub-swarms, one per sample_type, rather than mixing every
# sample_type into a single swarm -- alternating background shading marks
# each combination's group. A swarm plot (rather than a boxplot) is used
# because gene counts are small non-negative integers -- quartiles computed
# over a handful of discrete values are not a meaningful summary, and a
# swarm shows the actual per-sample counts (including exact ties) instead
# of implying a continuous distribution that isn't there. The
# combination-membership matrix is drawn as a second panel below the swarm,
# centered under each combination's group, in the usual UpSet style.

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams

rcParams['pdf.fonttype'] = 42

CATEGORIES = ['variant', 'ASE', 'junction']

# Same palette + assignment scheme as scripts/plot_on_target_rates.py's
# per-sample_type bar coloring: unique sample_types in order of first
# appearance, cycled through this 5-color list. Reused here (rather than a
# new palette) so sample_type colors stay visually consistent across the
# pipeline's plots.
_SAMPLE_TYPE_COLORS = [
    '#6997B9',  # blue
    '#BB6A68',  # red
    '#70A677',  # green
    '#D48653',  # orange
    '#A783A3',  # purple
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot an UpSet-style swarm plot of per-sample gene counts across every non-empty "
                    "combination of variant/ASE/junction hit categories.")
    parser.add_argument("--infile", required=True, help="Path to merged_all_hits.tsv.")
    parser.add_argument("--samples", nargs="+", required=True,
        help="Every sample in this BED panel's cohort (same order as --sample-types), so samples with "
             "zero genes in a given combination still show up as a 0 rather than being omitted.")
    parser.add_argument("--sample-types", nargs="+", required=True,
        help="Each --samples entry's sample_type, same order/length as --samples.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--title", default="Candidate Gene Hit Categories by Sample", help="Plot title.")
    args = parser.parse_args()
    if len(args.samples) != len(args.sample_types):
        parser.error("--samples and --sample-types must have the same number of entries")
    return args


def _sample_type_color_map(sample_types_in_order):
    """Unique sample_types in order of first appearance -> hex color,
    cycling through _SAMPLE_TYPE_COLORS -- same logic as
    plot_on_target_rates.py's color_dict."""
    seen = list(dict.fromkeys(sample_types_in_order))  # unique, first-appearance order
    return {st: _SAMPLE_TYPE_COLORS[i % len(_SAMPLE_TYPE_COLORS)] for i, st in enumerate(seen)}


def _combo_label(combo):
    """('variant', 'ASE') -> 'variant & ASE'"""
    return " & ".join(combo) if combo else "(none)"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sample_type_by_sample = dict(zip(args.samples, args.sample_types))
    color_map = _sample_type_color_map(args.sample_types)

    df = pd.read_csv(args.infile, sep='\t', keep_default_na=False)
    df = df.astype(object)
    df.fillna('.', inplace=True)

    # 'variant'/'ASE' are already real bool (pandas auto-infers a
    # pure-"True"/"False" column as bool on read, same as
    # plot_candidate_hits.py relies on). 'outlier_junction'/
    # 'cohort_outlier_junction' are Strong/Moderate/Weak/None/'.' strings --
    # "hit" means anything other than no-evidence ('None') or not-annotated
    # ('.', e.g. when cohort junction data wasn't available for this run).
    df['_variant_hit'] = df['variant'].astype(bool)
    df['_ase_hit'] = df['ASE'].astype(bool)
    df['_junction_hit'] = (
        ~df['outlier_junction'].isin(['None', '.'])
        | ~df['cohort_outlier_junction'].isin(['None', '.'])
    )

    hit_cols = {
        'variant': '_variant_hit',
        'ASE': '_ase_hit',
        'junction': '_junction_hit',
    }

    def row_combo(row):
        return tuple(cat for cat in CATEGORIES if row[hit_cols[cat]])

    df['_combo'] = df.apply(row_combo, axis=1)
    # Drop genes with no hit in any of the three categories -- standard
    # UpSet convention of only showing non-empty intersections. A gene
    # only appears in merged_all_hits.tsv at all if it matched via at
    # least one of variant/ASE/junction, so this should only ever drop a
    # small number of rows, if any.
    df = df[df['_combo'].apply(len) > 0]

    # Per-sample, per-combo gene counts, reindexed over every sample in
    # --samples (not just ones appearing in df) and every combo that
    # occurs anywhere in the cohort, so a sample with zero genes in a
    # given combo shows up as an explicit 0 dot rather than being silently
    # skipped in that box's distribution.
    counts = df.groupby(['sample', '_combo'])['gene'].nunique().unstack(fill_value=0)
    counts = counts.reindex(index=args.samples, fill_value=0)

    if counts.shape[1] == 0:
        print("WARNING: no non-empty category combinations found in the input -- nothing to plot.")
        # Still write an (empty) counts file and an (empty) placeholder PDF so
        # this rule's declared outputs exist.
        pd.DataFrame(columns=['sample', 'sample_type', 'combo', 'gene_count']).to_csv(
            f"{args.outdir}/hits_upset_counts.tsv", sep='\t', index=False)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No candidate hits to plot", ha='center', va='center')
        ax.axis('off')
        plt.savefig(f"{args.outdir}/hits_upset_plot.pdf")
        plt.close(fig)
        return

    # Classic UpSet ordering: largest intersections (by total gene count
    # summed across every sample) first.
    combo_order = counts.sum(axis=0).sort_values(ascending=False).index.tolist()
    counts = counts[combo_order]

    # Long-format counts tsv, for reference/debugging alongside the plot.
    long_df = counts.reset_index().melt(id_vars='sample', var_name='combo', value_name='gene_count')
    long_df['combo'] = long_df['combo'].apply(_combo_label)
    long_df['sample_type'] = long_df['sample'].map(sample_type_by_sample)
    long_df = long_df[['sample', 'sample_type', 'combo', 'gene_count']]
    long_df.to_csv(f"{args.outdir}/hits_upset_counts.tsv", sep='\t', index=False)
    print(f"Saved counts tsv to {args.outdir}/hits_upset_counts.tsv")

    # ---------- Figure: swarm plot on top, combination matrix below ----------
    n_combos = len(combo_order)
    sample_types_unique = list(dict.fromkeys(args.sample_types))
    n_st = len(sample_types_unique)
    st_index = {st: i for i, st in enumerate(sample_types_unique)}

    # Long-format frame: one row per (sample, combo).
    swarm_df = counts.reset_index().melt(id_vars='sample', var_name='_combo', value_name='gene_count')
    swarm_df['sample_type'] = swarm_df['sample'].map(sample_type_by_sample)

    # Each combo gets its own sub-swarm per sample_type (rather than one
    # swarm per combo mixing every sample_type together) -- one column of
    # points per (combo, sample_type), grouped side-by-side under that
    # combo. Swarm plots pack same-valued points side-by-side within a
    # column at a fixed marker size -- given enough ties (e.g. many
    # samples with 0 genes), a too-narrow column silently drops points
    # rather than overlapping them. Scale marker size and per-sub-column
    # width to the worst-case tie count (the largest number of samples
    # sharing one exact (combo, sample_type, gene_count) value) so every
    # point has room, instead of assuming a fixed width regardless of
    # cohort size.
    tie_counts = swarm_df.groupby(['_combo', 'sample_type', 'gene_count'])['sample'].transform('count')
    max_tie = int(tie_counts.max())
    swarm_size = 3
    per_subcol_width = max(0.9, 0.07 * max_tie)
    group_gap = 0.5  # inches between adjacent combos' groups
    group_width = n_st * per_subcol_width
    group_spacing = group_width + group_gap

    # Sub-swarm x-offsets within a combo's group, one per sample_type
    # (same order as the legend/color map), centered on the group.
    offsets = (np.arange(n_st) - (n_st - 1) / 2) * per_subcol_width
    offset_by_st = {st: offsets[st_index[st]] for st in sample_types_unique}

    group_center = {combo: i * group_spacing for i, combo in enumerate(combo_order)}
    swarm_df['x'] = swarm_df['_combo'].apply(lambda c: group_center[c]) + swarm_df['sample_type'].map(offset_by_st)

    margin = 0.6
    total_span = (n_combos - 1) * group_spacing + group_width + 2 * margin
    width = min(80, max(6, total_span))

    fig, (ax_box, ax_matrix) = plt.subplots(
        2, 1, figsize=(width, 8), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]},
    )

    # Faint alternating background bands, one per combo group, so
    # neighboring combos' sub-swarms read as visually distinct clusters
    # even though there's no per-sub-column text label.
    for i, combo in enumerate(combo_order):
        if i % 2 == 1:
            band_lo = group_center[combo] - group_width / 2 - group_gap / 4
            band_hi = group_center[combo] + group_width / 2 + group_gap / 4
            ax_box.axvspan(band_lo, band_hi, color='#F0F0F0', zorder=0)
            ax_matrix.axvspan(band_lo, band_hi, color='#F0F0F0', zorder=0)

    with warnings.catch_warnings():
        # We've already sized each sub-column to fit its worst-case tie
        # count above -- suppress seaborn's "points cannot be placed"
        # warning rather than let it surface as pipeline log noise on
        # every run.
        warnings.filterwarnings("ignore", message=".*cannot be placed.*")
        sns.swarmplot(
            data=swarm_df, x='x', y='gene_count', hue='sample_type',
            palette=color_map, ax=ax_box, size=swarm_size, native_scale=True, legend=False,
        )

    ax_box.set_xlabel('')
    ax_box.set_ylabel("Gene count")
    ax_box.set_title(args.title)
    ax_box.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Legend: one entry per sample_type, in the same color-assignment order.
    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color=color_map[st], label=st, markersize=7)
        for st in sample_types_unique
    ]
    ax_box.legend(handles=handles, title="sample_type", bbox_to_anchor=(1.02, 1), loc="upper left")

    # Combination matrix (standard UpSet style): one row per category, one
    # column per combo (at that combo's group center, spanning all its
    # sample_type sub-swarms above), filled dot if that category is in
    # the combo, faint open dot otherwise, with a connecting line through
    # the filled dots.
    n_cats = len(CATEGORIES)
    cat_y = {cat: n_cats - i for i, cat in enumerate(CATEGORIES)}
    for combo in combo_order:
        pos = group_center[combo]
        in_combo_y = [cat_y[cat] for cat in CATEGORIES if cat in combo]
        if len(in_combo_y) > 1:
            ax_matrix.plot([pos, pos], [min(in_combo_y), max(in_combo_y)], color='black', lw=1.5, zorder=2)
        for cat in CATEGORIES:
            filled = cat in combo
            ax_matrix.scatter(
                [pos], [cat_y[cat]],
                s=80, zorder=3,
                color='black' if filled else '#DDDDDD',
            )
    ax_matrix.set_yticks(list(cat_y.values()))
    ax_matrix.set_yticklabels(list(cat_y.keys()))
    ax_matrix.set_ylim(0.5, n_cats + 0.5)
    ax_matrix.set_xticks([group_center[c] for c in combo_order])
    ax_matrix.set_xticklabels([])
    xlim_lo = group_center[combo_order[0]] - group_width / 2 - margin
    xlim_hi = group_center[combo_order[-1]] + group_width / 2 + margin
    ax_matrix.set_xlim(xlim_lo, xlim_hi)
    for spine in ('top', 'right', 'bottom'):
        ax_matrix.spines[spine].set_visible(False)

    plt.tight_layout()
    outfile = f"{args.outdir}/hits_upset_plot.pdf"
    plt.savefig(outfile)
    plt.close(fig)
    print(f"Saved plot to {outfile}")


if __name__ == "__main__":
    main()