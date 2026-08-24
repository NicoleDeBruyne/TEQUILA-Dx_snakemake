#!/usr/bin/env python3

# scripts/plot_hits_upset.py
# UpSet-style breakdown of merged_all_hits.tsv's four hit categories --
# variant, ASE, outlier_junction (GTEx comparison), cohort_outlier_junction
# (cohort comparison) -- across a BED panel's whole cohort (all sample_types
# pooled, same scope as merged_all_hits.tsv itself).
#
# Unlike a standard UpSet plot (one bar per category combination, height =
# total element count), this shows the *distribution across samples* for
# each combination: a boxplot of "how many genes did this sample have in
# this exact combination of categories", with one dot per sample (colored
# by sample_type) overlaid on each box. The combination-membership matrix
# is drawn as a second panel below the boxplot, in the usual UpSet style.

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

rcParams['pdf.fonttype'] = 42

CATEGORIES = ['variant', 'ASE', 'junction_outlier', 'cohort_junction_outlier']

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
        description="Plot an UpSet-style boxplot of per-sample gene counts across every non-empty "
                    "combination of variant/ASE/outlier_junction/cohort_outlier_junction hit categories.")
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
    df['_junction_hit'] = ~df['outlier_junction'].isin(['None', '.'])
    df['_cohort_junction_hit'] = ~df['cohort_outlier_junction'].isin(['None', '.'])

    hit_cols = {
        'variant': '_variant_hit',
        'ASE': '_ase_hit',
        'junction_outlier': '_junction_hit',
        'cohort_junction_outlier': '_cohort_junction_hit',
    }

    def row_combo(row):
        return tuple(cat for cat in CATEGORIES if row[hit_cols[cat]])

    df['_combo'] = df.apply(row_combo, axis=1)
    # Drop genes with no hit in any of the four categories -- standard
    # UpSet convention of only showing non-empty intersections. A gene
    # only appears in merged_all_hits.tsv at all if it matched via at
    # least one of variant/ASE/junction/cohort_junction, so this should
    # only ever drop a small number of rows, if any.
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

    # ---------- Figure: boxplot+dots on top, combination matrix below ----------
    n_combos = len(combo_order)
    width = min(30, 3 + 0.6 * n_combos)
    fig, (ax_box, ax_matrix) = plt.subplots(
        2, 1, figsize=(width, 8), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]},
    )

    positions = np.arange(1, n_combos + 1)
    box_data = [counts[combo].values for combo in combo_order]
    ax_box.boxplot(
        box_data, positions=positions, widths=0.5,
        patch_artist=True, boxprops=dict(facecolor='none', color='black'),
        medianprops=dict(color='black'), whiskerprops=dict(color='black'),
        capprops=dict(color='black'), flierprops=dict(marker=''),
    )

    rng = np.random.default_rng(0)
    for pos, combo in zip(positions, combo_order):
        vals = counts[combo].values
        jitter = rng.normal(pos, 0.06, size=len(vals))
        dot_colors = [color_map[sample_type_by_sample[s]] for s in counts.index]
        ax_box.scatter(jitter, vals, c=dot_colors, s=30, alpha=0.8, zorder=3, edgecolors='none')

    ax_box.set_ylabel("Gene count")
    ax_box.set_title(args.title)
    ax_box.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Legend: one entry per sample_type, in the same color-assignment order.
    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color=color_map[st], label=st, markersize=7)
        for st in dict.fromkeys(args.sample_types)
    ]
    ax_box.legend(handles=handles, title="sample_type", bbox_to_anchor=(1.02, 1), loc="upper left")

    # Combination matrix (standard UpSet style): one row per category, one
    # column per combo, filled dot if that category is in the combo, faint
    # open dot otherwise, with a connecting line through the filled dots.
    n_cats = len(CATEGORIES)
    cat_y = {cat: n_cats - i for i, cat in enumerate(CATEGORIES)}
    for pos, combo in zip(positions, combo_order):
        in_combo_y = [cat_y[cat] for cat in CATEGORIES if cat in combo]
        if len(in_combo_y) > 1:
            ax_matrix.plot([pos, pos], [min(in_combo_y), max(in_combo_y)], color='black', lw=1.5, zorder=1)
        for cat in CATEGORIES:
            filled = cat in combo
            ax_matrix.scatter(
                [pos], [cat_y[cat]],
                s=80, zorder=2,
                color='black' if filled else '#DDDDDD',
            )
    ax_matrix.set_yticks(list(cat_y.values()))
    ax_matrix.set_yticklabels(list(cat_y.keys()))
    ax_matrix.set_ylim(0.5, n_cats + 0.5)
    ax_matrix.set_xticks(positions)
    ax_matrix.set_xticklabels([])
    ax_matrix.set_xlim(0.5, n_combos + 0.5)
    for spine in ('top', 'right', 'bottom'):
        ax_matrix.spines[spine].set_visible(False)

    plt.tight_layout()
    outfile = f"{args.outdir}/hits_upset_plot.pdf"
    plt.savefig(outfile)
    plt.close(fig)
    print(f"Saved plot to {outfile}")


if __name__ == "__main__":
    main()
