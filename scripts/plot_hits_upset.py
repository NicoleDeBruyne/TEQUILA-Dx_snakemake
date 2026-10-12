import argparse
import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

from hits_io import read_hit_rows
from palette import default_type_colors

rcParams['pdf.fonttype'] = 42
rcParams['font.size'] = 13
rcParams['axes.titlesize'] = 15
rcParams['axes.labelsize'] = 13
rcParams['xtick.labelsize'] = 11
rcParams['ytick.labelsize'] = 11
rcParams['legend.fontsize'] = 11
rcParams['legend.title_fontsize'] = 12

CATEGORIES = ['variant', 'ASE', 'junction']

CATEGORY_DISPLAY_ORDER = ['variant', 'junction', 'ASE']

# Layout (x-axis data units): one box per sample type inside each group (combination or tier)
BOX_WIDTH = 0.6          # width of one sample type's box
SUBCOL_GAP = 0.15        # gap between sample types' boxes within a group
GROUP_GAP = 0.6          # gap between groups
MARKER_SIZE = 5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot an UpSet-style raincloud plot of per-sample gene counts across every non-empty "
                    "combination of variant/ASE/junction hit categories.")
    parser.add_argument("--infile", required=True, help="Path to merged_all_hits.tsv.")
    parser.add_argument("--samples", nargs="+", required=True,
        help="Every sample in this BED panel's cohort (same order as --sample-types), so samples with "
             "zero genes in a given combination still show up as a 0 rather than being omitted.")
    parser.add_argument("--sample-types", nargs="+", required=True,
        help="Each --samples entry's sample_type, same order/length as --samples.")
    parser.add_argument("--colors", nargs="*", default=[],
        help="Each --samples entry's colour (its sample type's), same order as --samples. Without it, "
             "sample types get colours in the palette's select order.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--title", default="Candidate Gene Hit Categories by Sample", help="Plot title.")
    args = parser.parse_args()
    if len(args.samples) != len(args.sample_types):
        parser.error("--samples and --sample-types must have the same number of entries")
    return args


def _sample_type_color_map(sample_types_in_order, configured=None):
    return default_type_colors(list(dict.fromkeys(sample_types_in_order)), configured)


def _combo_label(combo):
    return " & ".join(combo) if combo else "(none)"


def _group_layout(keys, sample_types):
    """x positions: group centres, and each sample type's offset from its group centre."""
    n_st = len(sample_types)
    group_width = n_st * BOX_WIDTH + (n_st - 1) * SUBCOL_GAP
    group_spacing = group_width + GROUP_GAP
    offsets = {st: -group_width / 2 + BOX_WIDTH / 2 + i * (BOX_WIDTH + SUBCOL_GAP)
               for i, st in enumerate(sample_types)}
    centers = {k: i * group_spacing for i, k in enumerate(keys)}
    return centers, offsets, group_width


def _draw_box_strip(ax, long_df, keycol, keys, sample_types, color_map, centers, offsets, group_width):
    """For each group x sample type: a box plot of every sample's count (zeros included) with
    the non-zero counts overlaid as a jittered strip. Alternate groups get a light band."""
    rng = np.random.default_rng(0)
    for i, k in enumerate(keys):
        if i % 2 == 1:
            ax.axvspan(centers[k] - group_width / 2 - GROUP_GAP / 4,
                       centers[k] + group_width / 2 + GROUP_GAP / 4, color='#F0F0F0', zorder=0)
    for (k, st), g in long_df.groupby([keycol, 'sample_type'], sort=False):
        if k not in centers or st not in offsets:
            continue
        x = centers[k] + offsets[st]
        vals = g['gene_count'].to_numpy(dtype=float)
        line = dict(color='black', linewidth=0.9)
        ax.boxplot([vals], positions=[x], widths=BOX_WIDTH, showfliers=False, manage_ticks=False, zorder=2,
                   boxprops=line, whiskerprops=line, capprops=line, medianprops=line)
        nz = vals[vals > 0]
        if len(nz):
            ax.scatter(x + rng.uniform(-BOX_WIDTH * 0.35, BOX_WIDTH * 0.35, len(nz)), nz,
                       color=color_map[st], s=MARKER_SIZE ** 2, alpha=0.7, linewidths=0, zorder=3)
    ax.set_ylim(bottom=0)


def _sample_type_legend(ax, sample_types, color_map, st_total_counts):
    handles = [plt.Line2D([0], [0], marker='o', linestyle='', color=color_map[st], markersize=7,
                          label=f"{st} (n={st_total_counts[st]})")
               for st in sample_types]
    ax.legend(handles=handles, title="sample_type", bbox_to_anchor=(1.02, 1), loc="upper left")


def _figure_width(keys, group_width):
    return min(20, max(6, len(keys) * (group_width + GROUP_GAP) + 1.5))


def _plot_tier_density(df, args, sample_type_by_sample, sample_types_unique, st_total_counts, color_map):
    """Genes per tier (1-6) for every sample, split by sample type: box plot of all samples
    (zeros included) with the non-zero counts overlaid as a strip."""
    if 'tier' not in df.columns:
        print("WARNING: no 'tier' column found in the input -- skipping tier gene-count plot.")
        return

    TIERS = [1, 2, 3, 4, 5, 6]
    tier_num = pd.to_numeric(df['tier'], errors='coerce')
    tier_df = df.loc[tier_num.notna()].copy()
    tier_df['tier'] = tier_num.loc[tier_df.index].astype(int)

    counts = tier_df.groupby(['sample', 'tier'])['gene'].nunique().unstack(fill_value=0)
    counts = counts.reindex(index=args.samples, fill_value=0)
    counts = counts.reindex(columns=TIERS, fill_value=0)

    long_df = counts.reset_index().melt(id_vars='sample', var_name='tier', value_name='gene_count')
    long_df['sample_type'] = long_df['sample'].map(sample_type_by_sample)
    long_df = long_df[['sample', 'sample_type', 'tier', 'gene_count']]
    long_df.to_csv(f"{args.outdir}/tier_gene_counts.tsv", sep='\t', index=False)
    print(f"Saved counts tsv to {args.outdir}/tier_gene_counts.tsv")

    centers, offsets, group_width = _group_layout(TIERS, sample_types_unique)
    fig, ax = plt.subplots(figsize=(_figure_width(TIERS, group_width), 6))
    _draw_box_strip(ax, long_df, 'tier', TIERS, sample_types_unique, color_map, centers, offsets, group_width)
    ax.set_xlim(centers[TIERS[0]] - group_width / 2 - GROUP_GAP, centers[TIERS[-1]] + group_width / 2 + GROUP_GAP)
    ax.set_xticks([centers[t] for t in TIERS])
    ax.set_xticklabels([f"Tier {t}" for t in TIERS])
    ax.tick_params(axis='x', length=0)
    ax.set_ylabel("Gene count")
    ax.set_title(f"{args.title}: Genes per Tier by Sample")
    _sample_type_legend(ax, sample_types_unique, color_map, st_total_counts)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    outfile = f"{args.outdir}/tier_gene_counts_density.pdf"
    plt.savefig(outfile)
    plt.close(fig)
    print(f"Saved plot to {outfile}")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sample_type_by_sample = dict(zip(args.samples, args.sample_types))

    df = read_hit_rows(args.infile, sep='\t', keep_default_na=False)
    df = df.astype(object)
    df.fillna('.', inplace=True)

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
    df = df[df['_combo'].apply(len) > 0]

    counts = df.groupby(['sample', '_combo'])['gene'].nunique().unstack(fill_value=0)
    counts = counts.reindex(index=args.samples, fill_value=0)

    if counts.shape[1] == 0:
        print("WARNING: no non-empty category combinations found in the input -- nothing to plot.")
        pd.DataFrame(columns=['sample', 'sample_type', 'combo', 'gene_count']).to_csv(
            f"{args.outdir}/hits_upset_counts.tsv", sep='\t', index=False)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No candidate hits to plot", ha='center', va='center')
        ax.axis('off')
        plt.savefig(f"{args.outdir}/hits_upset_density.pdf")
        plt.close(fig)
        return

    combo_lookup = {frozenset(c): c for c in counts.columns}
    combo_order = [
        combo_lookup[frozenset(combo)]
        for size in range(1, len(CATEGORY_DISPLAY_ORDER) + 1)
        for combo in itertools.combinations(CATEGORY_DISPLAY_ORDER, size)
        if frozenset(combo) in combo_lookup
    ]
    counts = counts[combo_order]

    long_df = counts.reset_index().melt(id_vars='sample', var_name='combo', value_name='gene_count')
    long_df['combo'] = long_df['combo'].apply(_combo_label)
    long_df['sample_type'] = long_df['sample'].map(sample_type_by_sample)
    long_df = long_df[['sample', 'sample_type', 'combo', 'gene_count']]
    long_df.to_csv(f"{args.outdir}/hits_upset_counts.tsv", sep='\t', index=False)
    print(f"Saved counts tsv to {args.outdir}/hits_upset_counts.tsv")

    n_combos = len(combo_order)
    st_total_counts_all = {st: args.sample_types.count(st) for st in set(args.sample_types)}
    sample_types_unique = sorted(st_total_counts_all, key=lambda st: st_total_counts_all[st], reverse=True)
    n_st = len(sample_types_unique)
    st_index = {st: i for i, st in enumerate(sample_types_unique)}
    st_total_counts = st_total_counts_all
    configured = {}
    for st, c in zip(args.sample_types, args.colors or []):
        configured.setdefault(st, c.strip("'\""))
    color_map = _sample_type_color_map(sample_types_unique, configured)

    swarm_df = counts.reset_index().melt(id_vars='sample', var_name='_combo', value_name='gene_count')
    swarm_df['sample_type'] = swarm_df['sample'].map(sample_type_by_sample)

    centers, offsets, group_width = _group_layout(combo_order, sample_types_unique)
    fig, (ax_box, ax_matrix) = plt.subplots(
        2, 1, figsize=(_figure_width(combo_order, group_width), 8), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]},
    )
    _draw_box_strip(ax_box, swarm_df, '_combo', combo_order, sample_types_unique, color_map,
                    centers, offsets, group_width)
    for i, combo in enumerate(combo_order):
        if i % 2 == 1:
            ax_matrix.axvspan(centers[combo] - group_width / 2 - GROUP_GAP / 4,
                              centers[combo] + group_width / 2 + GROUP_GAP / 4, color='#F0F0F0', zorder=0)
    ax_box.set_xlabel('')
    ax_box.set_ylabel("Gene count")
    ax_box.set_title(args.title)
    _sample_type_legend(ax_box, sample_types_unique, color_map, st_total_counts)
    for spine in ('top', 'right'):
        ax_box.spines[spine].set_visible(False)

    n_cats = len(CATEGORIES)
    cat_y = {cat: n_cats - i for i, cat in enumerate(CATEGORIES)}
    for combo in combo_order:
        pos = centers[combo]
        in_combo_y = [cat_y[cat] for cat in CATEGORIES if cat in combo]
        if len(in_combo_y) > 1:
            ax_matrix.plot([pos, pos], [min(in_combo_y), max(in_combo_y)], color='black', lw=1.5, zorder=2)
        for cat in CATEGORIES:
            ax_matrix.scatter([pos], [cat_y[cat]], s=80, zorder=3,
                              color='black' if cat in combo else '#DDDDDD')
    ax_matrix.set_yticks(list(cat_y.values()))
    ax_matrix.set_yticklabels(list(cat_y.keys()))
    ax_matrix.set_ylim(0.5, n_cats + 0.5)
    ax_matrix.set_xticks([centers[c] for c in combo_order])
    ax_matrix.set_xticklabels([])
    ax_matrix.tick_params(axis='x', length=0)
    ax_box.tick_params(axis='x', length=0)
    ax_matrix.set_xlim(centers[combo_order[0]] - group_width / 2 - GROUP_GAP,
                       centers[combo_order[-1]] + group_width / 2 + GROUP_GAP)
    for spine in ('top', 'right', 'bottom'):
        ax_matrix.spines[spine].set_visible(False)

    plt.tight_layout()
    outfile = f"{args.outdir}/hits_upset_density.pdf"
    plt.savefig(outfile)
    plt.close(fig)
    print(f"Saved plot to {outfile}")

    _plot_tier_density(df, args, sample_type_by_sample, sample_types_unique, st_total_counts, color_map)


if __name__ == "__main__":
    main()
