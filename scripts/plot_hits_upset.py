

import argparse
import itertools
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, NullFormatter

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

_SAMPLE_TYPE_COLORS = [
    '#6997B9',
    '#BB6A68',
    '#70A677',
    '#D48653',
    '#A783A3',
]


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
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--title", default="Candidate Gene Hit Categories by Sample", help="Plot title.")
    args = parser.parse_args()
    if len(args.samples) != len(args.sample_types):
        parser.error("--samples and --sample-types must have the same number of entries")
    return args


def _sample_type_color_map(sample_types_in_order):
    seen = list(dict.fromkeys(sample_types_in_order))
    return {st: _SAMPLE_TYPE_COLORS[i % len(_SAMPLE_TYPE_COLORS)] for i, st in enumerate(seen)}


def _combo_label(combo):
    return " & ".join(combo) if combo else "(none)"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sample_type_by_sample = dict(zip(args.samples, args.sample_types))

    df = pd.read_csv(args.infile, sep='\t', keep_default_na=False)
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
    color_map = _sample_type_color_map(sample_types_unique)

    swarm_df = counts.reset_index().melt(id_vars='sample', var_name='_combo', value_name='gene_count')
    swarm_df['sample_type'] = swarm_df['sample'].map(sample_type_by_sample)

    MARKER_SIZE = 5
    marker_diameter_in = 2 * np.sqrt((MARKER_SIZE ** 2) / np.pi) / 72

    def _strip_width_factor(n):
        n_clamped = min(max(n, 1), 1000)
        return 3 * (2 ** np.log10(n_clamped))

    BOX_WIDTH = 0.25
    BOX_STRIP_GAP = 0.06
    SUBCOL_GAP = 0.4
    group_gap = 0.35

    strip_width_by_st = {st: _strip_width_factor(st_total_counts[st]) * marker_diameter_in for st in sample_types_unique}
    col_width_by_st = {st: BOX_WIDTH + BOX_STRIP_GAP + strip_width_by_st[st] for st in sample_types_unique}

    group_width = sum(col_width_by_st.values()) + (n_st - 1) * SUBCOL_GAP
    group_spacing = group_width + group_gap

    cum = -group_width / 2
    box_left_by_st = {}
    strip_center_by_st = {}
    for st in sample_types_unique:
        box_left_by_st[st] = cum
        strip_center_by_st[st] = cum + BOX_WIDTH + BOX_STRIP_GAP + strip_width_by_st[st] / 2
        cum += col_width_by_st[st] + SUBCOL_GAP

    group_center = {combo: i * group_spacing for i, combo in enumerate(combo_order)}

    margin = 0.4
    total_span = (n_combos - 1) * group_spacing + group_width + 2 * margin
    width = min(20, max(6, total_span))

    def _render(log_scale, outfile):
        fig, (ax_box, ax_matrix) = plt.subplots(
            2, 1, figsize=(width, 8), sharex=True,
            gridspec_kw={'height_ratios': [3, 1]},
        )

        for i, combo in enumerate(combo_order):
            if i % 2 == 1:
                band_lo = group_center[combo] - group_width / 2 - group_gap / 4
                band_hi = group_center[combo] + group_width / 2 + group_gap / 4
                ax_box.axvspan(band_lo, band_hi, color='#F0F0F0', zorder=0)
                ax_matrix.axvspan(band_lo, band_hi, color='#F0F0F0', zorder=0)

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*cannot be placed.*")

            rng = np.random.default_rng(0)

            BOX_MIN_N = 5


            BAND_GAP = 2 if not log_scale else None
            BAND_HEIGHT = 5.0 if not log_scale else None

            LINEAR_CAP = None
            if not log_scale and swarm_df['gene_count'].max() >= 50:
                per_strip_cutoffs = []
                for _, g in swarm_df.groupby(['_combo', 'sample_type']):
                    vals_desc = np.sort(g['gene_count'].to_numpy(dtype=float))[::-1]
                    n = len(vals_desc)
                    allowed_above = int(np.floor(0.05 * n))
                    idx = min(allowed_above, n - 1)
                    per_strip_cutoffs.append(vals_desc[idx])
                raw_cap = max(per_strip_cutoffs)
                LINEAR_CAP = max(10, int(np.ceil(raw_cap / 10.0) * 10))
                CATCHALL_BOUNDARY = LINEAR_CAP + BAND_GAP
                CATCHALL_Y = CATCHALL_BOUNDARY + BAND_HEIGHT / 2
                CATCHALL_JITTER = BAND_HEIGHT * 0.3

            ZERO_GAP = BAND_GAP
            ZERO_JITTER = BAND_HEIGHT * 0.3 if not log_scale else None

            def _shift(v):
                if np.isscalar(v):
                    return v + ZERO_GAP if v > 0 else v
                arr = np.asarray(v, dtype=float)
                return np.where(arr > 0, arr + ZERO_GAP, arr)

            local_swarm_df = swarm_df.copy()
            strip_x = pd.Series(index=local_swarm_df.index, dtype=float)
            strip_y = pd.Series(index=local_swarm_df.index, dtype=float)
            for (combo, st), g in local_swarm_df.groupby(['_combo', 'sample_type']):
                x_center = group_center[combo] + strip_center_by_st[st]
                half_width = strip_width_by_st[st] * 0.9 / 2
                n = len(g)
                strip_x.loc[g.index] = x_center + rng.uniform(-half_width, half_width, size=n)
                if log_scale:
                    base = g['gene_count'].to_numpy(dtype=float)
                    additive = base + rng.uniform(-0.15, 0.15, size=n)
                    multiplicative = base * (10 ** rng.uniform(-0.04, 0.04, size=n))
                    strip_y.loc[g.index] = np.where(base <= 1, additive, multiplicative)
                elif LINEAR_CAP is not None:
                    raw = g['gene_count'].to_numpy(dtype=float)
                    over = raw > LINEAR_CAP
                    zero = raw == 0
                    y_base = np.where(over, CATCHALL_Y, raw)
                    y_base = _shift(y_base)
                    jitter = np.where(
                        over,
                        rng.uniform(-CATCHALL_JITTER, CATCHALL_JITTER, size=n),
                        np.where(
                            zero,
                            rng.uniform(-ZERO_JITTER, ZERO_JITTER, size=n),
                            rng.uniform(-0.15, 0.15, size=n),
                        ),
                    )
                    strip_y.loc[g.index] = y_base + jitter
                else:
                    raw = g['gene_count'].to_numpy(dtype=float)
                    zero = raw == 0
                    jitter = np.where(
                        zero,
                        rng.uniform(-ZERO_JITTER, ZERO_JITTER, size=n),
                        rng.uniform(-0.15, 0.15, size=n),
                    )
                    strip_y.loc[g.index] = _shift(raw) + jitter
            local_swarm_df['strip_x'] = strip_x
            local_swarm_df['strip_y'] = strip_y

            for (combo, st), g in local_swarm_df.groupby(['_combo', 'sample_type']):
                if len(g) < BOX_MIN_N:
                    continue
                box_left = group_center[combo] + box_left_by_st[st]
                vals = g['gene_count'].to_numpy(dtype=float)
                q1, med, q3 = np.percentile(vals, [25, 50, 75])
                if not log_scale:
                    if LINEAR_CAP is not None:
                        q1, med, q3 = (min(v, LINEAR_CAP) for v in (q1, med, q3))
                    q1, med, q3 = (_shift(v) for v in (q1, med, q3))
                ax_box.add_patch(Rectangle(
                    (box_left, q1), BOX_WIDTH, max(q3 - q1, 0),
                    facecolor=color_map[st], alpha=0.25, edgecolor=color_map[st],
                    linewidth=1.3, zorder=1,
                ))
                ax_box.plot([box_left, box_left + BOX_WIDTH], [med, med],
                            color=color_map[st], linewidth=2.2, zorder=1.5)

            for st in sample_types_unique:
                sub = local_swarm_df[local_swarm_df['sample_type'] == st]
                zero_mask = sub['gene_count'] == 0
                ax_box.scatter(
                    sub.loc[zero_mask, 'strip_x'], sub.loc[zero_mask, 'strip_y'],
                    facecolors='none', edgecolors=color_map[st],
                    s=MARKER_SIZE ** 2, alpha=0.6, linewidth=1.0, zorder=2,
                )
                ax_box.scatter(
                    sub.loc[~zero_mask, 'strip_x'], sub.loc[~zero_mask, 'strip_y'],
                    color=color_map[st],
                    s=MARKER_SIZE ** 2, alpha=0.6, linewidth=0, zorder=2,
                )

            if log_scale:
                ax_box.set_yscale('symlog', linthresh=1)
                ylim_lo = -0.5
                ax_box.set_ylim(ylim_lo, swarm_df['gene_count'].max() * 1.4)

                max_val = swarm_df['gene_count'].max()
                n_decades = int(np.floor(np.log10(max(max_val, 1)))) + 1
                decade_ticks = [10 ** i for i in range(n_decades)]
                ax_box.set_yticks([0] + decade_ticks)
                ax_box.set_yticklabels(['0'] + [str(t) for t in decade_ticks])
                ax_box.yaxis.set_minor_locator(LogLocator(base=10, subs=[2, 4, 6, 8]))
                ax_box.yaxis.set_minor_formatter(NullFormatter())

                ax_box.axhspan(ylim_lo, 0.5, color='#D5D5D5', zorder=0.5)
                ax_box.axhline(0.5, color='#888888', lw=1.0, zorder=0.6)
            else:
                zero_band_half_height = BAND_HEIGHT / 2
                ylim_lo = -zero_band_half_height

                zero_band_top = zero_band_half_height
                ax_box.axhspan(ylim_lo, zero_band_top, color='#D5D5D5', zorder=0.5)
                ax_box.axhline(zero_band_top, color='#888888', lw=1.0, zorder=0.6)


                def _round_number_ticks(max_val):
                    step_candidates = [5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000]
                    target_ticks = 7
                    raw_step = max(max_val, 1) / target_ticks
                    step = next((c for c in step_candidates if c >= raw_step), step_candidates[-1])
                    return sorted(set([0, 1] + list(range(step, int(max_val) + 1, step))))

                if LINEAR_CAP is not None:
                    ylim_hi = _shift(CATCHALL_BOUNDARY + BAND_HEIGHT)
                    ax_box.set_ylim(ylim_lo, ylim_hi)

                    catchall_boundary = _shift(CATCHALL_BOUNDARY)
                    catchall_y = _shift(CATCHALL_Y)
                    ax_box.axhspan(catchall_boundary, ylim_hi, color='#D5D5D5', zorder=0.5)
                    ax_box.axhline(catchall_boundary, color='#888888', lw=1.0, zorder=0.6)
                    ax_box.text(
                        0.01, catchall_y, f'>{LINEAR_CAP}', color='#555555', fontsize=10,
                        ha='left', va='center', zorder=4, transform=ax_box.get_yaxis_transform(),
                    )

                    true_ticks = _round_number_ticks(LINEAR_CAP)
                    ax_box.set_yticks([_shift(t) for t in true_ticks])
                    ax_box.set_yticklabels([str(t) for t in true_ticks])
                else:
                    max_val = swarm_df['gene_count'].max()
                    ax_box.set_ylim(ylim_lo, _shift(max_val) + 0.5)

                    true_ticks = _round_number_ticks(max_val)
                    ax_box.set_yticks([_shift(t) for t in true_ticks])
                    ax_box.set_yticklabels([str(t) for t in true_ticks])


        ax_box.set_xlabel('')
        ax_box.set_ylabel("Gene count")
        ax_box.set_title(args.title + ("  (log scale)" if log_scale else ""))

        handles = [
            plt.Line2D([0], [0], marker='o', linestyle='', color=color_map[st],
                       label=f"{st} (n={st_total_counts[st]})", markersize=7)
            for st in sample_types_unique
        ]
        ax_box.legend(handles=handles, title="sample_type", bbox_to_anchor=(1.02, 1), loc="upper left")

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
        ax_matrix.tick_params(axis='x', length=0)
        ax_box.tick_params(axis='x', length=0)
        xlim_lo = group_center[combo_order[0]] - group_width / 2 - margin
        xlim_hi = group_center[combo_order[-1]] + group_width / 2 + margin
        ax_matrix.set_xlim(xlim_lo, xlim_hi)
        for spine in ('top', 'right', 'bottom'):
            ax_matrix.spines[spine].set_visible(False)

        plt.tight_layout()
        plt.savefig(outfile)
        plt.close(fig)
        print(f"Saved plot to {outfile}")

    _render(log_scale=False, outfile=f"{args.outdir}/hits_upset_density.pdf")
    _render(log_scale=True, outfile=f"{args.outdir}/hits_upset_density_log.pdf")


if __name__ == "__main__":
    main()