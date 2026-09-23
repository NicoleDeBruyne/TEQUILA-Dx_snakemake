import argparse
import itertools
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde

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

# Per-strip density shading (a smooth KDE over each column's true, unjittered
# values, rendered as stacked translucent rectangles behind the points) --
# opacity is normalized per column across its whole y-axis range, from
# DENSITY_ALPHA_MIN at the least-dense point to DENSITY_ALPHA_MAX at the
# most-dense point.
DENSITY_ALPHA_MIN = 0.05
DENSITY_ALPHA_MAX = 0.40
KDE_GRID_N = 250

# The density band is drawn slightly wider than the strip's own jitter
# spread, so it peeks out from behind the points rather than being exactly
# clipped to them.
DENSITY_WIDTH_FACTOR = 1.15


def _fit_kde(values):
    """A gaussian_kde over a column's true (unjittered) values, falling back
    to a delta-like spike when every value is identical (gaussian_kde can't
    fit a singular covariance)."""
    values = np.asarray(values, dtype=float)
    if np.ptp(values) == 0:
        v0 = values[0]
        return lambda x: np.where(np.isclose(np.asarray(x, dtype=float), v0), 1.0, 0.0)
    kde = gaussian_kde(values, bw_method='scott')
    return lambda x: kde(np.asarray(x, dtype=float))


def _normalize_alpha(dens, dens_min, dens_max):
    """Map raw density values to [DENSITY_ALPHA_MIN, DENSITY_ALPHA_MAX],
    linearly, using the caller's own min/max (computed across that column's
    *entire* y-axis range)."""
    dens = np.asarray(dens, dtype=float)
    if dens_max - dens_min < 1e-12:
        return np.full_like(dens, DENSITY_ALPHA_MAX)
    frac = (dens - dens_min) / (dens_max - dens_min)
    return DENSITY_ALPHA_MIN + frac * (DENSITY_ALPHA_MAX - DENSITY_ALPHA_MIN)


def _draw_density_gradient(ax, color_hex, alphas, x0, x1, y0, y1):
    """A smooth, continuously-interpolated vertical gradient (bilinear
    imshow) instead of a stack of flat-alpha rectangles -- stacked rects
    show a visible band at every boundary; a single interpolated image
    doesn't. Only valid where the y-axis is linear across [y0, y1] (true
    here since the log-scale render was removed and the "shift" bands are
    just constant offsets, not a real nonlinear transform)."""
    rgb = to_rgb(color_hex)
    alphas = np.asarray(alphas, dtype=float)
    img = np.zeros((len(alphas), 1, 4))
    img[:, 0, 0], img[:, 0, 1], img[:, 0, 2] = rgb
    img[:, 0, 3] = alphas
    ax.imshow(
        img, origin='lower', aspect='auto',
        extent=(x0, x1, y0, y1),
        interpolation='bilinear', zorder=0.8,
    )


# Shared "zero band" scheme, used by both plots: every non-zero true value is
# pushed up the y-axis by BAND_GAP, opening up a dedicated gap at the bottom
# of the axis (of height ZERO_BAND_HEIGHT, so from -ZERO_BAND_HEIGHT/2 to
# +ZERO_BAND_HEIGHT/2) that only exact zeros ever land in -- e.g. a true
# value of 1 always displays at y = 1 + BAND_GAP, never inside that band.
BAND_GAP = 2
ZERO_BAND_HEIGHT = 4.0
# Height of the pooled ">cap" catch-all bucket in the combo plot only (the
# tier plot has no catch-all bucket).
CATCHALL_BAND_HEIGHT = 5.0


def _shift(v, gap=BAND_GAP):
    """True value -> display y-position: 0 stays at 0, anything else moves
    up by `gap` to clear the zero band."""
    if np.isscalar(v):
        return v + gap if v > 0 else v
    arr = np.asarray(v, dtype=float)
    return np.where(arr > 0, arr + gap, arr)


def _round_number_ticks(max_val):
    """Clean, evenly-spaced tick candidates for a true-value axis (0, 1, then
    steps of 5/10/20/...), later mapped through _shift for display."""
    step_candidates = [5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000]
    target_ticks = 7
    raw_step = max(max_val, 1) / target_ticks
    step = next((c for c in step_candidates if c >= raw_step), step_candidates[-1])
    return sorted(set([0, 1] + list(range(step, int(max_val) + 1, step))))


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


def _plot_tier_density(df, args, sample_type_by_sample, sample_types_unique, st_total_counts, color_map):
    """Same visual language as the upset plot's density-shaded strips (KDE
    band + median bar, split by sample_type, plus the shared zero-band
    scheme), but keyed on gene tier (1-6) instead of variant/ASE/junction
    combo. No catch-all bucket -- tier counts don't need one."""
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

    swarm_df = long_df.copy()

    MARKER_SIZE = 5
    marker_diameter_in = 2 * np.sqrt((MARKER_SIZE ** 2) / np.pi) / 72

    def _strip_width_factor(n):
        n_clamped = min(max(n, 1), 1000)
        return 3 * (2 ** np.log10(n_clamped))

    SUBCOL_GAP = 0.4
    group_gap = 0.35
    n_st = len(sample_types_unique)

    strip_width_by_st = {st: _strip_width_factor(st_total_counts[st]) * marker_diameter_in for st in sample_types_unique}
    group_width = sum(strip_width_by_st.values()) + (n_st - 1) * SUBCOL_GAP
    group_spacing = group_width + group_gap

    cum = -group_width / 2
    strip_center_by_st = {}
    for st in sample_types_unique:
        strip_center_by_st[st] = cum + strip_width_by_st[st] / 2
        cum += strip_width_by_st[st] + SUBCOL_GAP

    group_center = {tier: i * group_spacing for i, tier in enumerate(TIERS)}

    margin = 0.4
    total_span = (len(TIERS) - 1) * group_spacing + group_width + 2 * margin
    width = min(20, max(6, total_span))

    fig, ax = plt.subplots(figsize=(width, 6))

    for i, tier in enumerate(TIERS):
        if i % 2 == 1:
            band_lo = group_center[tier] - group_width / 2 - group_gap / 4
            band_hi = group_center[tier] + group_width / 2 + group_gap / 4
            ax.axvspan(band_lo, band_hi, color='#F0F0F0', zorder=0)

    rng = np.random.default_rng(0)
    BOX_MIN_N = 5

    ZERO_JITTER = ZERO_BAND_HEIGHT * 0.3
    zero_band_half_height = ZERO_BAND_HEIGHT / 2
    ylim_lo = -zero_band_half_height
    zero_band_top = zero_band_half_height

    max_val = max(swarm_df['gene_count'].max(), 1)
    ylim_hi = _shift(max_val) + 0.5

    local_swarm_df = swarm_df.copy()
    strip_x = pd.Series(index=local_swarm_df.index, dtype=float)
    strip_y = pd.Series(index=local_swarm_df.index, dtype=float)
    for (tier, st), g in local_swarm_df.groupby(['tier', 'sample_type']):
        x_center = group_center[tier] + strip_center_by_st[st]
        half_width = strip_width_by_st[st] * 0.9 / 2
        n = len(g)
        strip_x.loc[g.index] = x_center + rng.uniform(-half_width, half_width, size=n)
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

    for (tier, st), g in local_swarm_df.groupby(['tier', 'sample_type']):
        if len(g) < BOX_MIN_N:
            continue
        x_center = group_center[tier] + strip_center_by_st[st]
        half_width = strip_width_by_st[st] * 0.9 / 2
        x0, x1 = x_center - half_width, x_center + half_width
        dens_half_width = half_width * DENSITY_WIDTH_FACTOR
        dx0, dx1 = x_center - dens_half_width, x_center + dens_half_width
        vals = g['gene_count'].to_numpy(dtype=float)
        kde = _fit_kde(vals)

        # Zero band (a single point mass at 0) + the continuous region above
        # it, normalized together as one column -- same scheme as the upset
        # plot, minus its catch-all bucket.
        main_grid = np.linspace(1e-6, max_val, KDE_GRID_N)
        dens_main = kde(main_grid)
        dens_zero = float(kde(np.array([0.0]))[0])
        dens_all = np.concatenate([dens_main, [dens_zero]])
        dens_min, dens_max = dens_all.min(), dens_all.max()

        alpha_zero = float(_normalize_alpha(np.array([dens_zero]), dens_min, dens_max)[0])
        ax.add_patch(Rectangle(
            (dx0, ylim_lo), dx1 - dx0, zero_band_top - ylim_lo,
            facecolor=color_map[st], edgecolor='none', alpha=alpha_zero, zorder=0.8,
        ))

        main_grid_display = _shift(main_grid)
        alphas_main = _normalize_alpha(dens_main, dens_min, dens_max)
        _draw_density_gradient(
            ax, color_map[st], alphas_main,
            dx0, dx1, main_grid_display[0], main_grid_display[-1],
        )

        med_display = _shift(np.median(vals))
        ax.plot([x0, x1], [med_display, med_display], color=color_map[st], linewidth=2.2, zorder=1.5)

    for st in sample_types_unique:
        sub = local_swarm_df[local_swarm_df['sample_type'] == st]
        zero_mask = sub['gene_count'] == 0
        ax.scatter(
            sub.loc[zero_mask, 'strip_x'], sub.loc[zero_mask, 'strip_y'],
            facecolors='none', edgecolors=color_map[st],
            s=MARKER_SIZE ** 2, alpha=0.6, linewidth=1.0, zorder=2,
        )
        ax.scatter(
            sub.loc[~zero_mask, 'strip_x'], sub.loc[~zero_mask, 'strip_y'],
            color=color_map[st],
            s=MARKER_SIZE ** 2, alpha=0.6, linewidth=0, zorder=2,
        )

    ax.axhspan(ylim_lo, zero_band_top, color='#D5D5D5', zorder=0.5)
    ax.axhline(zero_band_top, color='#888888', lw=1.0, zorder=0.6)

    ax.set_ylim(ylim_lo, ylim_hi)
    ax.set_xlim(
        group_center[TIERS[0]] - group_width / 2 - margin,
        group_center[TIERS[-1]] + group_width / 2 + margin,
    )
    ax.set_xticks([group_center[t] for t in TIERS])
    ax.set_xticklabels([f"Tier {t}" for t in TIERS])
    ax.set_xlabel('')
    ax.set_ylabel("Gene count")
    ax.set_title(f"{args.title}: Genes per Tier by Sample")

    true_ticks = _round_number_ticks(max_val)
    ax.set_yticks([_shift(t) for t in true_ticks])
    ax.set_yticklabels([str(t) for t in true_ticks])

    handles = [
        plt.Line2D([0], [0], marker='o', linestyle='', color=color_map[st],
                   label=f"{st} (n={st_total_counts[st]})", markersize=7)
        for st in sample_types_unique
    ]
    ax.legend(handles=handles, title="sample_type", bbox_to_anchor=(1.02, 1), loc="upper left")
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

    SUBCOL_GAP = 0.4
    group_gap = 0.35

    # One column per sample_type now -- density shading + median bar are
    # drawn behind/at the strip itself, so there's no separate box column.
    strip_width_by_st = {st: _strip_width_factor(st_total_counts[st]) * marker_diameter_in for st in sample_types_unique}
    col_width_by_st = strip_width_by_st

    group_width = sum(col_width_by_st.values()) + (n_st - 1) * SUBCOL_GAP
    group_spacing = group_width + group_gap

    cum = -group_width / 2
    strip_center_by_st = {}
    for st in sample_types_unique:
        strip_center_by_st[st] = cum + strip_width_by_st[st] / 2
        cum += col_width_by_st[st] + SUBCOL_GAP

    group_center = {combo: i * group_spacing for i, combo in enumerate(combo_order)}

    margin = 0.4
    total_span = (n_combos - 1) * group_spacing + group_width + 2 * margin
    width = min(20, max(6, total_span))

    def _render(outfile):
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

            LINEAR_CAP = None
            if swarm_df['gene_count'].max() >= 50:
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
                CATCHALL_Y = CATCHALL_BOUNDARY + CATCHALL_BAND_HEIGHT / 2
                CATCHALL_JITTER = CATCHALL_BAND_HEIGHT * 0.3

            ZERO_JITTER = ZERO_BAND_HEIGHT * 0.3

            # Axis-band boundaries, computed up front (rather than down in the
            # tick-setting code below) so the density-shading loop can place
            # its rectangles against the same boundaries the axis itself uses.
            zero_band_half_height = ZERO_BAND_HEIGHT / 2
            ylim_lo = -zero_band_half_height
            zero_band_top = zero_band_half_height
            if LINEAR_CAP is not None:
                catchall_boundary = _shift(CATCHALL_BOUNDARY)
                catchall_y = _shift(CATCHALL_Y)
                ylim_hi = _shift(CATCHALL_BOUNDARY + CATCHALL_BAND_HEIGHT)
            else:
                ylim_hi = _shift(swarm_df['gene_count'].max()) + 0.5

            local_swarm_df = swarm_df.copy()
            strip_x = pd.Series(index=local_swarm_df.index, dtype=float)
            strip_y = pd.Series(index=local_swarm_df.index, dtype=float)
            for (combo, st), g in local_swarm_df.groupby(['_combo', 'sample_type']):
                x_center = group_center[combo] + strip_center_by_st[st]
                half_width = strip_width_by_st[st] * 0.9 / 2
                n = len(g)
                strip_x.loc[g.index] = x_center + rng.uniform(-half_width, half_width, size=n)
                if LINEAR_CAP is not None:
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
                x_center = group_center[combo] + strip_center_by_st[st]
                half_width = strip_width_by_st[st] * 0.9 / 2
                x0, x1 = x_center - half_width, x_center + half_width
                dens_half_width = half_width * DENSITY_WIDTH_FACTOR
                dx0, dx1 = x_center - dens_half_width, x_center + dens_half_width
                vals = g['gene_count'].to_numpy(dtype=float)
                kde = _fit_kde(vals)

                # Three sub-regions, normalized together as one column: the
                # zero band (a single point mass at 0), the continuous region
                # up to LINEAR_CAP (or the true max if there's no cap), and --
                # if present -- the pooled catch-all bucket for values above
                # LINEAR_CAP.
                cap = LINEAR_CAP if LINEAR_CAP is not None else max(swarm_df['gene_count'].max(), 1)
                main_grid = np.linspace(1e-6, cap, KDE_GRID_N)
                dens_main = kde(main_grid)
                dens_zero = float(kde(np.array([0.0]))[0])
                dens_all = np.concatenate([dens_main, [dens_zero]])

                dens_tail = None
                if LINEAR_CAP is not None:
                    tail_vals = vals[vals > LINEAR_CAP]
                    if len(tail_vals) > 0:
                        dens_tail = float(np.mean(kde(tail_vals)))
                        dens_all = np.concatenate([dens_all, [dens_tail]])

                dens_min, dens_max = dens_all.min(), dens_all.max()

                # Zero band: flat shading over the whole gray zero-band rectangle.
                alpha_zero = float(_normalize_alpha(np.array([dens_zero]), dens_min, dens_max)[0])
                ax_box.add_patch(Rectangle(
                    (dx0, ylim_lo), dx1 - dx0, zero_band_top - ylim_lo,
                    facecolor=color_map[st], edgecolor='none', alpha=alpha_zero, zorder=0.8,
                ))

                # Main region: one smooth interpolated gradient, mapped
                # through the same shift the strip's own points use (a
                # constant offset, so the grid stays evenly spaced).
                main_grid_display = _shift(main_grid)
                alphas_main = _normalize_alpha(dens_main, dens_min, dens_max)
                _draw_density_gradient(
                    ax_box, color_map[st], alphas_main,
                    dx0, dx1, main_grid_display[0], main_grid_display[-1],
                )

                # Catch-all bucket: flat shading over the whole >cap rectangle.
                if dens_tail is not None:
                    alpha_tail = float(_normalize_alpha(np.array([dens_tail]), dens_min, dens_max)[0])
                    ax_box.add_patch(Rectangle(
                        (dx0, catchall_boundary), dx1 - dx0, ylim_hi - catchall_boundary,
                        facecolor=color_map[st], edgecolor='none', alpha=alpha_tail, zorder=0.8,
                    ))

                med_true = np.median(vals)
                med_display = _shift(min(med_true, cap)) if LINEAR_CAP is not None else _shift(med_true)

                # Single median bar, spanning the strip's own width (not a
                # separate fixed box width).
                ax_box.plot([x0, x1], [med_display, med_display],
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

            ax_box.axhspan(ylim_lo, zero_band_top, color='#D5D5D5', zorder=0.5)
            ax_box.axhline(zero_band_top, color='#888888', lw=1.0, zorder=0.6)

            if LINEAR_CAP is not None:
                ax_box.set_ylim(ylim_lo, ylim_hi)

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
        ax_box.set_title(args.title)

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

    _render(outfile=f"{args.outdir}/hits_upset_density.pdf")

    _plot_tier_density(df, args, sample_type_by_sample, sample_types_unique, st_total_counts, color_map)


if __name__ == "__main__":
    main()