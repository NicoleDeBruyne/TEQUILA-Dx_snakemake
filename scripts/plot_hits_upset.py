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
# each combination: a basic jittered strip plot of individual points, with
# a boxplot (just Q1-Q3 and the median, no whiskers) underlaid behind it.
# Each combination gets its own group of side-by-side sub-columns, one per
# sample_type -- every sample_type gets the SAME amount of horizontal
# space, regardless of its cohort size -- rather than mixing every
# sample_type into a single column. Alternating background shading marks
# each combination's group, and the y=0 region is shaded a darker gray so
# it stays immediately recognizable as the zero baseline. The
# combination-membership matrix is drawn as a second panel below, centered
# under each combination's group, in the usual UpSet style.

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

# Fixed left-to-right display order for combos on the plot: singles in this
# category order, then pairs, then the triple, always -- e.g. variant,
# junction, ASE, variant & junction, variant & ASE, junction & ASE,
# variant & junction & ASE. Independent of CATEGORIES above (which only
# controls the internal tuple representation of each combo) and NOT
# data-driven -- this order is used regardless of each combo's size/median
# in the data.
CATEGORY_DISPLAY_ORDER = ['variant', 'junction', 'ASE']

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
        # Still write an (empty) counts file and (empty) placeholder PDF so
        # this rule's declared outputs exist.
        pd.DataFrame(columns=['sample', 'sample_type', 'combo', 'gene_count']).to_csv(
            f"{args.outdir}/hits_upset_counts.tsv", sep='\t', index=False)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No candidate hits to plot", ha='center', va='center')
        ax.axis('off')
        plt.savefig(f"{args.outdir}/hits_upset_density.pdf")
        plt.close(fig)
        return

    # Fixed display order (see CATEGORY_DISPLAY_ORDER above): singles, then
    # pairs, then the triple, always in the same category priority --
    # never data-driven. Combos are matched by set membership since a
    # combo's actual tuple (a column key in `counts`) follows CATEGORIES'
    # own internal order, which may differ from the display order.
    combo_lookup = {frozenset(c): c for c in counts.columns}
    combo_order = [
        combo_lookup[frozenset(combo)]
        for size in range(1, len(CATEGORY_DISPLAY_ORDER) + 1)
        for combo in itertools.combinations(CATEGORY_DISPLAY_ORDER, size)
        if frozenset(combo) in combo_lookup
    ]
    counts = counts[combo_order]

    # Long-format counts tsv, for reference/debugging alongside the plot.
    long_df = counts.reset_index().melt(id_vars='sample', var_name='combo', value_name='gene_count')
    long_df['combo'] = long_df['combo'].apply(_combo_label)
    long_df['sample_type'] = long_df['sample'].map(sample_type_by_sample)
    long_df = long_df[['sample', 'sample_type', 'combo', 'gene_count']]
    long_df.to_csv(f"{args.outdir}/hits_upset_counts.tsv", sep='\t', index=False)
    print(f"Saved counts tsv to {args.outdir}/hits_upset_counts.tsv")

    # ---------- Layout: boxplot (fixed width) + strip (width scales with n), side by side ----------
    n_combos = len(combo_order)
    # Ordered by cohort size (largest first), not first-appearance in
    # --sample-types -- both the x-position order within each combo group
    # and the color assignment follow this order, so the biggest
    # sample_type is always leftmost/first-colored.
    st_total_counts_all = {st: args.sample_types.count(st) for st in set(args.sample_types)}
    sample_types_unique = sorted(st_total_counts_all, key=lambda st: st_total_counts_all[st], reverse=True)
    n_st = len(sample_types_unique)
    st_index = {st: i for i, st in enumerate(sample_types_unique)}
    st_total_counts = st_total_counts_all
    color_map = _sample_type_color_map(sample_types_unique)

    swarm_df = counts.reset_index().melt(id_vars='sample', var_name='_combo', value_name='gene_count')
    swarm_df['sample_type'] = swarm_df['sample'].map(sample_type_by_sample)

    MARKER_SIZE = 5  # basis for s=MARKER_SIZE**2 in scatter(); also used below to size the strip's allocated width
    marker_diameter_in = 2 * np.sqrt((MARKER_SIZE ** 2) / np.pi) / 72  # approx marker width, converted points -> inches

    # Strip width, in multiples of the marker's own width: 3x at n=1, and
    # doubling every time n increases by a factor of 10 (so 6x at n=10, 12x
    # at n=100, 24x at the n=1000 cap) -- i.e. 3 * 2**log10(n), n clamped to
    # [1, 1000] so a single lonely sample still gets a usable sliver and an
    # enormous cohort doesn't grow unbounded.
    def _strip_width_factor(n):
        n_clamped = min(max(n, 1), 1000)
        return 3 * (2 ** np.log10(n_clamped))

    BOX_WIDTH = 0.25       # fixed width for every boxplot, regardless of sample_type or n
    BOX_STRIP_GAP = 0.06   # gap between a sample_type's box and the start of its strip
    SUBCOL_GAP = 0.4       # inches between adjacent sample_type columns within a combo group
    group_gap = 0.35       # inches between adjacent combos' groups

    strip_width_by_st = {st: _strip_width_factor(st_total_counts[st]) * marker_diameter_in for st in sample_types_unique}
    col_width_by_st = {st: BOX_WIDTH + BOX_STRIP_GAP + strip_width_by_st[st] for st in sample_types_unique}

    group_width = sum(col_width_by_st.values()) + (n_st - 1) * SUBCOL_GAP
    group_spacing = group_width + group_gap

    # Each sample_type's column is boxplot-then-strip, left to right, placed
    # by cumulative width (not an even split) -- box_left_by_st is that
    # column's box's left edge, strip_center_by_st is its strip's jitter
    # center, both relative to the group's own center.
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

        # Faint alternating background bands, one per combo group, so
        # neighboring combos' groups read as visually distinct clusters even
        # though there's no per-sub-column text label.
        for i, combo in enumerate(combo_order):
            if i % 2 == 1:
                band_lo = group_center[combo] - group_width / 2 - group_gap / 4
                band_hi = group_center[combo] + group_width / 2 + group_gap / 4
                ax_box.axvspan(band_lo, band_hi, color='#F0F0F0', zorder=0)
                ax_matrix.axvspan(band_lo, band_hi, color='#F0F0F0', zorder=0)

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*cannot be placed.*")

            # Seeded fresh (not shared across the two _render calls) so the
            # linear and log-scaled versions get IDENTICAL jitter -- makes
            # the two files directly comparable rather than each having its
            # own independent random draw.
            rng = np.random.default_rng(0)

            BOX_MIN_N = 5     # minimum (dense, 0s included) sample count to draw a box at all

            # Adaptive top cap for the LINEAR render only -- if every value
            # in the whole dataset is under 50, no cap is needed at all and
            # the axis just scales to the real max as before. Otherwise,
            # find the smallest threshold such that NO (combo, sample_type)
            # strip has more than 5% of its own points above it (a strip's
            # own 5% cutoff is its value at rank floor(0.05*n) from the top;
            # the shared cap has to be at least as high as the largest such
            # cutoff across every strip, so no single strip ends up with
            # more than its allotted 5% pushed into the catch-all), then
            # round that up to the nearest multiple of 10 (rounding up, not
            # to nearest, so the 5%-or-fewer guarantee still holds after
            # rounding). Not applied to the log render, which already
            # handles wide ranges natively without needing a cap.

            # Shared sizing for the 0-band and the catch-all band, so the
            # two always match: BAND_GAP is the empty space between the
            # real data and each special zone (0<->1, cap<->catch-all);
            # BAND_HEIGHT is the total height of each zone itself. Not
            # used on the log render, where symlog gives 0 its own
            # honestly-scaled linear region natively and there's no cap/
            # catch-all at all.
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
                # Catch-all row sits a FIXED distance above the cap (not
                # scaled proportionally to it), matching the 0-band's fixed
                # size below -- both read as their own separate, consistently-
                # sized zone regardless of how big the cap itself is, and
                # both use the exact same BAND_GAP/BAND_HEIGHT as the 0-band.
                CATCHALL_BOUNDARY = LINEAR_CAP + BAND_GAP
                CATCHALL_Y = CATCHALL_BOUNDARY + BAND_HEIGHT / 2
                CATCHALL_JITTER = BAND_HEIGHT * 0.3

            # Extra visual gap pushed between 0 and the rest of the linear
            # scale: every nonzero value is shifted up by ZERO_GAP before
            # plotting (0 itself is untouched), so 0-valued points sit in
            # visibly empty space below the rest of the data instead of
            # blending into small nonzero values. Same size as the gap
            # above the catch-all band (BAND_GAP), so both gaps read as
            # the same fixed distance.
            ZERO_GAP = BAND_GAP
            # 0-valued points get their own, much wider vertical jitter
            # (rather than the same +/-0.15 used for real values) so a
            # dense cluster of them spreads out across the 0-band instead
            # of piling on top of each other into a solid blob. Kept
            # inside the 0-band's own height (BAND_HEIGHT), same margin
            # convention as the catch-all band's jitter above.
            ZERO_JITTER = BAND_HEIGHT * 0.3 if not log_scale else None

            def _shift(v):
                """True value(s) -> plotted y position, per the ZERO_GAP scheme above. Works on scalars or arrays."""
                if np.isscalar(v):
                    return v + ZERO_GAP if v > 0 else v
                arr = np.asarray(v, dtype=float)
                return np.where(arr > 0, arr + ZERO_GAP, arr)

            local_swarm_df = swarm_df.copy()
            strip_x = pd.Series(index=local_swarm_df.index, dtype=float)
            strip_y = pd.Series(index=local_swarm_df.index, dtype=float)
            for (combo, st), g in local_swarm_df.groupby(['_combo', 'sample_type']):
                x_center = group_center[combo] + strip_center_by_st[st]
                half_width = strip_width_by_st[st] * 0.9 / 2  # small margin so points don't sit flush against the edge
                n = len(g)
                strip_x.loc[g.index] = x_center + rng.uniform(-half_width, half_width, size=n)
                if log_scale:
                    # symlog handles 0 natively (it's a real, linear region
                    # near zero, not a log-transformed one), so gene_count
                    # is used as-is -- no placeholder substitution needed.
                    # Jitter is hybrid: additive (like the linear plot) for
                    # values inside the linear region (<=1, matching
                    # linthresh below), since an additive bump makes sense
                    # on a locally-linear scale; multiplicative (log-space)
                    # beyond that, since a fixed additive +/-0.15 would look
                    # huge relative to a value of 1 but invisible relative
                    # to a value of 100.
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

            # Boxplot, positioned to the LEFT of its strip (not underneath
            # it) -- just Q1-Q3 with a median line, no whiskers, on the
            # TRUE gene_count (0s included, no capping). On the log render
            # (symlog, see below) 0 is a real, directly-plottable position,
            # so no display floor/placeholder is needed there anymore; on
            # the linear render, if an adaptive cap is active, any bound
            # above it is capped to LINEAR_CAP, matching the points' own
            # catch-all treatment. Every box is the same fixed BOX_WIDTH
            # regardless of sample_type/n -- only the strip's width varies
            # with n, not the box's. Requires >=BOX_MIN_N samples in that
            # (combo, sample_type) sub-group (dense, 0s included) to draw at
            # all -- below that, points only.
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

            # Strip points, on top of the box underlay. 0-valued points are
            # drawn hollow (outline only, no fill) so they read as
            # categorically distinct from real nonzero counts, not just
            # smaller/fainter versions of them.
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
                # symlog rather than plain log: below linthresh=1 the axis
                # is LINEAR (0 is a real, honestly-positioned value there,
                # not a fabricated placeholder), and above it the axis is
                # log-compressed same as before. This replaces the old
                # LOG_ZERO_Y placeholder hack -- that approach made a box
                # or line spanning "0 vs 1" look like it spanned a full
                # order of magnitude, which is a real distortion (0 isn't
                # "a bit less than 0.1", it's categorically different from
                # any positive value, and log distance to it is undefined).
                ax_box.set_yscale('symlog', linthresh=1)
                ylim_lo = -0.5
                ax_box.set_ylim(ylim_lo, swarm_df['gene_count'].max() * 1.4)

                # Major ticks at 0/1/10/100/1000 (only as many as the
                # visible range actually needs) -- 0 is now a real tick, not
                # a relabeled placeholder. Minor ticks at 2,4,6,8 / 20,40,60,80
                # / etc with no labels, in the log-compressed region above
                # linthresh -- standard log-scale tick convention.
                max_val = swarm_df['gene_count'].max()
                n_decades = int(np.floor(np.log10(max(max_val, 1)))) + 1
                decade_ticks = [10 ** i for i in range(n_decades)]
                ax_box.set_yticks([0] + decade_ticks)
                ax_box.set_yticklabels(['0'] + [str(t) for t in decade_ticks])
                ax_box.yaxis.set_minor_locator(LogLocator(base=10, subs=[2, 4, 6, 8]))
                ax_box.yaxis.set_minor_formatter(NullFormatter())

                # Same "shade the 0 baseline" treatment as the linear plot
                # below -- 0 is real here too now, so the same convention
                # applies directly instead of needing a placeholder-specific
                # version of it.
                ax_box.axhspan(ylim_lo, 0.5, color='#D5D5D5', zorder=0.5)
                ax_box.axhline(0.5, color='#888888', lw=1.0, zorder=0.6)
            else:
                # 0-band spans the same total height (BAND_HEIGHT) as the
                # catch-all band above, split evenly above/below y=0 --
                # ZERO_JITTER (already sized to 0.3*BAND_HEIGHT, same
                # margin convention as the catch-all band's own jitter)
                # comfortably fits inside it with room to spare.
                zero_band_half_height = BAND_HEIGHT / 2
                ylim_lo = -zero_band_half_height

                # Shade up to just below where the shifted nonzero scale
                # begins -- the gap itself reads as visible empty space and
                # 0-valued points are unmistakably separated from the rest
                # of the scale rather than blending into small nonzero
                # values.
                zero_band_top = zero_band_half_height
                ax_box.axhspan(ylim_lo, zero_band_top, color='#D5D5D5', zorder=0.5)
                ax_box.axhline(zero_band_top, color='#888888', lw=1.0, zorder=0.6)


                def _round_number_ticks(max_val):
                    """0, 1, then evenly-spaced multiples of a 'nice' step (5, 10,
                    20, 25, 50, 100, ...) up to max_val -- e.g. 0,1,5,10,15,20 or
                    0,1,10,20,30 -- rather than MaxNLocator's arbitrary integer
                    steps."""
                    step_candidates = [5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000]
                    target_ticks = 7
                    raw_step = max(max_val, 1) / target_ticks
                    step = next((c for c in step_candidates if c >= raw_step), step_candidates[-1])
                    return sorted(set([0, 1] + list(range(step, int(max_val) + 1, step))))

                if LINEAR_CAP is not None:
                    # Axis ends exactly at the top of the catch-all band -- no
                    # extra proportional headroom beyond that fixed-size zone.
                    ylim_hi = _shift(CATCHALL_BOUNDARY + BAND_HEIGHT)
                    ax_box.set_ylim(ylim_lo, ylim_hi)

                    # Catch-all band for the adaptive cap -- same idea as
                    # the 0 band below, just at the top: values pushed up
                    # here are display-only (jittered around CATCHALL_Y,
                    # not their true magnitude), so it needs to read as
                    # clearly distinct from the real 1-LINEAR_CAP scale
                    # below it. Boundary/label positions are pushed through
                    # the same ZERO_GAP shift as everything else so they
                    # still land in the right place on the offset axis.
                    catchall_boundary = _shift(CATCHALL_BOUNDARY)
                    catchall_y = _shift(CATCHALL_Y)
                    ax_box.axhspan(catchall_boundary, ylim_hi, color='#D5D5D5', zorder=0.5)
                    ax_box.axhline(catchall_boundary, color='#888888', lw=1.0, zorder=0.6)
                    ax_box.text(
                        0.01, catchall_y, f'>{LINEAR_CAP}', color='#555555', fontsize=10,
                        ha='left', va='center', zorder=4, transform=ax_box.get_yaxis_transform(),
                    )

                    # y-ticks at true values (0, 1, and round-number
                    # multiples of 5/10 up to the cap), placed at their
                    # ZERO_GAP-shifted plotted positions so the axis labels
                    # still read correctly.
                    true_ticks = _round_number_ticks(LINEAR_CAP)
                    ax_box.set_yticks([_shift(t) for t in true_ticks])
                    ax_box.set_yticklabels([str(t) for t in true_ticks])
                else:
                    # Axis ends exactly at the max value -- just enough
                    # extra room (half a marker) so the topmost point isn't
                    # clipped by the axes frame.
                    max_val = swarm_df['gene_count'].max()
                    ax_box.set_ylim(ylim_lo, _shift(max_val) + 0.5)

                    # Same true-value tick remapping as above, sized to the
                    # uncapped data range.
                    true_ticks = _round_number_ticks(max_val)
                    ax_box.set_yticks([_shift(t) for t in true_ticks])
                    ax_box.set_yticklabels([str(t) for t in true_ticks])


        ax_box.set_xlabel('')
        ax_box.set_ylabel("Gene count")
        ax_box.set_title(args.title + ("  (log scale)" if log_scale else ""))

        # Legend: one entry per sample_type, in the same color-assignment
        # order, labeled with that sample_type's total cohort size.
        handles = [
            plt.Line2D([0], [0], marker='o', linestyle='', color=color_map[st],
                       label=f"{st} (n={st_total_counts[st]})", markersize=7)
            for st in sample_types_unique
        ]
        ax_box.legend(handles=handles, title="sample_type", bbox_to_anchor=(1.02, 1), loc="upper left")

        # Combination matrix (standard UpSet style): one row per category, one
        # column per combo (at that combo's group center, spanning all its
        # sample_type sub-columns above), filled dot if that category is in
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