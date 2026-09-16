#!/usr/bin/env python3
"""
scripts/get_full_length_ratio.py
Cohort-level merge step: combines every sample's own
{sample}_full_length_ratio.tsv (per-gene avgFLR + read_count, written
per-sample by scripts/get_full_length_ratio_sample.py via
rules/6_sample_qc.smk's _6C) into the cohort's avgFLR / read-count matrices
and heatmap. No BAM or GTF access here -- adding/removing a sample from a
cohort only reruns this cheap merge, not the per-sample BAM scan + GTF
parse. Invoked by rules/9_merge_results.smk's _9J.

Genes are split into two heatmaps by their median (across samples) read
count: genes with too few supporting reads have noisy/meaningless avgFLR,
so grouping them separately (rather than interleaving blank-looking low-
confidence cells throughout one big heatmap) keeps the well-supported genes
readable and makes clear which genes are low-confidence and why. Both live
as stacked subplots in one figure/PDF (high-read-count on top, low-read-
count below) so they're easy to compare side by side. Each panel carries
its own per-gene median-read-count bar (log-scaled, right), a per-sample
avgFLR boxplot (this panel's own genes, one jittered point per gene per
sample), and a per-sample total-read-count bar (log-scaled, bottom, same
value in both panels since it's summed across every gene, not just the
ones shown in that panel). The boxplot points and the total-reads bars are
both colored by sample_type (using the same colors as
_9I2_validate_sample_types -- resolved via the Snakefile's own
sample_type_color() and passed through --sample-types/--colors, not
recomputed here), with a shared legend on the figure -- no separate color
strip. Column width per sample shrinks (and the overall figure width is
capped) as sample count grows, so large cohorts (e.g. ~200 samples) stay a
print-friendly width instead of many feet wide. No sample-name labels are
drawn on the heatmap axis itself -- see git history / conversation notes.

Outputs:
  --outprefix + "_matrix.tsv"              -- genes x samples avgFLR matrix
                                              (every gene, unsplit).
  --outprefix + "_read_counts_matrix.tsv"  -- genes x samples read-count
                                              matrix backing the avgFLR
                                              values and the bar charts.
  --outprefix + "_heatmap.pdf"             -- one figure, two stacked
                                              subplots: avgFLR heatmap +
                                              per-sample avgFLR boxplot +
                                              read-count bar, for genes
                                              with median read count >=
                                              --min-reads (default 100) on
                                              top, and < --min-reads below.
"""

import argparse

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches
import seaborn as sns
from matplotlib import rcParams

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
rcParams['pdf.fonttype'] = 42


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge per-sample full-length-ratio TSVs into a cohort's avgFLR matrices + heatmap")
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Every sample's {sample}_full_length_ratio.tsv, in the desired plot order")
    parser.add_argument("--sample-types", nargs="*", default=[],
        help="Optional per-sample sample_type label (same order as --infiles), for boxplot/bar coloring")
    parser.add_argument("--colors", nargs="*", default=[],
        help="Optional per-sample pre-resolved hex color (same order/length as --sample-types -- same "
             "colors _9I2_validate_sample_types uses, resolved by the calling rule via the Snakefile's "
             "sample_type_color(), not recomputed here)")
    parser.add_argument("--outprefix", required=True, help="Prefix for all output files.")
    parser.add_argument("--title", default="Full-length transcript coverage ratio (FLR)", help="Plot title prefix.")
    parser.add_argument("--min-reads", type=int, default=100,
        help="Genes whose median (across samples) read count is below this are split into their own "
             "low-confidence heatmap instead of the main one. Default: 100")
    add_alias_map_arg(parser)
    return parser.parse_args()


def _panel_height(flr_df):
    """Height (inches) _draw_heatmap_panel will need for this gene set --
    computed upfront (pure function of gene count) so make_combined_heatmap
    can lay out both panels' absolute inch positions before any subplot
    exists. Includes the fixed row gaps below, so this is the panel's true
    total footprint, not just its content rows summed."""
    n_genes = len(flr_df.index)
    box_h = 1.6   # per-sample avgFLR boxplot row
    bar_h = 1.2   # per-sample total-read-count bar row
    if flr_df.empty or n_genes == 0:
        return 1.5
    row_height = 0.22 if n_genes <= 60 else max(0.05, 0.22 * (60 / n_genes))
    heat_h = min(max(4, row_height * n_genes), 36)
    return heat_h + _ROW_GAP + box_h + _ROW_GAP + bar_h


# Fixed, absolute (inches) gaps -- deliberately NOT expressed as a GridSpec
# hspace/wspace fraction. This panel mixes one huge row (the heatmap, up to
# 36in tall) with two small fixed-height rows (the boxplot and reads bar,
# 1.6in/1.2in): a fractional hspace is a fraction of the *average* row
# height across the whole grid, so with one giant row in the mix that
# average is dominated by the heatmap, and the same fraction blows up into
# a huge absolute gap between the two small rows -- which is what produced
# the excessive whitespace this layout replaces. Positioning every axes by
# its own absolute inch coordinates keeps every gap the same fixed size
# regardless of how tall the heatmap row ends up being.
_ROW_GAP = 0.15    # between heat/box and box/bar, inches
_COL_GAP = 0.15    # between heat/bar and bar/cbar, inches
_OUTER_GAP_MIN = 0.3   # between the two stacked (high/low) panels, inches
_BAR_W = 1.6       # right (gene median reads) bar width, inches
_CBAR_W = 0.25     # colorbar width, inches
_MARGIN = 0.15     # figure margin on each side, inches (bbox_inches="tight" trims any excess)
_TOP_MARGIN = 0.75 # figure top margin, inches -- reserves room for the suptitle + legend, which both
                   # sit in figure-level space above the top panel's own axes


def _draw_heatmap_panel(fig, x0, y_top, flr_df, gene_median_reads, sample_total_reads,
                         sample_colors, panel_title, heat_w, full_w, full_h):
    """Draws one avgFLR heatmap, a per-sample avgFLR boxplot (this panel's
    own genes, points colored by sample_type), and a per-sample total-
    read-count bar (also colored by sample_type), placed via absolute-inch
    axes positions rooted at (x0, y_top) -- see the gap constants above for
    why this isn't done via GridSpec height/width ratios. `x0`/`y_top` are
    in inches from the figure's bottom-left corner; `full_w`/`full_h` are
    the whole figure's size in inches, used to convert to the [0, 1]
    figure-fraction coordinates fig.add_axes() expects. `sample_total_reads`
    and `sample_colors` are expected to already cover every sample (not
    just the ones in `flr_df`'s columns), so they read the same across
    every panel."""
    n_genes = len(flr_df.index)
    box_h = 1.6
    bar_h = 1.2

    def rect(x, y, w, h):
        return [x / full_w, y / full_h, w / full_w, h / full_h]

    if flr_df.empty:
        # Still draw something -- this panel is part of a declared
        # Snakemake output, so the figure must exist even when a run's gene
        # panel happens to have nothing in this read-count bucket.
        ax = fig.add_axes(rect(x0, y_top - 1.5, heat_w + _COL_GAP + _BAR_W, 1.5))
        ax.axis("off")
        ax.text(0.5, 0.5, "No genes in this group.", ha="center", va="center")
        ax.set_title(panel_title)
        return

    gene_order = flr_df.mean(axis=1).sort_values(ascending=False).index
    sample_order = flr_df.mean(axis=0).sort_values(ascending=False).index
    ordered = flr_df.loc[gene_order, sample_order]
    n_samples = len(sample_order)
    gmed = gene_median_reads.reindex(gene_order)
    stot = sample_total_reads.reindex(sample_order)
    scolors = [sample_colors.get(s, "#555555") for s in sample_order]

    row_height = 0.22 if n_genes <= 60 else max(0.05, 0.22 * (60 / n_genes))
    heat_h = min(max(4, row_height * n_genes), 36)

    # Only individually label genes when there's room to read them; past
    # that, still-visible sort order (best/worst avgFLR) carries the
    # information without needing every gene name printed.
    show_gene_labels = (heat_h / max(n_genes, 1)) * 72 >= 5  # ~5pt/row floor
    gene_fontsize = min(8, max(4, (heat_h * 72) / max(n_genes, 1) - 1))

    y_heat_bottom = y_top - heat_h
    y_box_bottom = y_heat_bottom - _ROW_GAP - box_h
    y_bar_bottom = y_box_bottom - _ROW_GAP - bar_h
    x_bar = x0 + heat_w + _COL_GAP
    x_cbar = x_bar + _BAR_W + _COL_GAP

    ax_heat = fig.add_axes(rect(x0, y_heat_bottom, heat_w, heat_h))
    ax_gbar = fig.add_axes(rect(x_bar, y_heat_bottom, _BAR_W, heat_h), sharey=ax_heat)
    ax_cbar = fig.add_axes(rect(x_cbar, y_heat_bottom, _CBAR_W, heat_h))
    ax_box  = fig.add_axes(rect(x0, y_box_bottom, heat_w, box_h), sharex=ax_heat)
    ax_sbar = fig.add_axes(rect(x0, y_bar_bottom, heat_w, bar_h), sharex=ax_heat)

    sns.heatmap(ordered, ax=ax_heat, cmap="viridis", vmin=0, vmax=1,
                cbar=True, cbar_ax=ax_cbar, cbar_kws={"label": "avgFLR"},
                linewidths=0, yticklabels=show_gene_labels)
    if show_gene_labels:
        ax_heat.tick_params(axis="y", labelsize=gene_fontsize)
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("gene")
    ax_heat.set_title(panel_title, fontsize=10)
    ax_heat.tick_params(axis="x", labelbottom=False, bottom=False)

    # Per-gene median read count, right, log-scaled. Clipped to >=1 so
    # zero-read genes (common in the low-read-count panel) still show up as
    # a visible sliver rather than vanishing on a log axis.
    y = np.arange(n_genes) + 0.5
    ax_gbar.barh(y, gmed.clip(lower=1).to_numpy(), height=0.8, color="#555555")
    ax_gbar.set_xscale("log")
    ax_gbar.set_xlabel("median\nreads", fontsize=7)
    ax_gbar.tick_params(axis="y", left=False, labelleft=False)
    ax_gbar.set_ylim(ax_heat.get_ylim())
    ax_gbar.set_xlim(left=1)

    # Per-sample avgFLR distribution across this panel's own genes -- one
    # jittered point per gene colored by that sample's sample_type, with
    # the box outline drawn ON TOP of those points (zorder=3 vs the
    # scatter's zorder=2) so the box/whiskers/median stay visible even
    # where points are dense, rather than being buried under them. Column
    # positions (x = i+0.5) match the heatmap above it exactly, same
    # convention as every other row in this panel.
    x = np.arange(n_samples) + 0.5
    box_data = [ordered[s].dropna().to_numpy() for s in sample_order]
    rng = np.random.RandomState(0)
    for xi, vals, c in zip(x, box_data, scolors):
        if len(vals) == 0:
            continue
        jitter = (rng.rand(len(vals)) - 0.5) * 0.35
        ax_box.scatter(np.full(len(vals), xi) + jitter, vals, color=c, s=6, zorder=2, edgecolors="none")
    _black = dict(color="black", linewidth=0.7, zorder=3)
    bp = ax_box.boxplot(box_data, positions=x, widths=0.6, showfliers=False,
                         boxprops=_black, whiskerprops=_black, capprops=_black, medianprops=_black)
    ax_box.set_ylim(-0.02, 1.02)
    ax_box.set_ylabel("avgFLR", fontsize=7)
    ax_box.tick_params(axis="x", labelbottom=False, bottom=False)

    # Per-sample total read count (summed across every gene, not just this
    # panel's subset), bottom, log-scaled, colored the same as the boxplot
    # points above (i.e. by sample_type).
    ax_sbar.bar(x, stot.clip(lower=1).to_numpy(), width=0.8, color=scolors)
    ax_sbar.set_yscale("log")
    ax_sbar.set_ylim(bottom=1)
    ax_sbar.set_ylabel("total\nreads", fontsize=7)
    ax_sbar.tick_params(axis="x", labelbottom=False, bottom=False)


def make_combined_heatmap(flr_high, flr_low, gene_median_reads, sample_total_reads,
                           sample_type_map, sample_colors, out_pdf, title_high, title_low,
                           suptitle, n_samples):
    """One figure, two stacked panels (high-read-count panel on top,
    low-read-count panel below), each built by _draw_heatmap_panel via
    absolute-inch axes positions -- see the gap constants above
    _draw_heatmap_panel for why this isn't GridSpec height/width ratios.

    Column width shrinks as sample count grows (same idea as the row-height
    shrink for genes) and the total width is capped, so a ~200-sample
    cohort gets a print-friendly figure instead of a many-foot-wide PDF."""
    col_width = 0.28 if n_samples <= 50 else max(0.05, 0.28 * (50 / n_samples))
    heat_w = min(max(5, col_width * n_samples + 1), 20)
    full_w = _MARGIN + heat_w + _COL_GAP + _BAR_W + _COL_GAP + _CBAR_W + _MARGIN

    h_high = _panel_height(flr_high)
    h_low = _panel_height(flr_low)

    outer_gap = _OUTER_GAP_MIN

    full_h = _TOP_MARGIN + h_high + outer_gap + h_low + _MARGIN

    fig = plt.figure(figsize=(full_w, full_h))

    y_top_high = full_h - _TOP_MARGIN
    y_top_low = y_top_high - h_high - outer_gap

    _draw_heatmap_panel(fig, _MARGIN, y_top_high, flr_high, gene_median_reads, sample_total_reads,
                         sample_colors, title_high, heat_w, full_w, full_h)
    _draw_heatmap_panel(fig, _MARGIN, y_top_low, flr_low, gene_median_reads, sample_total_reads,
                         sample_colors, title_low, heat_w, full_w, full_h)

    # One legend for the whole figure -- sample_type -> color is the same
    # in both panels (and matches rules/9_merge_results.smk's
    # _9I2_validate_sample_types, which resolves colors via the Snakefile's
    # own sample_type_color()). Placed in the reserved top margin, above
    # the top panel's own axes, rather than overlapping it.
    seen = {}
    for s, t in sample_type_map.items():
        seen.setdefault(t, sample_colors.get(s, "#555555"))
    legend_handles = [matplotlib.patches.Patch(color=c, label=t) for t, c in seen.items()]
    fig.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(1, 1), title="sample_type",
               fontsize=7, title_fontsize=7, frameon=False)

    fig.suptitle(suptitle, y=1 - (_TOP_MARGIN / full_h) * 0.35)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()

    per_sample_dfs = [pd.read_csv(f, sep="\t") for f in args.infiles]
    samples = [d['sample'].iloc[0] for d in per_sample_dfs if not d.empty]

    sample_types = list(args.sample_types) if args.sample_types else []
    colors = list(args.colors) if args.colors else []
    if sample_types and len(sample_types) != len(args.infiles):
        raise ValueError("--sample-types must have one entry per --infiles entry, or be omitted entirely")
    if colors and len(colors) != len(args.infiles):
        raise ValueError("--colors must have one entry per --infiles entry, or be omitted entirely")
    sample_to_type = dict(zip(samples, sample_types)) if sample_types else {}
    sample_to_color = dict(zip(samples, colors)) if colors else {}

    long_df = pd.concat(per_sample_dfs, ignore_index=True)
    genes = list(long_df['gene'].drop_duplicates())

    matrix = long_df.pivot(index='gene', columns='sample', values='avgFLR').reindex(index=genes, columns=samples)
    matrix.index.name = "gene"

    count_matrix = (
        long_df.pivot(index='gene', columns='sample', values='read_count')
        .reindex(index=genes, columns=samples).fillna(0).astype(int)
    )
    count_matrix.index.name = "gene"

    out_matrix = args.outprefix + "_matrix.tsv"
    matrix.to_csv(out_matrix, sep="\t")
    print("Saved avgFLR matrix: " + out_matrix)

    # Alias-labeled copy: always produced (mirrors the real-ID matrix
    # verbatim when --alias-map is empty), so the rule's declared output
    # exists regardless of whether this bed panel actually has any aliases.
    alias_map = parse_alias_map(args.alias_map)
    alias_matrix = matrix.copy()
    alias_matrix.columns = resolve_all(alias_matrix.columns, alias_map)
    out_matrix_alias = args.outprefix + "_matrix_alias.tsv"
    alias_matrix.to_csv(out_matrix_alias, sep="\t")
    print("Saved alias-labeled avgFLR matrix: " + out_matrix_alias)

    out_counts = args.outprefix + "_read_counts_matrix.tsv"
    count_matrix.to_csv(out_counts, sep="\t")
    print("Saved read-count matrix: " + out_counts)

    # Per-gene median read count (across samples) decides which heatmap a
    # gene lands in; per-sample total read count (across every gene, not
    # just the genes in one heatmap) backs the bottom bar in both.
    gene_median_reads = count_matrix.median(axis=1)
    sample_total_reads = count_matrix.sum(axis=0)

    high_genes = gene_median_reads[gene_median_reads >= args.min_reads].index
    low_genes = gene_median_reads.index.difference(high_genes)

    out_heatmap = args.outprefix + "_heatmap.pdf"
    make_combined_heatmap(
        matrix.loc[high_genes], matrix.loc[low_genes], gene_median_reads, sample_total_reads,
        sample_to_type, sample_to_color,
        out_heatmap,
        title_high="median reads >= " + str(args.min_reads) + " (" + str(len(high_genes)) + " gene(s))",
        title_low="median reads < " + str(args.min_reads) + " (" + str(len(low_genes)) + " gene(s))",
        suptitle=args.title + " -- heatmap",
        n_samples=len(samples),
    )
    print("Saved avgFLR heatmap: " + out_heatmap +
          " (" + str(len(high_genes)) + " high-read-count, " + str(len(low_genes)) + " low-read-count gene(s))")


if __name__ == "__main__":
    main()
