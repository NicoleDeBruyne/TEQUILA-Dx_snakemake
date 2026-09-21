
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
    n_genes = len(flr_df.index)
    box_h = 1.6
    bar_h = 1.2
    if flr_df.empty or n_genes == 0:
        return 1.5
    row_height = 0.22 if n_genes <= 60 else max(0.05, 0.22 * (60 / n_genes))
    heat_h = min(max(4, row_height * n_genes), 36)
    return heat_h + _ROW_GAP + box_h + _ROW_GAP + bar_h


_ROW_GAP = 0.15
_COL_GAP = 0.15
_OUTER_GAP_MIN = 0.3
_BAR_W = 1.6
_CBAR_W = 0.25
_MARGIN = 0.15
_TOP_MARGIN = 0.75


def _draw_heatmap_panel(fig, x0, y_top, flr_df, gene_median_reads, sample_total_reads,
                         sample_colors, panel_title, heat_w, full_w, full_h):
    n_genes = len(flr_df.index)
    box_h = 1.6
    bar_h = 1.2

    def rect(x, y, w, h):
        return [x / full_w, y / full_h, w / full_w, h / full_h]

    if flr_df.empty:
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

    show_gene_labels = (heat_h / max(n_genes, 1)) * 72 >= 5
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

    y = np.arange(n_genes) + 0.5
    ax_gbar.barh(y, gmed.clip(lower=1).to_numpy(), height=0.8, color="#555555")
    ax_gbar.set_xscale("log")
    ax_gbar.set_xlabel("median\nreads", fontsize=7)
    ax_gbar.tick_params(axis="y", left=False, labelleft=False)
    ax_gbar.set_ylim(ax_heat.get_ylim())
    ax_gbar.set_xlim(left=1)

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

    ax_sbar.bar(x, stot.clip(lower=1).to_numpy(), width=0.8, color=scolors)
    ax_sbar.set_yscale("log")
    ax_sbar.set_ylim(bottom=1)
    ax_sbar.set_ylabel("total\nreads", fontsize=7)
    ax_sbar.tick_params(axis="x", labelbottom=False, bottom=False)


def make_combined_heatmap(flr_high, flr_low, gene_median_reads, sample_total_reads,
                           sample_type_map, sample_colors, out_pdf, title_high, title_low,
                           suptitle, n_samples):
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

    alias_map = parse_alias_map(args.alias_map)
    alias_matrix = matrix.copy()
    alias_matrix.columns = resolve_all(alias_matrix.columns, alias_map)
    out_matrix_alias = args.outprefix + "_matrix_alias.tsv"
    alias_matrix.to_csv(out_matrix_alias, sep="\t")
    print("Saved alias-labeled avgFLR matrix: " + out_matrix_alias)

    out_counts = args.outprefix + "_read_counts_matrix.tsv"
    count_matrix.to_csv(out_counts, sep="\t")
    print("Saved read-count matrix: " + out_counts)

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
