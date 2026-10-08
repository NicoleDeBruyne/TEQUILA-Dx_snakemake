"""
Merge every sample's {sample}_avg_coverage.tsv (from get_avg_coverage_sample.py) into cohort-level
gene x sample matrices and heatmaps -- one set each for average_coverage, median_coverage and
breadth_coverage.

Outputs (per metric, with <metric> = avg_coverage / median_coverage / breadth_coverage):
    <outprefix>_<metric>_matrix.tsv         -- genes (rows) x samples (columns)
    <outprefix>_<metric>_matrix_alias.tsv   -- same, sample columns relabeled by --alias-map
    <outprefix>_<metric>_heatmap.pdf        -- heatmap + per-sample on-target read count bar chart
    <outprefix>_<metric>_heatmap_alias.pdf  -- same, sample labels relabeled by --alias-map

Heatmap colors run red -> white -> blue. For the depth metrics (average/median): darkest red at 0,
white at --center (default 50), darkest blue at the matrix maximum; if no value exceeds --center, the
scale still tops out just above it so white always means the same depth. For breadth (fraction of
canonical exonic bases at depth >= the per-sample --min-depth, 50 by default): a fixed 0-1 scale,
darkest red at 0, white at 0.5, darkest blue at 1. Missing values (no transcript in the GTF) are
drawn grey.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches
import seaborn as sns
from matplotlib import rcParams

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
rcParams['pdf.fonttype'] = 42


METRICS = [
    # (output name, per-sample TSV column, label, kind)
    ("avg_coverage", "average_coverage", "average coverage", "depth"),
    ("median_coverage", "median_coverage", "median coverage", "depth"),
    ("breadth_coverage", "breadth_coverage", "breadth of coverage", "fraction"),
]

_ROW_GAP = 0.15
_COL_GAP = 0.15
_CBAR_W = 0.25
_MARGIN = 0.15
_TOP_MARGIN = 0.75
_BAR_H = 1.2
_LEGEND_W = 1.3


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge per-sample canonical-transcript coverage TSVs into cohort matrices + heatmaps")
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Every sample's {sample}_avg_coverage.tsv, in the same order as --samples")
    parser.add_argument("--on-target-files", nargs="+", required=True,
        help="Every sample's {sample}_on_target.tsv (from _6A), in the same order as --samples; "
             "its 'target' column is plotted as the per-sample on-target read count bar chart")
    parser.add_argument("--samples", nargs="+", required=True, help="Sample IDs, one per --infiles entry")
    parser.add_argument("--sample-types", nargs="*", default=[],
        help="Optional per-sample sample_type label (same order as --samples), for the bar chart legend")
    parser.add_argument("--colors", nargs="*", default=[],
        help="Optional per-sample pre-resolved hex color (same order as --samples), for the bar chart")
    parser.add_argument("--outprefix", required=True, help="Prefix for all output files.")
    parser.add_argument("--title", default="Canonical transcript coverage", help="Plot title prefix.")
    parser.add_argument("--center", type=float, default=50.0,
        help="Depth drawn as white in the heatmap (red below, blue above). Default: 50")
    add_alias_map_arg(parser)
    return parser.parse_args()


def _check_len(name, values, n):
    if values and len(values) != n:
        raise ValueError(name + " must have one entry per --samples entry, or be omitted entirely")


def make_heatmap(matrix, on_target, sample_colors, sample_types, sample_labels, center,
                 metric_label, title, out_pdf, kind="depth"):
    """matrix: genes x samples (already ordered). sample_labels: x tick labels in the same column order."""
    n_genes, n_samples = matrix.shape

    col_width = 0.28 if n_samples <= 50 else max(0.05, 0.28 * (50 / n_samples))
    heat_w = min(max(5, col_width * n_samples + 1), 20)
    row_height = 0.22 if n_genes <= 60 else max(0.05, 0.22 * (60 / n_genes))
    heat_h = min(max(4, row_height * n_genes), 36)

    show_sample_labels = (heat_w / max(n_samples, 1)) * 72 >= 5
    sample_fontsize = min(8, max(4, (heat_w * 72) / max(n_samples, 1) - 1))
    label_h = (0.12 * max((len(str(s)) for s in sample_labels), default=0) * sample_fontsize / 8
               if show_sample_labels else 0)

    full_w = _MARGIN + heat_w + _COL_GAP + max(_CBAR_W + 0.6, _LEGEND_W) + _MARGIN
    full_h = _TOP_MARGIN + heat_h + _ROW_GAP + _BAR_H + label_h + _MARGIN

    def rect(x, y, w, h):
        return [x / full_w, y / full_h, w / full_w, h / full_h]

    fig = plt.figure(figsize=(full_w, full_h))
    y_heat_bottom = full_h - _TOP_MARGIN - heat_h
    y_bar_bottom = y_heat_bottom - _ROW_GAP - _BAR_H

    if matrix.empty:
        ax = fig.add_axes(rect(_MARGIN, y_heat_bottom, heat_w, heat_h))
        ax.axis("off")
        ax.text(0.5, 0.5, "No genes / samples to plot.", ha="center", va="center")
        fig.suptitle(title)
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        return

    ax_heat = fig.add_axes(rect(_MARGIN, y_heat_bottom, heat_w, heat_h))
    ax_cbar = fig.add_axes(rect(_MARGIN + heat_w + _COL_GAP, y_heat_bottom, _CBAR_W, heat_h))
    ax_sbar = fig.add_axes(rect(_MARGIN, y_bar_bottom, heat_w, _BAR_H), sharex=ax_heat)

    if kind == "fraction":
        center, vmax = 0.5, 1.0
    else:
        data_max = np.nanmax(matrix.to_numpy()) if np.isfinite(matrix.to_numpy()).any() else center
        vmax = max(float(data_max), center + 1.0)
    norm = matplotlib.colors.TwoSlopeNorm(vmin=0.0, vcenter=center, vmax=vmax)
    cmap = matplotlib.colormaps["RdBu"].copy()
    cmap.set_bad("#BBBBBB")

    show_gene_labels = (heat_h / max(n_genes, 1)) * 72 >= 5
    gene_fontsize = min(8, max(4, (heat_h * 72) / max(n_genes, 1) - 1))

    sns.heatmap(matrix, ax=ax_heat, cmap=cmap, norm=norm,
                cbar=True, cbar_ax=ax_cbar, cbar_kws={"label": metric_label + (" (fraction of bases)" if kind == "fraction" else " (reads/base)")},
                linewidths=0, yticklabels=show_gene_labels, xticklabels=False)
    cbar = ax_heat.collections[0].colorbar
    ticks = sorted({0.0, center / 2, center, center + (vmax - center) / 2, vmax})
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.0%}" if kind == "fraction" else f"{t:,.0f}" for t in ticks])
    cbar.ax.tick_params(labelsize=7)
    if show_gene_labels:
        ax_heat.tick_params(axis="y", labelsize=gene_fontsize)
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("gene")
    ax_heat.tick_params(axis="x", labelbottom=False, bottom=False)

    x = np.arange(n_samples) + 0.5
    scolors = [sample_colors.get(s, "#555555") for s in matrix.columns]
    ax_sbar.bar(x, on_target.reindex(matrix.columns).fillna(0).clip(lower=1).to_numpy(),
                width=0.8, color=scolors)
    ax_sbar.set_yscale("log")
    ax_sbar.set_ylim(bottom=1)
    ax_sbar.set_ylabel("on-target\nreads", fontsize=7)
    ax_sbar.tick_params(axis="y", labelsize=7)
    if show_sample_labels:
        ax_sbar.set_xticks(x)
        ax_sbar.set_xticklabels(sample_labels, rotation=90, fontsize=sample_fontsize)
    else:
        ax_sbar.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_sbar.set_xlim(0, n_samples)

    if sample_types:
        seen = {}
        for s in matrix.columns:
            if s in sample_types:
                seen.setdefault(sample_types[s], sample_colors.get(s, "#555555"))
        handles = [matplotlib.patches.Patch(color=c, label=t) for t, c in seen.items()]
        legend_x = (_MARGIN + heat_w + _COL_GAP) / full_w
        legend_y = (y_bar_bottom + _BAR_H) / full_h
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(legend_x, legend_y),
                   title="sample_type", fontsize=7, title_fontsize=7, frameon=False)

    fig.suptitle(title, y=1 - (_TOP_MARGIN / full_h) * 0.35)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    n = len(args.samples)
    if len(args.infiles) != n or len(args.on_target_files) != n:
        raise ValueError("--infiles, --on-target-files and --samples must all have the same length")
    _check_len("--sample-types", args.sample_types, n)
    _check_len("--colors", args.colors, n)
    samples = list(args.samples)
    sample_types = dict(zip(samples, args.sample_types)) if args.sample_types else {}
    sample_colors = dict(zip(samples, args.colors)) if args.colors else {}
    alias_map = parse_alias_map(args.alias_map)

    per_sample = {s: pd.read_csv(f, sep="\t") for s, f in zip(samples, args.infiles)}
    genes = list(pd.concat([d["gene"] for d in per_sample.values()]).drop_duplicates())

    on_target = pd.Series(
        {s: float(pd.read_csv(f, sep="\t")["target"].iloc[0]) for s, f in zip(samples, args.on_target_files)}
    )

    for name, column, label, kind in METRICS:
        missing = [s for s, d in per_sample.items() if column not in d.columns]
        if missing:
            raise ValueError(f"Column '{column}' missing from the coverage file of: {missing}. "
                             "Rerun _6D_get_avg_coverage to regenerate them.")
        matrix = pd.DataFrame(
            {s: d.set_index("gene")[column] for s, d in per_sample.items()}
        ).reindex(index=genes, columns=samples)
        matrix.index.name = "gene"

        out_matrix = args.outprefix + "_" + name + "_matrix.tsv"
        matrix.to_csv(out_matrix, sep="\t", float_format="%.4f")
        alias_matrix = matrix.copy()
        alias_matrix.columns = resolve_all(alias_matrix.columns, alias_map)
        alias_matrix.to_csv(args.outprefix + "_" + name + "_matrix_alias.tsv", sep="\t", float_format="%.4f")
        print("Saved " + label + " matrix: " + out_matrix)

        gene_order = matrix.mean(axis=1).sort_values(ascending=False, na_position="last").index
        sample_order = matrix.mean(axis=0).sort_values(ascending=False, na_position="last").index
        ordered = matrix.loc[gene_order, sample_order]

        white = "50%" if kind == "fraction" else f"{args.center:g}"
        title = (args.title + " -- " + label + " (" + str(len(genes)) + " gene(s), "
                 + str(n) + " sample(s); white = " + white + ")")
        for suffix, labels in (("", list(sample_order)),
                               ("_alias", resolve_all(sample_order, alias_map))):
            out_pdf = args.outprefix + "_" + name + "_heatmap" + suffix + ".pdf"
            make_heatmap(ordered, on_target, sample_colors, sample_types, labels, args.center,
                         label, title, out_pdf, kind=kind)
            print("Saved " + label + " heatmap: " + out_pdf)


if __name__ == "__main__":
    main()
