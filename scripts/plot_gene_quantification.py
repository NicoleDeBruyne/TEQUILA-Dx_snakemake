"""
plot_gene_quantification.py

Cohort-level per-gene expression box plots for one quantification method (read count, max
coverage, splice-site assignment or AMALGAM). One PDF per gene (US Letter, landscape) in --outdir:
CPTM on the top row and MOTR on the bottom row, one panel per sample type, each with its own
y-axis. Every sample is a dot in its sample type's colour under a plain black box; the highest and
lowest sample of each panel are larger, in a darker shade, and labelled.

Inputs are the per-sample-type matrices written by quantify_gene_expression.py /
quantify_gene_by_assignment.py / normalize_amalgam_matrix.py (gene x sample, first column 'gene').
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

from palette import darker

rcParams["pdf.fonttype"] = 42
PAGE_W, PAGE_H = 11.0, 8.5
_BLACK = dict(color="black", linewidth=0.9)


def parse_args():
    p = argparse.ArgumentParser(description="Per-gene CPTM / MOTR box plots by sample type, one PDF per gene.")
    p.add_argument("--cptm", nargs="+", required=True, metavar="SAMPLE_TYPE:PATH",
                   help="Each sample type's CPTM matrix (gene x sample).")
    p.add_argument("--motr", nargs="+", required=True, metavar="SAMPLE_TYPE:PATH",
                   help="Each sample type's MOTR matrix (gene x sample), same sample types.")
    p.add_argument("--colors", nargs="*", default=[], metavar="SAMPLE_TYPE:COLOR",
                   help="Each sample type's colour.")
    p.add_argument("--metric-label", default="", help="Quantification method, for titles (e.g. 'read count').")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def _pairs(items):
    out = {}
    for it in items:
        k, _, v = it.partition(":")
        if v:
            out[k] = v.strip("'\"")
    return out


def _read(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def _panel(ax, vals, color, title, ylabel):
    vals = vals.dropna()
    if vals.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes, color="#777")
        ax.set_xticks([]); ax.set_yticks([])
    else:
        rng = np.random.RandomState(0)
        jitter = (rng.rand(len(vals)) - 0.5) * 0.3
        ax.scatter(1 + jitter, vals.to_numpy(), color=color, s=18, alpha=0.8, linewidths=0, zorder=1)
        ax.boxplot([vals.to_numpy()], widths=0.5, showfliers=False, zorder=2,
                   boxprops=_BLACK, whiskerprops=_BLACK, capprops=_BLACK, medianprops=_BLACK)
        dark = darker(color)
        idx = list(vals.index)
        for s in dict.fromkeys([vals.idxmax(), vals.idxmin()]):
            i = idx.index(s)
            ax.scatter([1 + jitter[i]], [vals[s]], color=dark, s=40, zorder=3, linewidths=0)
            ax.annotate(s, (1 + jitter[i], vals[s]), textcoords="offset points", xytext=(7, 0),
                        fontsize=7, color=dark, va="center")
        ax.set_xticks([])
        ax.set_xlim(0.5, 1.9)
    ax.set_title(title, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main():
    args = parse_args()
    cptm = {t: _read(p) for t, p in _pairs(args.cptm).items() if os.path.isfile(p)}
    motr = {t: _read(p) for t, p in _pairs(args.motr).items() if os.path.isfile(p)}
    colors = _pairs(args.colors)
    types = [t for t in _pairs(args.cptm) if t in cptm or t in motr]
    genes = sorted(set().union(*[set(df.index) for df in list(cptm.values()) + list(motr.values())]))
    os.makedirs(args.outdir, exist_ok=True)
    label = f" ({args.metric_label})" if args.metric_label else ""

    for gene in genes:
        fig, axes = plt.subplots(2, len(types), figsize=(PAGE_W, PAGE_H), squeeze=False)
        fig.suptitle(f"{gene}{label}", fontsize=13, fontweight="bold", x=0.07, ha="left")
        for j, t in enumerate(types):
            col = colors.get(t, "#7f7f7f")
            for i, (name, mats) in enumerate((("CPTM", cptm), ("MOTR", motr))):
                df = mats.get(t)
                vals = df.loc[gene] if df is not None and gene in df.index else pd.Series(dtype=float)
                n = int(vals.notna().sum())
                _panel(axes[i, j], vals, col, f"{t} (n={n})" if i == 0 else "",
                       name if j == 0 else None)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(os.path.join(args.outdir, f"{gene}.pdf"))
        plt.close(fig)
    print(f"Saved {len(genes)} per-gene PDF(s) to {args.outdir}")


if __name__ == "__main__":
    main()
