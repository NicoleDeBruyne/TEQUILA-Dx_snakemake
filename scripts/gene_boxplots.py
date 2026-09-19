#!/usr/bin/env python3
"""
scripts/gene_boxplots.py
Shared per-gene boxplot renderer, used by every quantification method
(scripts/quantify_gene_expression.py, scripts/quantify_gene_by_assignment.py,
scripts/normalize_amalgam_matrix.py) to write one <gene>_cptm.pdf and one
<gene>_motr.pdf per targeted-panel gene -- same plot convention (boxplot +
jittered per-sample points, highest/lowest sample labeled) regardless of
which normalization's matrix is being plotted, so call this once per
normalization with that normalization's own value_df and filename_suffix.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def make_gene_boxplots(value_df, outdir, y_label, filename_suffix):
    """One page per gene: a boxplot of every sample's value (CPTM, MOTR, or
    whichever normalization value_df holds) for that gene, with individual
    sample points overlaid and the highest- and lowest-value sample
    labeled by name. filename_suffix ("_cptm" or "_motr") picks which of
    the two per-gene PDFs (<gene>_cptm.pdf / <gene>_motr.pdf) this call
    writes -- call once per normalization."""
    os.makedirs(outdir, exist_ok=True)
    for gene in value_df.index:
        vals = value_df.loc[gene].dropna()
        out_pdf = os.path.join(outdir, gene + filename_suffix + ".pdf")
        if vals.empty:
            continue

        fig, ax = plt.subplots(figsize=(3.6, 4.2))
        _black = dict(color="black")
        ax.boxplot([vals.to_numpy()], showfliers=False, widths=0.5,
                   boxprops=_black, whiskerprops=_black, capprops=_black, medianprops=_black)

        jitter = (np.random.RandomState(0).rand(len(vals)) - 0.5) * 0.15
        ax.scatter(1 + jitter, vals.to_numpy(), color="#2c7fb8", zorder=3, s=18)

        max_sample = vals.idxmax()
        min_sample = vals.idxmin()
        ax.scatter([1 + jitter[list(vals.index).index(max_sample)]], [vals[max_sample]],
                   color="#c0392b", zorder=4, s=30)
        ax.annotate(max_sample, (1 + jitter[list(vals.index).index(max_sample)], vals[max_sample]),
                    textcoords="offset points", xytext=(6, 0), fontsize=7, color="#c0392b", va="center")
        if min_sample != max_sample:
            ax.scatter([1 + jitter[list(vals.index).index(min_sample)]], [vals[min_sample]],
                       color="#c0392b", zorder=4, s=30)
            ax.annotate(min_sample, (1 + jitter[list(vals.index).index(min_sample)], vals[min_sample]),
                        textcoords="offset points", xytext=(6, 0), fontsize=7, color="#c0392b", va="center")

        ax.set_xticks([])
        ax.set_ylabel(y_label)
        ax.set_title(gene, fontsize=10)
        fig.tight_layout()
        fig.savefig(out_pdf)
        plt.close(fig)
