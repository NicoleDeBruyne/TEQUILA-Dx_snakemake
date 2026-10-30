# 10R - merge the per-sample phasing vs WGS checks (3C) for a panel into one table + a quick plot

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

from palette import GREY, SELECT_ORDER, use_fonts
from sample_alias import add_alias_map_arg, parse_alias_map, resolve

use_fonts()
rcParams["pdf.fonttype"] = 42

MODE_ORDER = ["wgs_blinded", "wgs_aided"]
MODE_LABELS = {"wgs_blinded": "RNA only (WGS-blinded)", "wgs_aided": "WGS-aided"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine the per-sample _genes.tsv tables from _3C_check_phasing_vs_wgs.py and summarize "
                    "pairwise phase concordance with WGS per sample and phasing mode.")
    parser.add_argument("--gene-tables", nargs="+", required=True, help="Every WGS sample's *_genes.tsv.")
    parser.add_argument("--outprefix", required=True,
        help="Writes <outprefix>_genes.tsv, <outprefix>_summary.tsv and <outprefix>_summary.pdf.")
    add_alias_map_arg(parser)
    return parser.parse_args()


def summarize(genes):
    rows = []
    for (sample, mode), g in genes.groupby(["sample", "mode"], sort=False):
        pairs, adj = int(g["pairs"].sum()), int(g["adjacent_pairs"].sum())
        rows.append(dict(sample=sample, mode=mode,
                         phased_genes=len(g),
                         genes_with_pairs=int((g["pairs"] > 0).sum()),
                         wgs_phased_sites=int(g["wgs_phased_sites"].sum()),
                         assessed_sites=int(g["assessed_sites"].sum()),
                         pairs=pairs,
                         concordant_pairs=int(g["concordant_pairs"].sum()),
                         pair_concordance=g["concordant_pairs"].sum() / pairs if pairs else np.nan,
                         adjacent_pairs=adj,
                         adjacent_discordance=g["adjacent_discordant"].sum() / adj if adj else np.nan))
    return pd.DataFrame(rows)


def plot(summary, out):
    samples = list(dict.fromkeys(summary["label"]))
    modes = [m for m in MODE_ORDER if m in set(summary["mode"])] + \
            sorted(set(summary["mode"]) - set(MODE_ORDER))
    width = 0.8 / max(len(modes), 1)
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * len(samples) + 1.5), 3))
    for k, mode in enumerate(modes):
        sub = summary[summary["mode"] == mode].set_index("label").reindex(samples)
        xs = np.arange(len(samples)) + (k - (len(modes) - 1) / 2) * width
        ax.bar(xs, sub["pair_concordance"].fillna(0), width=width * 0.9, color=SELECT_ORDER[k % len(SELECT_ORDER)],
               label=MODE_LABELS.get(mode, mode))
        for x, n in zip(xs, sub["pairs"].fillna(0)):
            ax.text(x, 1.01, str(int(n)), ha="center", va="bottom", fontsize=5, color=GREY, rotation=90)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0, 1.15)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_ylabel("Site pairs phased as in WGS", fontsize=8)
    ax.set_title("Read phasing vs WGS phase (numbers = site pairs)", fontsize=8, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=6, frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    args = parse_args()
    tables = [pd.read_csv(p, sep="\t") for p in args.gene_tables if os.path.getsize(p) > 0]
    genes = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    if genes.empty:
        raise SystemExit("ERROR: all gene tables are empty, nothing to merge")

    alias_map = parse_alias_map(args.alias_map)
    genes.insert(1, "alias", [resolve(s, alias_map) for s in genes["sample"]])
    summary = summarize(genes)
    summary.insert(1, "alias", [resolve(s, alias_map) for s in summary["sample"]])
    summary["label"] = summary["alias"]

    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)), exist_ok=True)
    genes.to_csv(f"{args.outprefix}_genes.tsv", sep="\t", index=False)
    summary.drop(columns="label").to_csv(f"{args.outprefix}_summary.tsv", sep="\t", index=False)
    plot(summary, f"{args.outprefix}_summary.pdf")

    for mode, g in summary.groupby("mode"):
        pairs = g["pairs"].sum()
        rate = g["concordant_pairs"].sum() / pairs if pairs else float("nan")
        print(f"{mode}: {len(g)} samples, {int(pairs)} pairs, overall concordance {rate:.3f}")
    print(f"Wrote {args.outprefix}_summary.tsv, _genes.tsv and _summary.pdf")


if __name__ == "__main__":
    main()
