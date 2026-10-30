# 11B - merge soft-clip clusters across the group and flag the ones only one sample has

import argparse
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from palette import GREY, RED, use_fonts
from sample_alias import add_alias_map_arg, parse_alias_map, resolve


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge every sample's soft-clip clusters (rule _11A) for one group. A site is one gene, side "
                    "and exact breakpoint position. A sample has a site when its cluster there has a clipped "
                    "fraction of at least --present-fraction. A sample's cluster is an outlier when it has at least "
                    "--min-reads clipped reads, a clipped fraction of at least --min-fraction, a consensus identity "
                    "of at least --min-identity, and fewer than --max-sample-fraction of the group's samples other "
                    "than itself have the site.")
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--cluster-files", nargs="+", required=True, help="One _11A clusters TSV per --samples entry.")
    p.add_argument("--min-reads", type=int, default=10)
    p.add_argument("--min-fraction", type=float, default=0.2)
    p.add_argument("--present-fraction", type=float, default=0.05)
    p.add_argument("--min-identity", type=float, default=0.75)
    p.add_argument("--max-sample-fraction", type=float, default=0.05)
    p.add_argument("--title", default="")
    p.add_argument("--all-outfile", required=True)
    p.add_argument("--outliers-outfile", required=True)
    p.add_argument("--plot-outfile", required=True)
    add_alias_map_arg(p)
    a = p.parse_args()
    if len(a.samples) != len(a.cluster_files):
        p.error("--samples and --cluster-files must have the same length")
    return a


def main():
    a = parse_args()
    alias = parse_alias_map(a.alias_map)
    frames = []
    for s, path in zip(a.samples, a.cluster_files):
        if not os.path.isfile(path):
            print(f"WARNING: clusters file missing for {s}: {path}")
            continue
        df = pd.read_csv(path, sep="\t")
        if not df.empty:
            frames.append(df)
    n = len(a.samples)
    cols_out = ["site", "sample", "alias", "gene", "chrom", "strand", "side", "end", "pos", "location",
                "n_clipped", "depth", "clip_fraction", "median_clip_len", "max_clip_len", "consensus",
                "consensus_identity", "null_identity", "hap1_clipped", "hap2_clipped", "unphased_clipped",
                "n_samples_present", "n_other_samples_present", "n_samples", "outlier"]
    os.makedirs(os.path.dirname(a.all_outfile) or ".", exist_ok=True)
    if not frames:
        pd.DataFrame(columns=cols_out).to_csv(a.all_outfile, sep="\t", index=False)
        pd.DataFrame(columns=cols_out).to_csv(a.outliers_outfile, sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "No soft-clip clusters", ha="center", va="center")
        ax.axis("off")
        fig.savefig(a.plot_outfile)
        print("No soft-clip clusters in any sample.")
        return
    df = pd.concat(frames, ignore_index=True)
    df["sample"] = df["sample"].astype(str)
    # a site = gene + position + side
    df["site"] = df["gene"].astype(str) + ":" + df["chrom"].astype(str) + ":" + df["pos"].astype(str) + ":" + \
        df["side"].astype(str)
    df["alias"] = [resolve(s, alias) for s in df["sample"]]
    present = df["clip_fraction"] >= a.present_fraction
    df["n_samples_present"] = df["site"].map(df[present].groupby("site")["sample"].nunique()).fillna(0).astype(int)
    df["n_other_samples_present"] = df["n_samples_present"] - present.astype(int)
    df["n_samples"] = n
    df["consensus_identity"] = pd.to_numeric(df["consensus_identity"], errors="coerce")
    df["null_identity"] = pd.to_numeric(df["null_identity"], errors="coerce")
    # outlier needs enough clipped reads, a high clipped fraction, a consistent clip sequence,
    # and has to be rare in the rest of the group
    max_other = math.ceil(a.max_sample_fraction * n) - 1
    df["outlier"] = ((df["n_clipped"] >= a.min_reads) & (df["clip_fraction"] >= a.min_fraction)
                     & (df["consensus_identity"] >= a.min_identity) & (df["n_other_samples_present"] <= max_other))
    df = df[cols_out].sort_values(["outlier", "n_clipped"], ascending=[False, False])
    df.to_csv(a.all_outfile, sep="\t", index=False)
    out = df[df["outlier"]].sort_values(["sample", "n_clipped"], ascending=[True, False])
    out.to_csv(a.outliers_outfile, sep="\t", index=False)
    print(f"{len(df):,} sample cluster(s) at {df['site'].nunique():,} site(s); {len(out):,} outlier(s) "
          f"in {out['sample'].nunique()} of {n} sample(s) (at most {max_other} other sample(s) may have the site)")
    if not out.empty:
        print(out.groupby(["location", "end"]).size().rename("outliers").to_string())

    # clipped reads vs fraction, and clip identity vs a random-clip null
    use_fonts()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw=dict(width_ratios=[1.4, 1]))
    rec = df[~df["outlier"]]
    ax.scatter(rec["n_clipped"], rec["clip_fraction"], s=6, color=GREY, alpha=0.5, linewidths=0,
               label=f"not an outlier (n={len(rec):,})")
    ax.scatter(out["n_clipped"], out["clip_fraction"], s=14, color=RED, linewidths=0,
               label=f"outlier (n={len(out):,})")
    for _, r in out.nlargest(15, "n_clipped").iterrows():
        ax.annotate(f"{r['alias']} {r['gene']}", (r["n_clipped"], r["clip_fraction"]), fontsize=6,
                    xytext=(3, 2), textcoords="offset points")
    ax.axvline(a.min_reads, color="black", lw=0.6, ls="--")
    ax.axhline(a.min_fraction, color="black", lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("Clipped reads at cluster")
    ax.set_ylabel("Clipped fraction of depth")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.set_title(a.title or "Unexplained soft-clip clusters", loc="left", fontsize=9)
    bins = [i / 40 for i in range(41)]
    for col, color, lab in (("null_identity", GREY, "random clips (null)"), ("consensus_identity", RED, "clusters")):
        v = df[col].dropna()
        if not v.empty:
            bx.hist(v, bins=bins, color=color, alpha=0.6, label=f"{lab} (n={len(v):,})")
    bx.axvline(a.min_identity, color="black", lw=0.6, ls="--")
    bx.set_xlabel("Median identity of clips to their consensus")
    bx.set_ylabel("Clusters")
    bx.legend(frameon=False, fontsize=7)
    for x in (ax, bx):
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.plot_outfile)
    print(f"Saved {a.all_outfile}, {a.outliers_outfile}, {a.plot_outfile}")


if __name__ == "__main__":
    main()
