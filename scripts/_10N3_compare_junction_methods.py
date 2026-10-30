# 10N3 - GTEx-based vs in-house cohort junction calls for each sample (Sankey of junction status, UpSet of genes).
# most of the logic lives in junction_method_comparison.py so the final report uses the exact same code

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages

from palette import BLUE, GREY, use_fonts
from junction_method_comparison import (
    GT_OTHER,
    GT_READS,
    METRIC_EVENTS,
    NOT_SIG,
    NOT_TESTED,
    NT_OTHER,
    NT_SAMPLES,
    SIG_EVENT,
    SIG_NO_EVENT,
    cohort_status,
    computed_metrics,
    draw_sankey,
    event_set,
    gtex_status,
    panel_genes,
    read_tsv,
)

use_fonts()
rcParams["pdf.fonttype"] = 42
rcParams["font.size"] = 10


COHORT_SETS = [
    ("Cohort PSI_approx", ["junction_PSI_approx"]),
    ("Cohort PSI", ["junction_PSI"]),
    ("Cohort IR", ["5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio"]),
    ("Cohort IPA", ["junction_IPA_ratio"]),
]

MAX_COMBOS = 25


def parse_args():
    p = argparse.ArgumentParser(
        description="Compare the GTEx junction analysis, the cohort junction analysis and the cohort transcript "
                    "analysis for one sample-type group, using the calls that reach the merged hits. Each page is "
                    "one sample: a flow diagram of the junctions that are a PSI_approx event hit in GTEx or in the "
                    "cohort, from their GTEx status to their cohort status (with the reason a junction was not "
                    "tested), and an UpSet plot of the genes with hits by the analyses that found them. "
                    "The first page pools all samples.")
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--bed", required=True, help="Panel BED; column 4 is the gene symbol.")
    p.add_argument("--tissues", nargs="+", required=True, help="GTEx tissues of the group.")
    p.add_argument("--gtex-final-files", nargs="+", required=True,
                   help="10M final outlier junctions, one per --tissues entry (same order).")
    p.add_argument("--gtex-significant-files", nargs="+", required=True,
                   help="10M junctions passing coverage, padj and delta PSI (before the event filter), one per "
                        "--tissues entry.")
    p.add_argument("--gtex-sample-thresholds", nargs="+", type=int, required=True,
                   help="10M recurrence limit (max samples) for each --tissues entry.")
    p.add_argument("--gtex-all-junction-files", nargs="+", required=True,
                   help="Per-sample, per-tissue GTEx test results (_all_junctions.tsv).")
    p.add_argument("--gtex-all-junction-samples", nargs="+", required=True,
                   help="Sample of each --gtex-all-junction-files entry (same order).")
    p.add_argument("--gtex-all-junction-tissues", nargs="+", required=True,
                   help="Tissue of each --gtex-all-junction-files entry (same order).")
    p.add_argument("--cohort-results-dir", required=True, help="8B per-gene results directory.")
    p.add_argument("--cohort-outliers", required=True, help="8C outliers.tsv (all significant rows).")
    p.add_argument("--cohort-outliers-filtered", required=True, help="8C outliers_filtered.tsv (rows in events).")
    p.add_argument("--isoform-outliers", default=None, help="9G isoform outliers.")
    p.add_argument("--landscape-outliers", default=None, help="9G isoform-mix outliers.")
    p.add_argument("--min-jxn-reads", type=int, default=20, help="jxn_min_reads: reads a GTEx-novel junction needs in the sample to be tested.")
    p.add_argument("--coverage-threshold", type=int, default=50, help="8B coverage threshold (for labels).")
    p.add_argument("--n-threshold", type=int, default=10, help="8B minimum samples per junction (for labels).")
    p.add_argument("--title", default="Junction method comparison")
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", required=True)
    a = p.parse_args()
    for name in ("gtex_final_files", "gtex_significant_files", "gtex_sample_thresholds"):
        if len(getattr(a, name)) != len(a.tissues):
            p.error(f"--{name.replace('_', '-')} needs one entry per --tissues entry")
    if not (len(a.gtex_all_junction_files) == len(a.gtex_all_junction_samples) == len(a.gtex_all_junction_tissues)):
        p.error("--gtex-all-junction-files, --gtex-all-junction-samples and --gtex-all-junction-tissues must match")
    return a


# hit genes per method (GTEx and each cohort metric)
def gene_sets(a, gtex_sig_event, metrics):
    sets = {"GTEx PSI_approx": {(s, g) for s, g, _ in gtex_sig_event}}
    filt = read_tsv(a.cohort_outliers_filtered, ["sample", "gene", "event_type"])
    evs = event_set(filt["event_type"])
    for label, mets in COHORT_SETS:
        mets = [m for m in mets if m in metrics]
        if not mets:
            continue
        allowed = set().union(*(METRIC_EVENTS[m] for m in mets))
        hit = filt[evs.apply(lambda e: bool(e & allowed))]
        sets[label] = set(zip(hit["sample"], hit["gene"]))
    for label, path in (("Isoform usage", a.isoform_outliers), ("Isoform mix", a.landscape_outliers)):
        if path:
            df = read_tsv(path, ["sample", "gene", "qc"])
            df = df[df["qc"].astype(str).str.lower() != "fail"]
            sets[label] = set(zip(df["sample"], df["gene"]))
    return sets


def draw_upset(fig, spec, sample, sets):
    names = list(sets)
    genes = {}
    for name in names:
        for s, g in sets[name]:
            if s == sample:
                genes.setdefault(g, set()).add(name)
    gs = spec.subgridspec(2, 3, height_ratios=[2.2, 1], width_ratios=[0.7, 1.3, 3.2], hspace=0.05, wspace=0.02)
    ax_bar = fig.add_subplot(gs[0, 2])
    ax_mat = fig.add_subplot(gs[1, 2], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_mat)
    ax_lab = fig.add_subplot(gs[1, 1], sharey=ax_mat)
    ax_lab.axis("off")
    ax_bar.set_title(f"Genes with hits in {sample} (n={len(genes):,})", fontsize=11)

    counts = {}
    for v in genes.values():
        key = tuple(n for n in names if n in v)
        counts[key] = counts.get(key, 0) + 1
    order = sorted(counts, key=lambda c: (-counts[c], len(c), [names.index(n) for n in c]))
    hidden = order[MAX_COMBOS:]
    order = order[:MAX_COMBOS]

    y = {n: len(names) - 1 - i for i, n in enumerate(names)}
    if not order:
        ax_bar.text(0.5, 0.5, "No genes with hits", ha="center", va="center", transform=ax_bar.transAxes)
    x = np.arange(len(order))
    colors = [GREY if "GTEx PSI_approx" in c else BLUE for c in order]
    ax_bar.bar(x, [counts[c] for c in order], color=colors, width=0.7)
    for xi, c in zip(x, order):
        ax_bar.text(xi, counts[c], str(counts[c]), ha="center", va="bottom", fontsize=7)
    ax_bar.set_ylabel("Genes")
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ("top", "right"):
        ax_bar.spines[sp].set_visible(False)
    if order:
        ax_bar.set_ylim(0, max(counts[c] for c in order) * 1.15)
    if hidden:
        ax_bar.text(1, 1, f"{len(hidden)} smaller combination(s) not shown ({sum(counts[c] for c in hidden)} genes)",
                    ha="right", va="bottom", fontsize=7, transform=ax_bar.transAxes)
    if order:
        handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=GREY)]
        ax_bar.legend(handles, ["missed by GTEx", "found by GTEx"], frameon=False, fontsize=8, loc="upper right")
    else:
        ax_bar.set_ylim(0, 1)
        ax_bar.set_yticks([])

    for i in range(len(names)):
        if i % 2 == 0:
            ax_mat.axhspan(i - 0.5, i + 0.5, color="#f2f2f2", zorder=0)
    for xi, c in zip(x, order):
        ax_mat.scatter([xi] * len(names), [y[n] for n in names], s=40, color="#dddddd", zorder=2)
        ys = [y[n] for n in c]
        ax_mat.scatter([xi] * len(ys), ys, s=40, color="black", zorder=3)
        if len(ys) > 1:
            ax_mat.plot([xi, xi], [min(ys), max(ys)], color="black", lw=1.5, zorder=2)
    ax_mat.set_xlim(-0.6, max(len(order), 1) - 0.4)
    ax_mat.set_ylim(-0.5, len(names) - 0.5)
    ax_mat.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for n in names:
        ax_lab.text(0.5, y[n], n, ha="center", va="center", fontsize=8)
    for sp in ax_mat.spines.values():
        sp.set_visible(False)

    sizes = [sum(1 for v in genes.values() if n in v) for n in names]
    ax_set.barh([y[n] for n in names], sizes, color="#606060", height=0.6)
    ax_set.invert_xaxis()
    ax_set.tick_params(axis="y", left=False, labelleft=False)
    ax_set.set_xlabel("Genes", fontsize=8)
    ax_set.tick_params(axis="x", labelsize=7)
    for n, v in zip(names, sizes):
        ax_set.text(v, y[n], f"{v} ", ha="right", va="center", fontsize=7)
    if not any(sizes):
        ax_set.set_xlim(1, 0)
        ax_set.set_xticks([0])
    for sp in ("top", "left", "right"):
        ax_set.spines[sp].set_visible(False)
    return genes


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genes = panel_genes(a.bed)
    print("Reading GTEx junction results ...")
    gtex, gtex_sig_event = gtex_status(a, a.samples)
    print("Reading cohort junction results ...")
    cohort = cohort_status(a, a.samples, genes)
    metrics = computed_metrics(a)
    sets = gene_sets(a, gtex_sig_event, metrics)
    print(f"Gene sets: {', '.join(sets)}")

    labels = {
        "gtex": {SIG_EVENT: "significant, event\n(≥1 tissue)", SIG_NO_EVENT: "significant, no event\n(≥1 tissue)",
                 NOT_SIG: "not significant\n(all tissues)", NOT_TESTED: "not tested\n(≥1 tissue)",
                 GT_READS: f"not tested (any tissue):\n<{a.min_jxn_reads} reads in this sample",
                 GT_OTHER: "not tested (any tissue):\nother reason"},
        "cohort": {SIG_EVENT: "significant, event", SIG_NO_EVENT: "significant, no event",
                   NOT_SIG: "not significant",
                   NT_SAMPLES: f"not tested: <{a.n_threshold} samples\nwith coverage ≥{a.coverage_threshold}",
                   NT_OTHER: "not tested: other reason"},
    }
    # junctions that are an event hit by either method
    keys = sorted(set(gtex) | set(cohort))
    jdf = pd.DataFrame(keys, columns=["sample", "gene", "junction"])
    jdf["gtex_status"] = [gtex.get(k, GT_READS) for k in keys]
    jdf["cohort_status"] = [cohort.get(k, NT_OTHER) for k in keys]
    jdf = jdf[(jdf["gtex_status"] == SIG_EVENT) | (jdf["cohort_status"] == SIG_EVENT)].reset_index(drop=True)
    jpath = os.path.join(a.outdir, f"{a.prefix}_junction_status.tsv.gz")
    jdf.to_csv(jpath, sep="\t", index=False)
    print(f"Junction statuses -> {jpath}")

    def flows_of(df):
        return {k: int(v) for k, v in df.groupby(["gtex_status", "cohort_status"]).size().items()}

    note = ("Junctions that are a PSI_approx event hit in GTEx or in the cohort. Significant = reaches the merged "
            "hits (coverage, padj and delta PSI thresholds,\nreliable haplotypes, recurrence limit); "
            "no event = passes those thresholds but is not part of a splicing event. "
            "The dark blue flow is junctions that are an event hit in both.")
    flow_rows = []
    rows = []
    pdf_path = os.path.join(a.outdir, f"{a.prefix}_method_comparison.pdf")
    # a page per sample, plus the flow and gene-hit tables
    with PdfPages(pdf_path) as pdf:
        pooled = flows_of(jdf)
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_axes([0.04, 0.08, 0.92, 0.8])
        draw_sankey(ax, pooled, f"All {len(a.samples)} samples: {len(jdf):,} junction(s)", labels, fontsize=11)
        fig.suptitle(f"{a.title}: GTEx vs cohort junction analysis", fontsize=14, y=0.97)
        fig.text(0.5, 0.02, note, ha="center", fontsize=7, color="#555555")
        pdf.savefig(fig)
        plt.close(fig)
        for (l, r), v in pooled.items():
            flow_rows.append(("all", l, r, v))
        for sample in a.samples:
            sub = jdf[jdf["sample"] == sample]
            flows = flows_of(sub)
            for (l, r), v in flows.items():
                flow_rows.append((sample, l, r, v))
            fig = plt.figure(figsize=(11, 8.5))
            outer = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.2], left=0.0, right=0.98, top=0.86, bottom=0.08,
                                     wspace=0.08)
            draw_sankey(fig.add_subplot(outer[0, 0]), flows, f"Junctions with an event hit (n={len(sub):,})", labels)
            found = draw_upset(fig, outer[0, 1], sample, sets)
            fig.suptitle(f"{a.title}: {sample}", fontsize=14, y=0.97)
            fig.text(0.5, 0.015, note, ha="center", fontsize=7, color="#555555")
            pdf.savefig(fig)
            plt.close(fig)
            for g in sorted(found):
                row = {"sample": sample, "gene": g}
                row.update({n: n in found[g] for n in sets})
                row["combination"] = " & ".join(n for n in sets if n in found[g])
                rows.append(row)
    print(f"Plots -> {pdf_path}")
    fpath = os.path.join(a.outdir, f"{a.prefix}_junction_flows.tsv")
    pd.DataFrame(flow_rows, columns=["sample", "gtex_status", "cohort_status", "n_junctions"]).to_csv(
        fpath, sep="\t", index=False)
    print(f"Flow counts -> {fpath}")

    gdf = pd.DataFrame(rows, columns=["sample", "gene"] + list(sets) + ["combination"])
    for n in sets:
        gdf[n] = gdf[n].map({True: "TRUE", False: "FALSE"})
    gpath = os.path.join(a.outdir, f"{a.prefix}_gene_hits_by_method.tsv")
    gdf.to_csv(gpath, sep="\t", index=False)
    print(f"Gene hits by method -> {gpath}")


if __name__ == "__main__":
    main()
