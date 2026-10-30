# GTEx vs in-house junction method comparison, used by 10N3 and the final report

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

from palette import BLUE, GREY

# which event types each method can actually produce, and the status categories
GTEX_EVENTS = {"exon_skipping_approx", "exon_inclusion_approx", "alt_5ss_approx", "alt_3ss_approx", "complex_approx"}

METRIC_EVENTS = {
    "junction_PSI_approx":    {"alt_5ss_approx", "alt_3ss_approx", "exon_skipping_approx", "exon_inclusion_approx", "complex_approx"},
    "junction_PSI":           {"alt_5ss", "alt_3ss", "exon_skipping", "exon_inclusion", "complex"},
    "5ss_IR_ratio":           {"5ss_IR"},
    "3ss_IR_ratio":           {"3ss_IR"},
    "junction_full_IR_ratio": {"full_IR"},
    "junction_IPA_ratio":     {"IPA"},
}

COHORT_SETS = [
    ("Cohort PSI_approx", ["junction_PSI_approx"]),
    ("Cohort PSI", ["junction_PSI"]),
    ("Cohort IR", ["5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio"]),
    ("Cohort IPA", ["junction_IPA_ratio"]),
]

SIG_EVENT, SIG_NO_EVENT, NOT_SIG, NOT_TESTED = "significant, event", "significant, no event", "not significant", "not tested"
NT_SAMPLES = "not tested - too few samples with coverage"
NT_OTHER = "not tested - other reason"
GT_READS = "not tested in any tissue - too few reads in this sample"
GT_OTHER = "not tested in any tissue - other reason"
GTEX_STATUSES = [SIG_EVENT, SIG_NO_EVENT, NOT_SIG, NOT_TESTED, GT_READS, GT_OTHER]
COHORT_STATUSES = [SIG_EVENT, SIG_NO_EVENT, NOT_SIG, NT_SAMPLES, NT_OTHER]
COLORS = {
    SIG_EVENT: BLUE, SIG_NO_EVENT: "#42b4b5", NOT_SIG: GREY, NOT_TESTED: "#926d17",
    GT_READS: "#c69528", GT_OTHER: "#f4d77e", NT_SAMPLES: "#c69528", NT_OTHER: "#f4d77e",
}
MAX_COMBOS = 25


def read_tsv(path, cols):
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=lambda c: c in cols)
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    return df[cols]


def numeric(s):
    return pd.to_numeric(s, errors="coerce")


def is_true(s):
    return s.astype(str).str.strip().str.upper() == "TRUE"


def event_set(s):
    return s.fillna("").astype(str).apply(lambda v: {e.strip() for e in v.replace(";", ",").split(",") if e.strip()})


def panel_genes(bed):
    genes = []
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[3] not in genes:
                genes.append(f[3])
    return genes


# GTEx status per junction: significant event, significant only, not significant, or why it wasn't tested
def gtex_status(a, samples):
    samples = set(samples)
    tissues_of = {s: set() for s in samples}
    tested, other = {}, {}
    seen = set()
    for path, s, t in zip(a.gtex_all_junction_files, a.gtex_all_junction_samples, a.gtex_all_junction_tissues):
        if s not in samples:
            continue
        tissues_of[s].add(t)
        df = read_tsv(path, ["sample", "gene", "phasing", "junction", "padj", "jxn_alignment_count", "flag"])
        if df.empty:
            if not os.path.isfile(path):
                print(f"WARNING: GTEx all-junctions file missing: {path}")
            continue
        df = df[df["sample"].astype(str) == s]
        seen.update(zip(df["sample"], df["gene"], df["junction"]))
        ok = df[numeric(df["padj"]).notna()]
        for g, j in set(zip(ok["gene"], ok["junction"])):
            tested.setdefault((s, g, j), set()).add(t)
        bulk = df[df["phasing"] == "bulk"]
        low = (bulk["flag"].astype(str).str.contains("not_detected_in_gtex", regex=False)
               & (numeric(bulk["jxn_alignment_count"]).fillna(0) < a.min_jxn_reads))
        for g, j in set(zip(bulk.loc[~low, "gene"], bulk.loc[~low, "junction"])):
            other.setdefault((s, g, j), set()).add(t)

    sig_event, sig_no_event = set(), set()
    for t, final, sig, n in zip(a.tissues, a.gtex_final_files, a.gtex_significant_files, a.gtex_sample_thresholds):
        fdf = read_tsv(final, ["sample", "gene", "junction"])
        sig_event.update(zip(fdf["sample"], fdf["gene"], fdf["junction"]))
        sdf = read_tsv(sig, ["sample", "gene", "junction", "event", "sample_count"])
        no_ev = sdf[~event_set(sdf["event"]).apply(lambda e: bool(e & GTEX_EVENTS))
                    & (numeric(sdf["sample_count"]) <= n)]
        sig_no_event.update(zip(no_ev["sample"], no_ev["gene"], no_ev["junction"]))

    keys = seen | set(tested) | sig_event | sig_no_event
    out = {}
    for k in keys:
        if k[0] not in samples:
            continue
        if k in sig_event:
            out[k] = SIG_EVENT
        elif k in sig_no_event:
            out[k] = SIG_NO_EVENT
        elif tested.get(k, set()) >= tissues_of[k[0]] and tissues_of[k[0]]:
            out[k] = NOT_SIG
        elif tested.get(k):
            out[k] = NOT_TESTED
        elif other.get(k):
            out[k] = GT_OTHER
        else:
            out[k] = GT_READS
    return out, sig_event


# same idea for the in-house cohort analysis
def cohort_status(a, samples, genes):
    samples = set(samples)
    want = ("p_value_junction_PSI_approx", "modz_junction_PSI_approx")
    tested, seen, label = set(), set(), {}
    for gene in genes:
        path = os.path.join(a.cohort_results_dir, f"{gene}.tsv")
        if not os.path.isfile(path):
            continue
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                         usecols=lambda c: c in ("sample", "gene", "phasing", "junction") + want)
        if df.empty:
            continue
        if "gene" not in df.columns:
            df["gene"] = gene
        df = df[df["sample"].astype(str).isin(samples)]
        seen.update(zip(df["sample"], df["gene"], df["junction"]))
        col = next((c for c in want if c in df.columns), None)
        if col is not None:
            ok = df[numeric(df[col]).notna()]
            tested.update(zip(ok["sample"], ok["gene"], ok["junction"]))
            bulk = df[df["phasing"] == "bulk"] if "phasing" in df.columns else df
            label.update(zip(zip(bulk["sample"], bulk["gene"], bulk["junction"]), bulk[col].astype(str)))

    cols = ["sample", "gene", "junction", "event_type", "outlier_junction_PSI_approx"]
    filt = read_tsv(a.cohort_outliers_filtered, cols)
    allo = read_tsv(a.cohort_outliers, cols)
    in_event = filt[event_set(filt["event_type"]).apply(lambda e: bool(e & METRIC_EVENTS["junction_PSI_approx"]))]
    sig_event = set(zip(in_event["sample"], in_event["gene"], in_event["junction"]))
    sig_any = allo[is_true(allo["outlier_junction_PSI_approx"])]
    sig_no_event = set(zip(sig_any["sample"], sig_any["gene"], sig_any["junction"])) - sig_event

    out = {}
    for k in seen | tested | sig_event | sig_no_event:
        if k[0] not in samples:
            continue
        if k in sig_event:
            out[k] = SIG_EVENT
        elif k in sig_no_event:
            out[k] = SIG_NO_EVENT
        elif k in tested:
            out[k] = NOT_SIG
        elif label.get(k) == "low_n":
            out[k] = NT_SAMPLES
        else:
            out[k] = NT_OTHER
    return out


def computed_metrics(a):
    path = os.path.join(a.cohort_results_dir, "computed_metrics.txt")
    if os.path.isfile(path):
        with open(path) as fh:
            return [l.strip() for l in fh if l.strip()]
    if os.path.isfile(a.cohort_outliers):
        cols = pd.read_csv(a.cohort_outliers, sep="\t", nrows=0).columns
        return [m for m in METRIC_EVENTS if f"outlier_{m}" in cols]
    return []


def gene_sets(a, gtex_sig_event, metrics):
    sets = {"GTEx PSI_approx": {(s, g) for s, g, _ in gtex_sig_event}}
    filt = read_tsv(a.cohort_outliers_filtered, ["sample", "gene", "event_type"])
    drop = set(getattr(a, "drop_events", None) or [])
    evs = event_set(filt["event_type"]).apply(lambda e: e - drop)
    for label, mets in COHORT_SETS:
        mets = [m for m in mets if m in metrics]
        if not mets:
            continue
        allowed = set().union(*(METRIC_EVENTS[m] for m in mets)) - drop
        if not allowed:
            continue
        hit = filt[evs.apply(lambda e: bool(e & allowed))]
        sets[label] = set(zip(hit["sample"], hit["gene"]))
    for label, path in (("Isoform usage", a.isoform_outliers), ("Isoform mix", a.landscape_outliers)):
        if path:
            df = read_tsv(path, ["sample", "gene", "qc"])
            df = df[df["qc"].astype(str).str.lower() != "fail"]
            sets[label] = set(zip(df["sample"], df["gene"]))
    return sets


def _spread(ys, gap):
    order = np.argsort(ys)
    out = np.array(ys, dtype=float)
    for a, b in zip(order[:-1], order[1:]):
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    shift = max(0.0, out.max() - 1.02) if len(out) else 0.0
    return out - shift


# Sankey from GTEx status to cohort status
def draw_sankey(ax, flows, title, labels, fontsize=8):
    ax.set_xlim(-0.85, 1.85)
    ax.set_ylim(-0.08, 1.12)
    ax.axis("off")
    total = int(sum(flows.values()))
    both = int(flows.get((SIG_EVENT, SIG_EVENT), 0))
    ax.set_title(title, fontsize=11)
    if total == 0:
        ax.text(0.5, 0.5, "No event hits in either analysis", ha="center", va="center", fontsize=10)
        return
    left = {k: sum(v for (l, r), v in flows.items() if l == k) for k in GTEX_STATUSES}
    right = {k: sum(v for (l, r), v in flows.items() if r == k) for k in COHORT_STATUSES}
    left_nz = [k for k in GTEX_STATUSES if left[k]]
    right_nz = [k for k in COHORT_STATUSES if right[k]]
    gap = 0.035
    scale = min((1 - gap * (len(left_nz) - 1)) / total, (1 - gap * (len(right_nz) - 1)) / total)

    def stack(keys, sizes):
        span = sum(sizes[k] for k in keys) * scale + gap * (len(keys) - 1)
        top = 0.5 + span / 2
        pos = {}
        for k in keys:
            h = sizes[k] * scale
            pos[k] = (top - h, top)
            top -= h + gap
        return pos

    lpos, rpos = stack(left_nz, left), stack(right_nz, right)
    w = 0.05
    for pos, x, counts in ((lpos, 0.0, left), (rpos, 1.0 - w, right)):
        for k, (y0, y1) in pos.items():
            ax.add_patch(Rectangle((x, y0), w, y1 - y0, color=COLORS[k], lw=0, zorder=3))

    lcur = {k: lpos[k][1] for k in left_nz}
    rcur = {k: rpos[k][1] for k in right_nz}
    lslots, rslots = {}, {}
    for l in left_nz:
        for r in COHORT_STATUSES:
            v = flows.get((l, r), 0)
            if v:
                h = v * scale
                lslots[(l, r)] = (lcur[l] - h, lcur[l])
                lcur[l] -= h
    for r in right_nz:
        for l in GTEX_STATUSES:
            v = flows.get((l, r), 0)
            if v:
                h = v * scale
                rslots[(l, r)] = (rcur[r] - h, rcur[r])
                rcur[r] -= h
    x0, x1 = w, 1.0 - w
    xm = (x0 + x1) / 2
    for (l, r), (a0, a1) in lslots.items():
        b0, b1 = rslots[(l, r)]
        concordant = l == SIG_EVENT and r == SIG_EVENT
        color = BLUE if concordant else COLORS[r] if l == SIG_EVENT else COLORS[l]
        verts = [(x0, a1), (xm, a1), (xm, b1), (x1, b1), (x1, b0), (xm, b0), (xm, a0), (x0, a0), (x0, a1)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
                 Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="none",
                               alpha=0.75 if concordant else 0.45, zorder=2))

    tops = [p[SIG_EVENT][1] for p in (lpos, rpos) if SIG_EVENT in p]
    ax.text(xm, (max(tops) if tops else 1.0) + 0.025, f"shared event hits: {both:,}", ha="center", va="bottom",
            fontsize=fontsize + 1, fontweight="bold", color=BLUE if both else "#555555")

    for pos, counts, xs, ha, lab in ((lpos, left, -0.03, "right", labels["gtex"]),
                                     (rpos, right, 1.03, "left", labels["cohort"])):
        keys = list(pos)
        mids = [(pos[k][0] + pos[k][1]) / 2 for k in keys]
        ys = _spread(mids, 0.075 * fontsize / 8)
        for k, m, y in zip(keys, mids, ys):
            n_txt = f"({counts[k]:,}; {both:,} shared)" if k == SIG_EVENT else f"({counts[k]:,})"
            ax.text(xs + (-0.04 if ha == "right" else 0.04), y, f"{lab[k]}  {n_txt}", ha=ha, va="center",
                    fontsize=fontsize, fontweight="bold" if k == SIG_EVENT else "normal")
            ax.plot([xs, xs + (-0.04 if ha == "right" else 0.04)], [m, y], color="#888888", lw=0.6)
    ax.text(w / 2, 1.09, "GTEx", ha="center", va="bottom", fontsize=fontsize + 2, fontweight="bold")
    ax.text(1 - w / 2, 1.09, "Cohort", ha="center", va="bottom", fontsize=fontsize + 2, fontweight="bold")


def draw_upset(fig, spec, sample, sets, shown=None):
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
    ax_bar.set_title(f"Genes with hits in {shown or sample} (n={len(genes):,})", fontsize=11)

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


def status_labels(a):
    return {
        "gtex": {SIG_EVENT: "significant, event\n(≥1 tissue)", SIG_NO_EVENT: "significant, no event\n(≥1 tissue)",
                 NOT_SIG: "not significant\n(all tissues)", NOT_TESTED: "not tested\n(≥1 tissue)",
                 GT_READS: f"not tested (any tissue):\n<{a.min_jxn_reads} reads in this sample",
                 GT_OTHER: "not tested (any tissue):\nother reason"},
        "cohort": {SIG_EVENT: "significant, event", SIG_NO_EVENT: "significant, no event",
                   NOT_SIG: "not significant",
                   NT_SAMPLES: f"not tested: <{a.n_threshold} samples\nwith coverage ≥{a.coverage_threshold}",
                   NT_OTHER: "not tested: other reason"},
    }


NOTE = ("Junctions that are a PSI_approx event hit in GTEx or in the cohort. Significant = reaches the merged "
        "hits (coverage, padj and delta PSI thresholds,\nreliable haplotypes, recurrence limit); "
        "no event = passes those thresholds but is not part of a splicing event. "
        "The dark blue flow is junctions that are an event hit in both.")


# flows + gene sets for one sample, and the combined figure for the report
def sample_comparison(a, sample):
    sample = str(sample)
    gtex, gtex_sig_event = gtex_status(a, [sample])
    cohort = cohort_status(a, [sample], panel_genes(a.bed))
    keys = sorted(set(gtex) | set(cohort))
    flows = {}
    n = 0
    for k in keys:
        g, c = gtex.get(k, GT_READS), cohort.get(k, NT_OTHER)
        if g == SIG_EVENT or c == SIG_EVENT:
            flows[(g, c)] = flows.get((g, c), 0) + 1
            n += 1
    sets = {name: {x for x in members if str(x[0]) == sample}
            for name, members in gene_sets(a, gtex_sig_event, computed_metrics(a)).items()}
    return flows, n, sets


def draw_comparison(fig, rect, a, sample, shown=None):
    flows, n, sets = sample_comparison(a, sample)
    x0, y0, x1, y1 = rect
    outer = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.2], left=x0, right=x1, top=y1, bottom=y0, wspace=0.08)
    draw_sankey(fig.add_subplot(outer[0, 0]), flows, f"Junctions with an event hit (n={n:,})", status_labels(a))
    draw_upset(fig, outer[0, 1], str(sample), sets, shown=shown)
    return n
