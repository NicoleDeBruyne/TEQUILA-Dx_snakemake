"""
plot_junction_distributions.py

For each junction that appears as a hit in merged_all_hits.tsv, plot the cohort-wide
distribution of the relevant metric across all samples, split by sample type.

All metric values (including junction_PSI_approx from GTEx hits) are sourced from
the per-gene raw TSVs produced by 8A (via the per-group manifest). GTEx hits and
cohort hits for junction_PSI_approx are unioned — the outlier labels reflect both.

Only bulk phasing rows are plotted. Only samples with denominator >= coverage_threshold
are included. Outlier samples (those appearing in the hits table for that junction×metric)
are labelled on the plot.

Output: one PDF per (metric × junction) at:
  {outdir}/{metric}_distributions/{junction}.pdf
where outdir is already the filter-set-specific directory ({cohort_outdir}/{bed_id}/output/{filter_set}).
A sentinel file {outdir}/junction_distributions.done is written when complete.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import traceback
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Column names used in merged_all_hits.tsv
# ---------------------------------------------------------------------------

# GTEx PSI_approx junction hits live in these columns (bulk + haplotypes)
_GTEX_JXN_COLS = ["bulk_jxns", "hap1_jxns", "hap2_jxns"]

# Cohort junction hit columns, keyed by delta column suffix → metric name
# The delta columns encode which metric fired; we use the jxns column to get
# the junction IDs and the delta column to confirm a non-empty hit.
_COHORT_PHASING_PREFIXES = ["cohort_bulk", "cohort_hap1", "cohort_hap2"]
_COHORT_DELTA_SUFFIX_TO_METRIC: Dict[str, str] = {
    "deltaPSI":       "junction_PSI",
    "deltaPSIapprox": "junction_PSI_approx",
    "delta5ssIR":     "5ss_IR_ratio",
    "delta3ssIR":     "3ss_IR_ratio",
    "deltaFullIR":    "junction_full_IR_ratio",
    "deltaIPA":       "junction_IPA_ratio",
}

# For each cohort metric: which column in the 8A per-gene TSV holds the rescaled value,
# and which column is the denominator (must be >= coverage_threshold for sample to count).
_METRIC_RESCALED_COL: Dict[str, str] = {
    "junction_PSI_approx":    "rescaled_junction_PSI_approx",
    "junction_PSI":           "rescaled_junction_PSI",
    "5ss_IR_ratio":           "rescaled_5ss_IR_ratio",
    "3ss_IR_ratio":           "rescaled_3ss_IR_ratio",
    "junction_full_IR_ratio": "rescaled_junction_full_IR_ratio",
    "junction_IPA_ratio":     "rescaled_junction_IPA_ratio",
}
_METRIC_DENOM_COL: Dict[str, str] = {
    "junction_PSI_approx":    "junction_coverage_approx",
    "junction_PSI":           "junction_coverage",
    "5ss_IR_ratio":           "5ss_coverage",
    "3ss_IR_ratio":           "3ss_coverage",
    "junction_full_IR_ratio": "junction_coverage",
    "junction_IPA_ratio":     "5ss_coverage",
}

_MISSING = {"", ".", "nan", "none", "n/a"}


def _nonempty(cell) -> bool:
    if pd.isna(cell):
        return False
    return str(cell).strip().lower() not in _MISSING


def _split_jxns(cell) -> List[str]:
    if not _nonempty(cell):
        return []
    return [j.strip() for j in re.split(r"[,;]", str(cell))
            if j.strip() and j.strip().lower() not in _MISSING]


def _safe_filename(jxn: str) -> str:
    """Make junction ID safe for use as a filename."""
    return re.sub(r"[^\w.\-]", "_", jxn)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot cohort junction metric distributions for hit junctions."
    )
    p.add_argument("--hits-tsv",         required=True,
                   help="merged_all_hits.tsv produced by rule _9F_final_merge.")
    p.add_argument("--outdir",           required=True,
                   help="Base output directory ({cohort_outdir}/{bed_id}/output/{filter_set}).")
    p.add_argument("--filter-set",       required=True,
                   help="'default' or 'stringent'.")
    # 8A manifests: space-separated list of "sample_type:manifest_path" pairs
    p.add_argument("--cohort-manifests", nargs="*", default=[],
                   metavar="SAMPLE_TYPE:PATH",
                   help="Per-group 8A manifest files: 'sample_type:path' pairs.")
    # Sample metadata
    p.add_argument("--samples",          nargs="+", required=True,
                   help="All sample IDs in this cohort (same order as --sample-types).")
    p.add_argument("--sample-types",     nargs="+", required=True,
                   help="Sample type for each sample in --samples.")
    p.add_argument("--coverage-threshold", type=int, default=50,
                   help="Minimum denominator coverage for a sample to be included in the plot.")
    p.add_argument("--sentinel",         required=True,
                   help="Path to write sentinel file on completion.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Extract hits from merged_all_hits.tsv
# ---------------------------------------------------------------------------

def extract_hits(
    hits_df: pd.DataFrame,
) -> Dict[Tuple[str, str], Set[str]]:
    """
    Returns a single dict mapping (metric, junction) -> set of outlier sample IDs.

    GTEx PSI_approx hits (bulk_deltaPSI / hap_deltaPSI columns) and cohort
    junction_PSI_approx hits are unioned under the same metric key so that
    outlier labels reflect both sources. All values are loaded from 8A raw TSVs.
    """
    hits: Dict[Tuple[str, str], Set[str]] = defaultdict(set)

    for _, row in hits_df.iterrows():
        sample = str(row["sample"])

        # GTEx PSI_approx hits → map to junction_PSI_approx (values from 8A)
        if _nonempty(row.get("bulk_deltaPSI", "")):
            for jxn in _split_jxns(row.get("bulk_jxns", "")):
                hits[("junction_PSI_approx", jxn)].add(sample)
        for hap in ("hap1", "hap2"):
            if _nonempty(row.get(f"{hap}_deltaPSI", "")):
                for jxn in _split_jxns(row.get(f"{hap}_jxns", "")):
                    hits[("junction_PSI_approx", jxn)].add(sample)

        # Cohort hits: all phasing prefixes × all delta columns
        for prefix in _COHORT_PHASING_PREFIXES:
            jxns_col = f"{prefix}_jxns"
            jxns = _split_jxns(row.get(jxns_col, ""))
            if not jxns:
                continue
            for delta_suffix, metric in _COHORT_DELTA_SUFFIX_TO_METRIC.items():
                delta_col = f"{prefix}_{delta_suffix}"
                if _nonempty(row.get(delta_col, "")):
                    for jxn in jxns:
                        hits[(metric, jxn)].add(sample)

    return hits


# ---------------------------------------------------------------------------
# Bulk load raw metric values (read each file once across all junctions)
# ---------------------------------------------------------------------------

def load_all_cohort_values(
    manifests: Dict[str, str],
    needed_pairs: Set[Tuple[str, str]],   # {(metric, junction)}
    coverage_threshold: int,
) -> Dict[Tuple[str, str], Dict[str, float]]:
    """
    Read every per-gene 8A raw TSV once and accumulate values for all requested
    (metric, junction) pairs.  Returns {(metric, junction): {sample: value}}.

    Only bulk phasing rows with denominator >= coverage_threshold are kept.
    """
    # Pre-compute which (rescaled_col, denom_col) pairs we need per junction
    # Group needed_pairs by junction so we can handle multiple metrics per file read.
    jxn_to_metric_cols: Dict[str, List[Tuple[str, str, str]]] = defaultdict(list)
    for (metric, jxn) in needed_pairs:
        rescaled_col = _METRIC_RESCALED_COL.get(metric)
        denom_col    = _METRIC_DENOM_COL.get(metric)
        if rescaled_col and denom_col:
            jxn_to_metric_cols[jxn].append((metric, rescaled_col, denom_col))

    needed_junctions: Set[str] = set(jxn_to_metric_cols.keys())
    result: Dict[Tuple[str, str], Dict[str, float]] = defaultdict(dict)

    for sample_type, manifest_path in manifests.items():
        if not os.path.isfile(manifest_path):
            continue
        try:
            manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
        except Exception as e:
            print(f"  [WARNING] Could not read manifest {manifest_path}: {e}")
            continue

        n_files = len(manifest)
        for i, (_, mrow) in enumerate(manifest.iterrows(), 1):
            result_path = mrow.get("result_path")
            if pd.isna(result_path) or str(result_path) in ("None", ""):
                continue
            if not os.path.isfile(str(result_path)):
                continue
            try:
                df = pd.read_csv(str(result_path), sep="\t", dtype=str)
            except Exception as e:
                print(f"  [WARNING] Could not read {result_path}: {e}")
                continue

            if "junction" not in df.columns or "phasing" not in df.columns or "sample" not in df.columns:
                continue

            # Filter to bulk rows and only junctions we care about
            df = df[df["phasing"] == "bulk"]
            df = df[df["junction"].isin(needed_junctions)]
            if df.empty:
                continue

            for jxn, grp in df.groupby("junction"):
                for (metric, rescaled_col, denom_col) in jxn_to_metric_cols.get(str(jxn), []):
                    if rescaled_col not in grp.columns or denom_col not in grp.columns:
                        continue
                    vals   = grp[rescaled_col].pipe(pd.to_numeric, errors="coerce")
                    denoms = grp[denom_col].pipe(pd.to_numeric, errors="coerce")
                    mask   = denoms >= coverage_threshold
                    for sample, val in zip(grp["sample"][mask], vals[mask]):
                        if sample and np.isfinite(val):
                            result[(metric, str(jxn))][str(sample)] = float(val)

            if i % 500 == 0:
                print(f"    [{sample_type}] {i}/{n_files} gene files scanned ...")

    return result




# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_junction_metric(
    metric:             str,
    junction:           str,
    sample_values:      Dict[str, float],   # {sample: value}
    sample_to_type:     Dict[str, str],     # {sample: sample_type}
    outlier_samples:    Set[str],
    outpath:            str,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Group values by sample type
    type_order = sorted(set(sample_to_type.get(s, "unknown") for s in sample_values))
    type_data:    Dict[str, List[float]] = {t: [] for t in type_order}
    type_samples: Dict[str, List[str]]  = {t: [] for t in type_order}

    for sample, val in sample_values.items():
        stype = sample_to_type.get(sample, "unknown")
        if stype not in type_data:
            type_data[stype]    = []
            type_samples[stype] = []
            type_order.append(stype)
        type_data[stype].append(val)
        type_samples[stype].append(sample)

    if not type_order:
        return

    n_types  = len(type_order)
    fig_w    = max(4.0, n_types * 1.8)
    fig, ax  = plt.subplots(figsize=(fig_w, 4.5))

    positions = list(range(1, n_types + 1))
    box_data  = [type_data[t] for t in type_order]

    bp = ax.boxplot(
        box_data,
        positions=positions,
        widths=0.5,
        showfliers=False,
        patch_artist=False,
        boxprops=dict(color="black"),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        medianprops=dict(color="#c0392b", linewidth=1.5),
    )

    rng = np.random.default_rng(seed=42)
    for pos, stype in zip(positions, type_order):
        vals    = type_data[stype]
        samples = type_samples[stype]
        if not vals:
            continue
        xs = rng.normal(pos, 0.06, size=len(vals))
        for x, y, s in zip(xs, vals, samples):
            is_hit = s in outlier_samples
            ax.scatter(
                x, y,
                color="#c0392b" if is_hit else "#2c7fb8",
                s=20,
                alpha=0.85,
                zorder=3,
                edgecolors="black" if is_hit else "none",
                linewidths=0.8 if is_hit else 0,
            )
            if is_hit:
                ax.annotate(
                    s, (x, y),
                    fontsize=5,
                    xytext=(4, 0),
                    textcoords="offset points",
                    va="center",
                )

    ax.set_xticks(positions)
    ax.set_xticklabels(type_order, fontsize=9)
    ax.set_ylabel(metric, fontsize=9)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(0.3, n_types + 0.7)

    ax.set_title(junction, fontsize=8, pad=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with PdfPages(outpath) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    if len(args.samples) != len(args.sample_types):
        sys.exit("ERROR: --samples and --sample-types must have equal length.")

    sample_to_type = dict(zip(args.samples, args.sample_types))

    # Parse manifest pairs: "sample_type:path"
    manifests: Dict[str, str] = {}
    for item in args.cohort_manifests:
        stype, _, path = item.partition(":")
        if path:
            manifests[stype.strip()] = path.strip()

    print(f"\n{'='*70}")
    print("  Junction Distribution Plotter")
    print(f"{'='*70}")
    print(f"  Filter set:          {args.filter_set}")
    print(f"  Coverage threshold:  {args.coverage_threshold}")
    print(f"  Cohort manifests:    {len(manifests)}")
    print(f"  Samples:             {len(args.samples)}")

    if not os.path.isfile(args.hits_tsv):
        print(f"\n[WARNING] hits TSV not found: {args.hits_tsv}")
        print("Writing empty sentinel.")
        os.makedirs(os.path.dirname(args.sentinel), exist_ok=True)
        open(args.sentinel, "w").close()
        return

    hits_df = pd.read_csv(args.hits_tsv, sep="\t", dtype=str)
    if hits_df.empty:
        print("\n[INFO] No hits in hits TSV. Nothing to plot.")
        os.makedirs(os.path.dirname(args.sentinel), exist_ok=True)
        open(args.sentinel, "w").close()
        return

    print(f"\nExtracting hits from {args.hits_tsv} ({len(hits_df)} rows)...")
    all_pairs = extract_hits(hits_df)

    print(f"Found {len(all_pairs)} (metric × junction) pairs to plot.")

    outdir_fs = args.outdir

    # Check which pairs already have a PDF on disk and skip them
    already_done: List[Tuple[str, str]] = []
    remaining: List[Tuple[str, str]] = []
    for (metric, junction) in all_pairs:
        safe_jxn = _safe_filename(junction)
        outpath = os.path.join(outdir_fs, f"{metric}_distributions", f"{safe_jxn}.pdf")
        if os.path.isfile(outpath):
            already_done.append((metric, junction))
        else:
            remaining.append((metric, junction))

    print(f"{len(already_done)} (metric × junction) pairs already plotted.")
    print(f"{len(remaining)} (metric × junction) pairs remain.\n")

    if not remaining:
        print("Nothing to do.")
    else:
        # --- Bulk load all values from 8A gene TSVs (one pass) ---
        print("Loading 8A raw values (one pass over all gene TSVs) ...")
        all_values = load_all_cohort_values(manifests, set(remaining), args.coverage_threshold)
        print(f"  Loaded values for {len(all_values)} (metric × junction) pairs.\n")

        # --- Plot ---
        n_written = 0
        n_failed  = 0

        for (metric, junction) in sorted(remaining):
            outlier_samples = all_pairs[(metric, junction)]
            try:
                sample_values = all_values.get((metric, junction), {})

                if not sample_values:
                    print(f"  [SKIP] No valid values for {metric} / {junction}")
                    continue

                safe_jxn = _safe_filename(junction)
                outpath  = os.path.join(outdir_fs, f"{metric}_distributions", f"{safe_jxn}.pdf")

                print(f"  Plotting {metric} / {junction} "
                      f"({len(sample_values)} samples, {len(outlier_samples)} outlier(s)) ...")

                plot_junction_metric(
                    metric=metric,
                    junction=junction,
                    sample_values=sample_values,
                    sample_to_type=sample_to_type,
                    outlier_samples=outlier_samples,
                    outpath=outpath,
                )
                n_written += 1

            except Exception as e:
                print(f"  [ERROR] {metric} / {junction}: {e}")
                traceback.print_exc()
                n_failed += 1

        print(f"\nDone. Wrote {n_written} PDF(s), {n_failed} failure(s).")

    os.makedirs(os.path.dirname(args.sentinel), exist_ok=True)
    open(args.sentinel, "w").close()
    print(f"Sentinel written: {args.sentinel}")


if __name__ == "__main__":
    main()
