#!/usr/bin/env python3
"""
scripts/identify_cohort_junction_outliers.py

Statistical outlier-detection stage for cohort-level splice-junction
analysis. Companion script scripts/cohort_junction_analysis.py computes the
upstream per-gene, per-sample junction coverage/usage metrics this script
reads.

Reads the gene -> result-file manifest produced by cohort_junction_analysis.py,
loads each gene's raw per-junction/sample/phasing metrics, fits a Beta
distribution per junction across the cohort's bulk samples, runs a
beta-binomial test per sample/junction against that distribution, applies
FDR (Benjamini-Hochberg) correction across the whole cohort, and identifies
+ classifies outlier junctions for one or more padj:delta threshold pairs.
"""

from __future__ import annotations

import os
import sys
import argparse
import warnings
import traceback
import time
import math
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import concurrent.futures
from scipy.stats import betabinom, beta
from statsmodels.stats.multitest import multipletests
from pandas.errors import PerformanceWarning

warnings.filterwarnings("ignore", category=PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*Glyph.*missing from font.*", category=UserWarning)
warnings.filterwarnings("ignore", message=".*Adding colorbar to a different Figure.*", category=UserWarning)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Maps each metric to the event types it can produce
_METRIC_EVENTS: Dict[str, List[str]] = {
    "junction_PSI_approx":    ["alt_5ss_approx", "alt_3ss_approx", "exon_skipping_approx", "exon_inclusion_approx"],
    "junction_PSI":           ["alt_5ss", "alt_3ss", "exon_skipping", "exon_inclusion"],
    "5ss_IR_ratio":           ["5ss_IR"],
    "3ss_IR_ratio":           ["3ss_IR"],
    "junction_full_IR_ratio": ["full_IR"],
    "junction_IPA_ratio":     ["IPA"],
}

# Columns in a raw per-gene TSV (from cohort_junction_analysis.py) that are
# booleans on write but come back as text after a TSV round-trip.
_BOOL_COLUMN_PREFIXES = ("low_phased_",)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Identifies splice junctions with unusual usage frequencies "
                     "from the per-gene metrics computed by cohort_junction_analysis.py."
    )
    p.add_argument("--manifest",                   required=True,
                   help="gene -> result-file manifest TSV from cohort_junction_analysis.py.")
    p.add_argument("--bed",                        required=True)
    p.add_argument("--outprefix",                  required=True)
    p.add_argument("--approx",                     action="store_true")
    p.add_argument("--has-ipa",                     action="store_true",
                   help="Set if --genome was provided to cohort_junction_analysis.py "
                        "(enables testing junction_IPA_ratio as a metric).")
    p.add_argument("--coverage-threshold",         type=int,   default=20)
    p.add_argument("--PSI-rescale-factor",         type=float, default=1e-3)
    p.add_argument("--n-threshold",                type=int,   default=30)
    method_group = p.add_mutually_exclusive_group(required=True)
    method_group.add_argument(
        "--bb-thresholds",
        nargs="*", default=None,
        metavar="PADJ:DELTA",
        help="Use beta-binomial testing (per-junction Beta distribution + "
             "beta-binomial test, FDR-corrected). One or more padj:delta "
             "threshold pairs, e.g. 0.05:0.1 0.01:0.1 0.01:0.2 -- each "
             "combination produces its own output subdirectory. Pass with "
             "no values to compute beta-binomial statistics (columns "
             "n/alpha/beta/expected/p1/p99/delta/p_value/padj per metric) "
             "without identifying outliers.")
    method_group.add_argument(
        "--z-thresholds",
        nargs="*", default=None,
        metavar="MODZ",
        help="Use modified z-score testing (per-junction median/MAD, no "
             "p-values or FDR correction). One or more |modZ| cutoffs, e.g. "
             "3.5 5 -- each produces its own output subdirectory. Defaults "
             "to 3.5 if the flag is given with no values.")
    p.add_argument("--no-ss-IR",                   action="store_true")
    p.add_argument("--gtf",                        default=None,
                   help="GTF/GTF.gz. If provided, adds junction_type column.")
    p.add_argument("--threads",                    type=int,   default=1)
    p.add_argument("--test-n-genes",               type=int,   default=None)
    return p.parse_args()


# ---------------------------------------------------------------------------
# BED / manifest I/O
# ---------------------------------------------------------------------------

def load_bed(path: str) -> Dict[str, Tuple[str, str, str]]:
    gene_info: Dict[str, Tuple[str, str, str]] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                print(f"[WARNING] BED line has <6 columns, skipping: {line!r}")
                continue
            chrom, start, end, gene, _, strand = parts[:6]
            gene   = gene.strip()
            strand = strand.strip()
            region = f"{chrom.strip()}:{start.strip()}-{end.strip()}"
            if gene in gene_info:
                print(f"[WARNING] Gene '{gene}' appears more than once in BED. Using last entry.")
            gene_info[gene] = (chrom.strip(), region, strand)
    print(f"BED file loaded: {len(gene_info)} gene(s)")
    return gene_info


def load_manifest(path: str) -> List[Tuple[str, Optional[str]]]:
    """
    Reads the gene -> result-file manifest from cohort_junction_analysis.py.
    Returns a list of (gene, path_or_None) in the manifest's own row order
    (which is the input BED file's order).
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    rows: List[Tuple[str, Optional[str]]] = []
    for _, row in df.iterrows():
        gene = row["gene"]
        path_val = row.get("result_path")
        rows.append((gene, None if (pd.isna(path_val) or path_val == "None") else path_val))
    n_with_data = sum(1 for _, p in rows if p is not None)
    print(f"Manifest loaded: {len(rows)} gene(s), {n_with_data} with results")
    return rows


def load_gene_raw_metrics(path: str) -> pd.DataFrame:
    """Loads one gene's raw metrics TSV (written by cohort_junction_analysis.py),
    restoring boolean dtype for low_phased_* columns lost in the TSV round-trip."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in df.columns:
        if any(col.startswith(p) for p in _BOOL_COLUMN_PREFIXES):
            df[col] = (df[col].astype(str).str.strip().str.lower()
                       .map({"true": True, "false": False}).fillna(False))
        else:
            # Try numeric conversion; columns that are genuinely text
            # (junction, 5ss, 3ss, phasing, sample, region, gene, and the
            # object-typed metric/rescaled columns with "low_coverage"
            # sentinels) simply fail conversion and stay as strings, which
            # downstream code already handles via pd.to_numeric(errors="coerce").
            # (errors="ignore" was removed in newer pandas -- try/except is
            # the version-portable equivalent.)
            try:
                df[col] = pd.to_numeric(df[col])
            except (ValueError, TypeError):
                pass
    return df


# ---------------------------------------------------------------------------
# Beta distribution fitting
# ---------------------------------------------------------------------------

def _fit_one_beta(x: np.ndarray, tol: float, n_threshold: int):
    x = x[np.isfinite(x)]
    n = len(x)
    if n < n_threshold: return n, "low_n", "low_n", "low_n"
    var = float(np.var(x))
    if var < tol:
        m = float(np.clip(np.median(x), tol, 1.0 - tol))
        k = max(m * (1.0 - m) / tol - 1.0, tol)
        a = m * k; b = (1.0 - m) * k
        return n, float(a), float(b), float(m)
    try:
        a, b, *_ = beta.fit(x, floc=0, fscale=1)
        return n, float(a), float(b), float(a / (a + b))
    except Exception:
        return n, "error", "error", "error"


def _fit_beta_rows(args):
    """Worker: fit a sub-block of rows. Returns list of result tuples."""
    mat_block, tol, n_threshold = args
    return [_fit_one_beta(mat_block[i], tol, n_threshold)
            for i in range(len(mat_block))]


def fit_beta_dist_chunk(mat, feat_names, tol, n_threshold, threads: int = 1):
    """
    Fit beta distributions for every row in mat.

    Parallelised across rows using a process pool when threads > 1.
    Results are collected via executor.map, which preserves submission order,
    so output is deterministic regardless of the number of workers.
    """
    n = len(mat)
    if n == 0:
        return pd.DataFrame(columns=["n", "alpha", "beta_param", "expected"])

    n_workers = min(threads, n)
    if n_workers <= 1:
        results = [_fit_one_beta(mat[i], tol, n_threshold) for i in range(n)]
    else:
        chunk_size = max(1, (n + n_workers - 1) // n_workers)
        chunks = [
            (mat[i : i + chunk_size], tol, n_threshold)
            for i in range(0, n, chunk_size)
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
            # map() preserves submission order — results are deterministic
            results = []
            for block in ex.map(_fit_beta_rows, chunks):
                results.extend(block)

    return pd.DataFrame(results, index=feat_names,
                        columns=["n", "alpha", "beta_param", "expected"])


# ---------------------------------------------------------------------------
# Beta-binomial test
# ---------------------------------------------------------------------------

def _betabinom_test_rows(args):
    """
    Worker: run betabinom test on a sub-DataFrame.
    Receives a dict of arrays (picklable) rather than a DataFrame.
    Returns a 1-D array of p-values (NaN where not testable).
    """
    usage, coverage, alpha_v, beta_v, val, cov_thresh = args
    valid = (np.isfinite(usage) & np.isfinite(coverage) &
             np.isfinite(alpha_v) & np.isfinite(beta_v) & np.isfinite(val) &
             (coverage >= cov_thresh))
    p = np.full(len(usage), np.nan)
    if valid.any():
        u   = np.round(usage[valid]).astype(int)
        cv  = np.round(coverage[valid]).astype(int)
        lte = betabinom.cdf(u, cv, alpha_v[valid], beta_v[valid])
        gte = betabinom.cdf(cv - u, cv, beta_v[valid], alpha_v[valid])
        p[valid] = np.clip(2.0 * np.minimum(lte, gte), 0.0, 1.0)
    return p


def beta_binomial_test_chunk(
    df: pd.DataFrame,
    coverage_threshold: int,
    metric_col: str,
    usage_col: str,
    coverage_col: str,
    p_col: str,
    threads: int = 1,
) -> pd.DataFrame:
    """
    Run betabinom tests for every row in df.

    Parallelised across rows using a process pool when threads > 1.
    Results are collected via executor.map (order-preserving) so p-values
    are identical regardless of worker count.
    """
    df = df.copy()
    usage    = df[usage_col].to_numpy(dtype=float)
    coverage = df[coverage_col].to_numpy(dtype=float)
    alpha_v  = pd.to_numeric(df[f"alpha_{metric_col}"], errors="coerce").to_numpy()
    beta_v   = pd.to_numeric(df[f"beta_{metric_col}"],  errors="coerce").to_numpy()
    val      = df[metric_col].to_numpy(dtype=float)

    n = len(df)
    n_workers = min(threads, n)

    if n_workers <= 1:
        p = _betabinom_test_rows(
            (usage, coverage, alpha_v, beta_v, val, coverage_threshold)
        )
    else:
        chunk_size = max(1, (n + n_workers - 1) // n_workers)
        chunks = [
            (
                usage   [i : i + chunk_size],
                coverage[i : i + chunk_size],
                alpha_v [i : i + chunk_size],
                beta_v  [i : i + chunk_size],
                val     [i : i + chunk_size],
                coverage_threshold,
            )
            for i in range(0, n, chunk_size)
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
            # map() preserves submission order — p-values are deterministic
            p = np.concatenate(list(ex.map(_betabinom_test_rows, chunks)))

    df[p_col] = p
    return df


# ---------------------------------------------------------------------------
# Modified z-score fitting (median/MAD)
# ---------------------------------------------------------------------------
#
# Alternative to the Beta-binomial approach above: for each junction, take
# the bulk cohort's median and median absolute deviation (MAD) of its
# rescaled metric value, then score every row (bulk and hap1/hap2) as a
# modified z-score: 0.6745 * (value - median) / MAD (the 0.6745 constant
# makes MAD-based spread comparable to a standard deviation for normally
# distributed data -- the standard Iglewicz & Hoaglin formulation).
#
# No p-value or FDR correction is computed for this method -- outliers are
# identified by a direct |modZ| cutoff (see --z-thresholds).

_MODZ_CONST = 0.6745        # MAD -> z-score scaling constant
_MEANAD_CONST = 0.7979      # mean-absolute-deviation -> z-score scaling constant


def _fit_one_modz(x: np.ndarray, tol: float, n_threshold: int):
    """Returns (n, median, effective_mad) for one junction's bulk values.

    effective_mad is always on the same scale as MAD (i.e. modZ = 0.6745 *
    (x - median) / effective_mad), even when the true MAD is degenerate:
      - if MAD is essentially zero but the mean absolute deviation isn't,
        falls back to a MAD-equivalent derived from the mean absolute
        deviation (rescaled so the same 0.6745 constant applies).
      - if there's no spread at all (every bulk value is ~identical), the
        junction isn't testable -- returns the "no_variance" sentinel,
        same treatment as "low_n"/"error" downstream.
    """
    x = x[np.isfinite(x)]
    n = len(x)
    if n < n_threshold:
        return n, "low_n", "low_n"
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    if mad < tol:
        meanad = float(np.mean(np.abs(x - med)))
        if meanad < tol:
            return n, "no_variance", "no_variance"
        mad = meanad * (_MODZ_CONST / _MEANAD_CONST)
    return n, med, mad


def _fit_modz_rows(args):
    """Worker: fit a sub-block of rows. Returns list of result tuples."""
    mat_block, tol, n_threshold = args
    return [_fit_one_modz(mat_block[i], tol, n_threshold)
            for i in range(len(mat_block))]


def fit_modz_dist_chunk(mat, feat_names, tol, n_threshold, threads: int = 1):
    """Fit median/effective-MAD for every row in mat. Parallelised across
    rows using a process pool when threads > 1, same pattern as
    fit_beta_dist_chunk."""
    n = len(mat)
    if n == 0:
        return pd.DataFrame(columns=["n", "median", "mad"])

    n_workers = min(threads, n)
    if n_workers <= 1:
        results = [_fit_one_modz(mat[i], tol, n_threshold) for i in range(n)]
    else:
        chunk_size = max(1, (n + n_workers - 1) // n_workers)
        chunks = [
            (mat[i : i + chunk_size], tol, n_threshold)
            for i in range(0, n, chunk_size)
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
            results = []
            for block in ex.map(_fit_modz_rows, chunks):
                results.extend(block)

    return pd.DataFrame(results, index=feat_names, columns=["n", "median", "mad"])


# ---------------------------------------------------------------------------
# Per-metric beta pipeline
# ---------------------------------------------------------------------------

def _run_one_metric(
    combined_df:        pd.DataFrame,
    metric_col:         str,
    rescaled_col:       str,
    usage_col:          str,
    coverage_col:       str,
    id_col:             str,
    coverage_threshold: int,
    PSI_rescale_factor: float,
    n_threshold:        int,
    threads:            int,
    method:             str,
) -> pd.DataFrame:
    """Fits a per-junction reference distribution across bulk samples and
    scores every row (bulk + hap1/hap2) against it, using either:
      method="beta_binomial"    -- n/alpha/beta/expected/p1/p99/delta/p_value
      method="modified_zscore"  -- n/median/mad/modz
    columns per metric. See the module docstring / --bb-thresholds vs.
    --z-thresholds for how these feed into outlier identification.
    """
    is_ss_metric = metric_col in ("5ss_IR_ratio", "3ss_IR_ratio")

    if is_ss_metric:
        ss_pos_col = "5ss" if metric_col == "5ss_IR_ratio" else "3ss"
        combined_df = combined_df.copy()
        combined_df["_ss_id"] = combined_df["region"].str.split(":").str[0] + "_" + \
                                  combined_df[ss_pos_col].astype(str)
        fit_id_col = "_ss_id"
    else:
        fit_id_col = id_col

    bulk_sub = combined_df[combined_df["phasing"] == "bulk"][
        ["sample", fit_id_col, coverage_col, rescaled_col]
    ].copy()

    low_cov = bulk_sub[coverage_col].to_numpy(dtype=float) < coverage_threshold
    rvals   = np.array(pd.to_numeric(bulk_sub[rescaled_col], errors="coerce"), dtype=float)
    rvals[low_cov] = np.nan

    bulk_wide = (
        bulk_sub.assign(**{rescaled_col: rvals})
        .drop_duplicates(subset=["sample", fit_id_col])
        .pivot(index=fit_id_col, columns="sample", values=rescaled_col)
        .astype(np.float32)
    )

    n_col = f"n_{metric_col}"
    if method == "beta_binomial":
        alpha_col    = f"alpha_{metric_col}"
        beta_col     = f"beta_{metric_col}"
        expected_col = f"expected_{metric_col}"
        delta_col    = f"delta_{metric_col}"
        p_col        = f"p_value_{metric_col}"
        p1_col       = f"p1_{metric_col}"
        p99_col      = f"p99_{metric_col}"
        fit_indicator_col = alpha_col
        empty_cols = (alpha_col, beta_col, expected_col, p1_col, p99_col, delta_col, p_col)
    else:
        median_col = f"median_{metric_col}"
        mad_col    = f"mad_{metric_col}"
        modz_col   = f"modz_{metric_col}"
        fit_indicator_col = median_col
        empty_cols = (median_col, mad_col, modz_col)

    if bulk_wide.empty:
        combined_df[n_col] = 0
        for col in empty_cols:
            combined_df[col] = "low_n"
        return combined_df

    mat        = bulk_wide.to_numpy(dtype=np.float32)
    feat_names = bulk_wide.index.tolist()

    # Count non-NaN values per junction (row) — only fit rows with enough data
    n_valid  = np.sum(~np.isnan(mat), axis=1)

    fit_mask    = n_valid >= n_threshold
    low_n_names = [feat_names[i] for i in range(len(feat_names)) if not fit_mask[i]]
    fit_names   = [feat_names[i] for i in range(len(feat_names)) if fit_mask[i]]

    if method == "beta_binomial":
        if fit_names:
            beta_df = fit_beta_dist_chunk(
                mat[fit_mask], fit_names, PSI_rescale_factor, n_threshold, threads
            ).reset_index().rename(
                columns={"index": fit_id_col, "n": n_col, "alpha": alpha_col,
                         "beta_param": beta_col, "expected": expected_col})
        else:
            beta_df = pd.DataFrame(columns=[fit_id_col, n_col, alpha_col, beta_col, expected_col])

        if low_n_names:
            low_n_df = pd.DataFrame({
                fit_id_col:   low_n_names,
                n_col:        [int(n_valid[i]) for i in range(len(feat_names)) if not fit_mask[i]],
                alpha_col:    "low_n",
                beta_col:     "low_n",
                expected_col: "low_n",
            })
            beta_df = pd.concat([beta_df, low_n_df], ignore_index=True)

        combined_df = combined_df.merge(beta_df, on=fit_id_col, how="left")
        combined_df[n_col] = combined_df[n_col].fillna(0)
        for col in (alpha_col, beta_col, expected_col):
            combined_df[col] = combined_df[col].fillna("low_n")
    else:
        if fit_names:
            modz_df = fit_modz_dist_chunk(
                mat[fit_mask], fit_names, PSI_rescale_factor, n_threshold, threads
            ).reset_index().rename(
                columns={"index": fit_id_col, "n": n_col, "median": median_col, "mad": mad_col})
        else:
            modz_df = pd.DataFrame(columns=[fit_id_col, n_col, median_col, mad_col])

        if low_n_names:
            low_n_df = pd.DataFrame({
                fit_id_col: low_n_names,
                n_col:      [int(n_valid[i]) for i in range(len(feat_names)) if not fit_mask[i]],
                median_col: "low_n",
                mad_col:    "low_n",
            })
            modz_df = pd.concat([modz_df, low_n_df], ignore_index=True)

        combined_df = combined_df.merge(modz_df, on=fit_id_col, how="left")
        combined_df[n_col] = combined_df[n_col].fillna(0)
        for col in (median_col, mad_col):
            combined_df[col] = combined_df[col].fillna("low_n")

    fit_vals     = combined_df[fit_indicator_col]
    is_low_n     = fit_vals == "low_n"
    is_error     = fit_vals == "error"          # beta_binomial only
    is_no_var    = fit_vals == "no_variance"    # modified_zscore only
    has_fit      = ~is_low_n & ~is_error & ~is_no_var
    coverage_arr = combined_df[coverage_col].to_numpy(dtype=float)
    has_cov      = coverage_arr >= coverage_threshold

    # Identify low-phased hap rows upfront (before any testing)
    lp_col_name = f"low_phased_{metric_col}"
    is_hap = combined_df["phasing"].isin(["hap1", "hap2"])
    if lp_col_name in combined_df.columns:
        is_low_phased = combined_df[lp_col_name].astype(bool)
    else:
        is_low_phased = pd.Series(False, index=combined_df.index)
    # low_phased only applies when the row also has coverage (otherwise low_coverage wins)
    is_low_phased_hap = is_hap & is_low_phased & has_cov

    if method == "beta_binomial":
        # Use ceil(n * percentile) as index for conservative bounds:
        # e.g. n=107: ceil(107*0.01)=2 → take 3rd value (0-based index 2),
        # excluding the bottom 2; p99 takes the (n - ceil(n*0.01) - 1)th value.
        p1_vals  = np.full(len(feat_names), np.nan)
        p99_vals = np.full(len(feat_names), np.nan)
        for row_i in range(len(mat)):
            row  = mat[row_i]
            vals = np.sort(row[~np.isnan(row)])
            n    = len(vals)
            if n == 0:
                continue
            if n <= 10:
                p1_vals[row_i]  = vals[0]
                p99_vals[row_i] = vals[-1]
            else:
                k = math.ceil(n * 0.01)
                p1_vals[row_i]  = vals[k]
                p99_vals[row_i] = vals[n - 1 - k]
        p1_series  = pd.Series(p1_vals,  index=feat_names, name=p1_col)
        p99_series = pd.Series(p99_vals, index=feat_names, name=p99_col)
        combined_df = combined_df.merge(
            p1_series.reset_index().rename(columns={"index": fit_id_col}),
            on=fit_id_col, how="left")
        combined_df = combined_df.merge(
            p99_series.reset_index().rename(columns={"index": fit_id_col}),
            on=fit_id_col, how="left")
        for pc in (p1_col, p99_col):
            combined_df[pc] = combined_df[pc].where(combined_df[pc].notna(), other="low_n")

        rescaled_v = pd.to_numeric(combined_df[rescaled_col], errors="coerce")
        p1_v       = pd.to_numeric(combined_df[p1_col],       errors="coerce")
        p99_v      = pd.to_numeric(combined_df[p99_col],      errors="coerce")
        delta_num  = np.where(
            rescaled_v > p99_v, rescaled_v - p99_v,
            np.where(rescaled_v < p1_v, rescaled_v - p1_v, 0.0)
        )
        delta_v = pd.Series(delta_num, index=combined_df.index).astype(object)
        delta_v[is_low_n]                       = "low_n"
        delta_v[is_error]                       = "error"
        delta_v[has_fit & ~has_cov]             = "low_coverage"
        delta_v[has_fit & is_low_phased_hap]    = "low_phased_coverage"
        combined_df[delta_col] = delta_v

        # Assign p_value sentinel strings in priority order:
        #   1. low_coverage  — row coverage < threshold (takes priority over everything)
        #   2. low_phased_coverage — hap row that fails phasing check for this metric
        #   3. low_n / error — distribution could not be fit
        #   4. actual test   — all other rows with a fit and sufficient coverage
        combined_df[p_col] = pd.Series(np.nan, index=combined_df.index, dtype=object)  # default; filled below
        combined_df.loc[is_low_n,                    p_col] = "low_n"
        combined_df.loc[is_error,                    p_col] = "error"
        combined_df.loc[has_fit & ~has_cov,          p_col] = "low_coverage"
        combined_df.loc[has_fit & is_low_phased_hap, p_col] = "low_phased_coverage"

        # Testable: has a fitted distribution, sufficient coverage, and is NOT low-phased
        testable = combined_df[has_fit & has_cov & ~is_low_phased_hap].copy()
        if len(testable) > 0:
            orig_index = testable.index
            tested = beta_binomial_test_chunk(
                testable.reset_index(drop=True),
                coverage_threshold, metric_col, usage_col, coverage_col, p_col,
                threads,
            )
            combined_df.loc[orig_index, p_col] = tested[p_col].values
    else:
        # modZ sentinel priority mirrors p_value's above, applied to modz_col.
        combined_df[modz_col] = pd.Series(np.nan, index=combined_df.index, dtype=object)
        combined_df.loc[is_low_n,                    modz_col] = "low_n"
        combined_df.loc[is_no_var,                    modz_col] = "no_variance"
        combined_df.loc[has_fit & ~has_cov,          modz_col] = "low_coverage"
        combined_df.loc[has_fit & is_low_phased_hap, modz_col] = "low_phased_coverage"

        testable = has_fit & has_cov & ~is_low_phased_hap
        if testable.any():
            rescaled_v = pd.to_numeric(combined_df.loc[testable, rescaled_col], errors="coerce")
            median_v   = pd.to_numeric(combined_df.loc[testable, median_col],   errors="coerce")
            mad_v      = pd.to_numeric(combined_df.loc[testable, mad_col],      errors="coerce")
            modz_v     = _MODZ_CONST * (rescaled_v - median_v) / mad_v
            combined_df.loc[testable, modz_col] = modz_v.astype(object)

    if is_ss_metric:
        combined_df = combined_df.drop(columns=["_ss_id"])

    return combined_df


def run_all_metrics(
    combined_df:        pd.DataFrame,
    coverage_threshold: int,
    PSI_rescale_factor: float,
    n_threshold:        int,
    threads:            int,
    has_ipa:            bool,
    method:             str,
    no_ss_ir:           bool = False,
) -> pd.DataFrame:
    metrics = [
        ("junction_PSI_approx", "rescaled_junction_PSI_approx", "junction_usage", "junction_coverage_approx", "junction"),
        ("junction_PSI",        "rescaled_junction_PSI",        "junction_usage", "junction_coverage",        "junction"),
        ("junction_full_IR_ratio", "rescaled_junction_full_IR_ratio", "junction_full_IR_count", "junction_coverage", "junction"),
    ]
    if not no_ss_ir:
        metrics.insert(2, ("5ss_IR_ratio", "rescaled_5ss_IR_ratio", "5ss_usage", "5ss_coverage", "junction"))
        metrics.insert(3, ("3ss_IR_ratio", "rescaled_3ss_IR_ratio", "3ss_usage", "3ss_coverage", "junction"))
    if has_ipa:
        metrics.append(
            ("junction_IPA_ratio", "rescaled_junction_IPA_ratio", "junction_IPA_count", "5ss_coverage", "junction")
        )
    for metric_col, rescaled_col, usage_col, coverage_col, id_col in metrics:
        if metric_col not in combined_df.columns:
            continue
        combined_df = _run_one_metric(
            combined_df, metric_col, rescaled_col, usage_col, coverage_col,
            id_col, coverage_threshold, PSI_rescale_factor, n_threshold, threads,
            method,
        )
    return combined_df


def run_gene_stats_pipeline(
    gene:                str,
    combined:            pd.DataFrame,
    approx_only:         bool,
    coverage_threshold:  int,
    PSI_rescale_factor:  float,
    n_threshold:         int,
    threads:             int,
    has_ipa:             bool,
    no_ss_ir:            bool,
    method:              str,
) -> pd.DataFrame:
    """Fits a per-junction reference distribution across this gene's bulk
    samples (Beta or median/MAD, depending on `method`) and scores every
    sample/junction against it, for every applicable metric. This is the
    per-gene statistical-testing step; it only needs this one gene's own
    data (no cross-gene information), unlike the FDR correction (beta_binomial
    only) that follows once every gene has been tested."""
    print(f"  Gene: {gene}  Fitting + scoring ({method}) ...")
    t0 = time.time()

    if approx_only:
        combined = _run_one_metric(
            combined, "junction_PSI_approx", "rescaled_junction_PSI_approx",
            "junction_usage", "junction_coverage_approx", "junction",
            coverage_threshold, PSI_rescale_factor, n_threshold, threads,
            method,
        )
    else:
        combined = run_all_metrics(combined, coverage_threshold, PSI_rescale_factor,
                                   n_threshold, threads, has_ipa, method, no_ss_ir)

    # Diagnostic fit/test counts. "Fit" columns are alpha_ (beta_binomial) or
    # median_ (modified_zscore); "test" columns are p_value_ or modz_.
    fit_prefix  = "alpha_" if method == "beta_binomial" else "median_"
    test_prefix = "p_value_" if method == "beta_binomial" else "modz_"
    not_fittable = ("low_n", "error") if method == "beta_binomial" else ("low_n", "no_variance")

    test_cols_present = [c for c in combined.columns if c.startswith(test_prefix)]
    n_fit = 0; n_tests = 0
    bulk_combined = combined[combined["phasing"] == "bulk"]
    for test_col in test_cols_present:
        mc    = test_col.replace(test_prefix, "")
        a_col = f"{fit_prefix}{mc}"
        if a_col not in combined.columns: continue
        is_ss = any(x in test_col for x in ("5ss_IR_ratio", "3ss_IR_ratio"))
        if is_ss:
            ss_col = "5ss" if "5ss" in test_col else "3ss"
            if ss_col in bulk_combined.columns:
                n_fit += int(
                    bulk_combined[bulk_combined[a_col].apply(
                        lambda x: x not in not_fittable
                    )][ss_col].nunique()
                )
        else:
            n_fit += int(
                bulk_combined[bulk_combined[a_col].apply(
                    lambda x: x not in not_fittable
                )]["junction"].nunique()
            )
        n_tests += int(pd.to_numeric(combined[test_col], errors="coerce").notna().sum())
    print(f"       → fit {n_fit:,} distributions and scored {n_tests:,} rows "
          f"({time.time()-t0:.2f}s)")
    return combined


# ---------------------------------------------------------------------------
# Output column definitions
# ---------------------------------------------------------------------------

# NOTE: junction_type inserted after 3ss — silently skipped if absent (no --gtf)
def _per_metric_output_cols(metric_col: str) -> List[str]:
    """Every column _run_one_metric can produce for one metric, across both
    methods. Only whichever set was actually computed for a given run will
    be present in the DataFrame -- select_output_columns() silently skips
    the rest."""
    return [
        f"n_{metric_col}",
        # beta_binomial columns
        f"alpha_{metric_col}", f"beta_{metric_col}", f"expected_{metric_col}",
        f"p1_{metric_col}", f"p99_{metric_col}",
        f"delta_{metric_col}", f"p_value_{metric_col}", f"padj_{metric_col}",
        # modified_zscore columns
        f"median_{metric_col}", f"mad_{metric_col}", f"modz_{metric_col}",
    ]


_OUTPUT_COLS = [
    "sample", "gene", "gene_rank", "region", "phasing", "junction", "5ss", "3ss",
    "junction_type",   # only present when --gtf provided
    "junction_usage",
    "junction_read_diversity",
    "junction_coverage_approx",
    "junction_PSI_approx", "rescaled_junction_PSI_approx",
    *_per_metric_output_cols("junction_PSI_approx"),
    "junction_coverage",
    "junction_PSI", "rescaled_junction_PSI",
    *_per_metric_output_cols("junction_PSI"),
    "5ss_usage", "5ss_coverage",
    "5ss_IR_ratio", "rescaled_5ss_IR_ratio",
    *_per_metric_output_cols("5ss_IR_ratio"),
    "3ss_usage", "3ss_coverage",
    "3ss_IR_ratio", "rescaled_3ss_IR_ratio",
    *_per_metric_output_cols("3ss_IR_ratio"),
    "junction_full_IR_count",   # singular
    "junction_full_IR_ratio", "rescaled_junction_full_IR_ratio",
    *_per_metric_output_cols("junction_full_IR_ratio"),
    "junction_IPA_count",       # singular
    "junction_IPA_ratio", "rescaled_junction_IPA_ratio",
    *_per_metric_output_cols("junction_IPA_ratio"),
]


def select_output_columns(df: pd.DataFrame, *_) -> pd.DataFrame:
    """Select and reorder to canonical output schema. Missing columns silently skipped."""
    present = [c for c in _OUTPUT_COLS if c in df.columns]
    return df[present]


# ---------------------------------------------------------------------------
# GTF parsing — junctions for QC figures and junction_type classification
# ---------------------------------------------------------------------------

def parse_gtf_junctions(
    gtf_path: str,
    gene_names: List[str],
) -> Dict[str, Dict]:
    """
    Returns {gene: {"canonical_junctions": set, "all_junctions": set}}
    Junction strings are "chrom_ss1_ss2".
    """
    import gzip as _gz
    import re as _re

    gene_set = set(gene_names)

    def _attr(attr_str, key):
        m = _re.search(rf'{key}' + r'\s+"([^"]+)"', attr_str)
        return m.group(1) if m else None

    def _strip_ver(s):
        return _re.sub(r'\.\d+$', '', s) if s else s

    tx_canonical: set = set()
    gene_exons: Dict[str, Dict[str, List[Tuple[str, int, int]]]] = defaultdict(lambda: defaultdict(list))

    open_fn = _gz.open if gtf_path.endswith(".gz") else open

    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            feat = parts[2]
            if feat not in ("transcript", "exon"): continue

            chrom = parts[0]
            start = int(parts[3]) - 1   # 0-based
            end   = int(parts[4])        # half-open
            attrs = parts[8]

            gname = _attr(attrs, "gene_name") or _attr(attrs, "gene_symbol")
            gid   = _strip_ver(_attr(attrs, "gene_id"))
            gene_key = None
            if gname and gname in gene_set:
                gene_key = gname
            elif gid and gid in gene_set:
                gene_key = gid
            if gene_key is None:
                continue

            tx_id = _attr(attrs, "transcript_id")
            if tx_id is None: continue

            if feat == "transcript":
                if "Ensembl_canonical" in attrs:
                    tx_canonical.add(tx_id)
            elif feat == "exon":
                gene_exons[gene_key][tx_id].append((chrom, start, end))

    result: Dict[str, Dict] = {}
    for gene_key, tx_dict in gene_exons.items():
        canonical_jxns: set = set()
        all_jxns_set:   set = set()
        for tx_id, exons in tx_dict.items():
            if len(exons) < 2: continue
            exons_sorted = sorted(exons, key=lambda x: x[1])
            for i in range(len(exons_sorted) - 1):
                chrom_e, _, end_e   = exons_sorted[i]
                chrom_n, start_n, _ = exons_sorted[i + 1]
                if chrom_e != chrom_n: continue
                jxn = f"{chrom_e}_{end_e + 1}_{start_n}"
                all_jxns_set.add(jxn)
                if tx_id in tx_canonical:
                    canonical_jxns.add(jxn)
        result[gene_key] = {
            "canonical_junctions": canonical_jxns,
            "all_junctions":       all_jxns_set,
        }
    return result


def assign_junction_types(
    df: pd.DataFrame,
    gtf_junctions: Dict[str, Dict],
) -> pd.DataFrame:
    """
    Add a 'junction_type' column with values 'canonical', 'annotated', or 'novel'.
    - canonical : junction is in the Ensembl canonical transcript for that gene
    - annotated : junction is in any transcript for that gene (but not canonical)
    - novel     : junction is not in any annotated transcript

    Operates per-gene. Returns df with junction_type column added.
    """
    df = df.copy()
    jtype = pd.Series("novel", index=df.index, dtype=object)

    for gene, gdf in df.groupby("gene"):
        if gene not in gtf_junctions:
            continue  # leave as "novel"
        canonical_set = gtf_junctions[gene]["canonical_junctions"]
        annotated_set = gtf_junctions[gene]["all_junctions"]
        jxns = df.loc[gdf.index, "junction"]
        # annotated (not canonical)
        is_annotated = jxns.isin(annotated_set) & ~jxns.isin(canonical_set)
        # canonical
        is_canonical = jxns.isin(canonical_set)
        jtype.loc[gdf.index[is_annotated]] = "annotated"
        jtype.loc[gdf.index[is_canonical]] = "canonical"

    df["junction_type"] = jtype
    return df


# ---------------------------------------------------------------------------
# Event classification
# ---------------------------------------------------------------------------

def _classify_events_for_metric(
    metric_col:         str,
    outlier_set:        set,
    grp_coords:         Dict[Tuple, Dict[str, Tuple[int, int]]],
    strand_map:         Dict[str, str],
    padj_threshold:     float,
    delta_threshold:    float,
    positive_delta_set: set,
    delta_map:          Dict[Tuple, float],
) -> Dict[Tuple, List[str]]:
    """
    Classify splicing events for one metric.

    delta_map: {(sample, gene, phasing, jxn): delta_value} for all rows with
    a numeric delta.  Used to determine direction for alt_5ss/alt_3ss/exon
    classification.

    positive_delta_set: subset of outlier_set where delta >= +delta_threshold.
    For IR/IPA events (5ss_IR, 3ss_IR, full_IR, IPA) the row must be in
    positive_delta_set to receive the event label.

    Event rules:
      alt_3ss: query junction J is an outlier; there exists another junction J'
               in the same (sample, gene, phasing) group sharing J's 5ss but
               with a different 3ss, and delta(J') has the opposite sign to
               delta(J).  J' need not be an outlier itself.
      alt_5ss: symmetric — same 3ss, different 5ss, opposite-sign delta.
      exon_skipping / exon_inclusion:
               Three outlier junctions form a cassette triple:
                 long:  ss1a → ss2b
                 left:  ss1a → ss2a  (ss2a < ss2b)
                 right: ss1b → ss2b  (ss1b > ss1a)
               with ss2a == ss1b (the skipped exon boundaries meet).
               If long goes up (delta > 0) and both short go down → exon_skipping
               If long goes down (delta < 0) and both short go up  → exon_inclusion
    """
    allowed_events = _METRIC_EVENTS.get(metric_col, [])
    if not allowed_events:
        return {}

    outlier_by_grp: Dict[Tuple, set] = defaultdict(set)
    for key in outlier_set:
        sample, gene, phasing, jxn = key
        outlier_by_grp[(sample, gene, phasing)].add(jxn)

    event_map: Dict[Tuple, List[str]] = defaultdict(list)

    def _add(sample, gene, phasing, jxn, ev):
        key = (sample, gene, phasing, jxn)
        if ev not in event_map[key]:
            event_map[key].append(ev)

    def _delta(sample, gene, phasing, jxn):
        return delta_map.get((sample, gene, phasing, jxn), float("nan"))

    def _sign(d):
        """Return 1, -1, or 0."""
        if d > 0: return 1
        if d < 0: return -1
        return 0

    for grp_key, outlier_jxns in outlier_by_grp.items():
        sample, gene, phasing = grp_key
        strand = strand_map.get(gene, "+")
        coords = grp_coords[grp_key]  # ALL junctions for this group

        for jxn in outlier_jxns:
            if jxn not in coords:
                continue
            ss1, ss2 = coords[jxn]
            five_ss  = ss1 if strand == "+" else ss2
            three_ss = ss2 if strand == "+" else ss1

            row_key  = (sample, gene, phasing, jxn)
            is_pos_delta = row_key in positive_delta_set

            # IR/IPA events: only assign when delta is positive
            if "5ss_IR"  in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "5ss_IR")
            if "3ss_IR"  in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "3ss_IR")
            if "full_IR" in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "full_IR")
            if "IPA"     in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "IPA")

            if "alt_5ss" in allowed_events or "alt_3ss" in allowed_events or \
               "alt_5ss_approx" in allowed_events or "alt_3ss_approx" in allowed_events:
                d_jxn = _delta(sample, gene, phasing, jxn)
                s_jxn = _sign(d_jxn)
                alt_5ss_label = next((e for e in allowed_events if e.startswith("alt_5ss")), None)
                alt_3ss_label = next((e for e in allowed_events if e.startswith("alt_3ss")), None)
                # Partner must also be an outlier
                for partner_jxn in outlier_jxns:
                    if partner_jxn == jxn or partner_jxn not in coords:
                        continue
                    p_ss1, p_ss2 = coords[partner_jxn]
                    p_five  = p_ss1 if strand == "+" else p_ss2
                    p_three = p_ss2 if strand == "+" else p_ss1
                    d_partner = _delta(sample, gene, phasing, partner_jxn)
                    s_partner = _sign(d_partner)
                    opposite = (s_jxn != 0 and s_partner != 0 and s_jxn != s_partner)

                    if alt_3ss_label:
                        if p_five == five_ss and p_three != three_ss and opposite:
                            _add(sample, gene, phasing, jxn,         alt_3ss_label)
                            _add(sample, gene, phasing, partner_jxn, alt_3ss_label)

                    if alt_5ss_label:
                        if p_three == three_ss and p_five != five_ss and opposite:
                            _add(sample, gene, phasing, jxn,         alt_5ss_label)
                            _add(sample, gene, phasing, partner_jxn, alt_5ss_label)

            skip_label    = next((e for e in allowed_events if e.startswith("exon_skipping")),  None)
            incl_label    = next((e for e in allowed_events if e.startswith("exon_inclusion")), None)
            if skip_label or incl_label:
                d_long = _delta(sample, gene, phasing, jxn)
                s_long = _sign(d_long)
                if s_long == 0:
                    continue
                # jxn is the long junction (ss1a → ss2b).
                # Find outlier short junctions sharing ss1 (left) and ss2 (right).
                # No requirement that left's ss2 == right's ss1 — the skipped exon
                # can have any size; all that matters is the outer coordinates match.
                left_candidates = [
                    (j, s1, s2) for j, (s1, s2) in coords.items()
                    if j in outlier_jxns and j != jxn
                    and s1 == ss1 and s2 < ss2
                ]
                right_candidates = [
                    (j, s1, s2) for j, (s1, s2) in coords.items()
                    if j in outlier_jxns and j != jxn
                    and s2 == ss2 and s1 > ss1
                ]
                for jl, ls1, ls2 in left_candidates:
                    for jr, rs1, rs2 in right_candidates:
                        d_left  = _delta(sample, gene, phasing, jl)
                        d_right = _delta(sample, gene, phasing, jr)
                        s_left  = _sign(d_left)
                        s_right = _sign(d_right)
                        # Both short junctions must go opposite to the long
                        if s_left == 0 or s_right == 0:
                            continue
                        if s_left == s_long or s_right == s_long:
                            continue
                        if skip_label and s_long > 0:
                            for j in (jxn, jl, jr):
                                _add(sample, gene, phasing, j, skip_label)
                        if incl_label and s_long < 0:
                            for j in (jxn, jl, jr):
                                _add(sample, gene, phasing, j, incl_label)

    return dict(event_map)


def classify_all_events(
    sig_df:              pd.DataFrame,
    final_df:            pd.DataFrame,
    strand_map:          Dict[str, str],
    effect_threshold:    float,
    threads:             int,
    has_ipa:             bool,
    computed_metrics:    List[str],
    effect_col_fn,
    unreliable_hap_outliers: Dict[str, set],
) -> pd.DataFrame:
    """
    effect_col_fn(mc) -> the column name holding this metric's signed
    "effect" value per row -- delta_{mc} for beta_binomial, modz_{mc} for
    modified_zscore. Positive values in that column are what count as
    "more of the event" for IR/IPA metrics and for alt-5ss/3ss/exon
    direction comparisons; effect_threshold is the corresponding magnitude
    cutoff (delta_threshold or a z-threshold).
    """
    if sig_df.empty:
        sig_df = sig_df.copy()
        sig_df["event_type"] = ""
        return sig_df

    active_metrics = [mc for mc in computed_metrics if _METRIC_EVENTS.get(mc)]
    if not active_metrics:
        sig_df = sig_df.copy(); sig_df["event_type"] = "other"; return sig_df

    # Split grp_coords by gene to keep per-task payloads small
    all_genes  = sorted(final_df["gene"].unique())
    n_genes    = len(all_genes)
    batch_size = max(1, n_genes // (threads * 4)) if threads > 1 else n_genes
    gene_batches = [all_genes[i:i + batch_size]
                    for i in range(0, n_genes, batch_size)]

    # grp_coords_by_gene: {gene: {(sample, gene, phasing): {jxn: (ss1, ss2)}}}
    grp_coords_by_gene: Dict[str, Dict] = {}
    for g in all_genes:
        grp_coords_by_gene[g] = {}
    for _, row in final_df[["sample", "gene", "phasing", "junction"]].iterrows():
        g   = row["gene"]
        key = (row["sample"], g, row["phasing"])
        jxn = row["junction"]
        parts = jxn.split("_")
        grp_coords_by_gene[g][key] = grp_coords_by_gene[g].get(key, {})
        grp_coords_by_gene[g][key][jxn] = (int(parts[-2]), int(parts[-1]))

    # Per-metric sets / maps (full, then sliced per batch at dispatch time)
    metric_outlier_sets:   Dict[str, set] = {}
    metric_pos_delta_sets: Dict[str, set] = {}
    metric_delta_maps:     Dict[str, Dict[Tuple, float]] = {}

    _ir_ipa_metrics = frozenset((
        "5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio", "junction_IPA_ratio"
    ))
    _has_jxn_type = "junction_type" in sig_df.columns

    for mc in computed_metrics:
        ocol    = f"outlier_{mc}"
        delta_c = effect_col_fn(mc)
        if ocol not in sig_df.columns:
            metric_outlier_sets[mc]   = set()
            metric_pos_delta_sets[mc] = set()
        else:
            mask = sig_df[ocol].astype(bool)

            # For IR/IPA metrics: additionally require positive delta and,
            # if --gtf was provided, canonical or annotated junction_type.
            if mc in _ir_ipa_metrics:
                if delta_c in sig_df.columns:
                    dv   = pd.to_numeric(sig_df[delta_c], errors="coerce")
                    mask = mask & dv.ge(effect_threshold)
                if _has_jxn_type:
                    mask = mask & sig_df["junction_type"].isin(["canonical", "annotated"])

            metric_outlier_sets[mc] = set(zip(
                sig_df.loc[mask, "sample"], sig_df.loc[mask, "gene"],
                sig_df.loc[mask, "phasing"], sig_df.loc[mask, "junction"],
            ))
            if delta_c in sig_df.columns:
                dv       = pd.to_numeric(sig_df[delta_c], errors="coerce")
                pos_mask = mask & dv.ge(effect_threshold)
            else:
                pos_mask = mask
            metric_pos_delta_sets[mc] = set(zip(
                sig_df.loc[pos_mask, "sample"], sig_df.loc[pos_mask, "gene"],
                sig_df.loc[pos_mask, "phasing"], sig_df.loc[pos_mask, "junction"],
            ))
        if delta_c not in final_df.columns:
            metric_delta_maps[mc] = {}
        else:
            dv    = pd.to_numeric(final_df[delta_c], errors="coerce")
            valid = dv.notna()
            metric_delta_maps[mc] = dict(zip(
                zip(final_df.loc[valid, "sample"], final_df.loc[valid, "gene"],
                    final_df.loc[valid, "phasing"], final_df.loc[valid, "junction"]),
                dv[valid],
            ))

    def _slice_batch(gene_batch, mc):
        gene_set = set(gene_batch)
        gc = {}
        for g in gene_batch:
            gc.update(grp_coords_by_gene.get(g, {}))
        os_ = {k for k in metric_outlier_sets.get(mc, set())   if k[1] in gene_set}
        ps_ = {k for k in metric_pos_delta_sets.get(mc, set()) if k[1] in gene_set}
        dm_ = {k: v for k, v in metric_delta_maps.get(mc, {}).items() if k[1] in gene_set}
        return gc, os_, ps_, dm_

    all_event_maps: List[Dict[Tuple, List[str]]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futs = {}
        for mc in active_metrics:
            for batch in gene_batches:
                gc, os_, ps_, dm_ = _slice_batch(batch, mc)
                if not os_:
                    continue
                fut = ex.submit(
                    _classify_events_for_metric, mc,
                    os_, gc, strand_map, effect_threshold, effect_threshold, ps_, dm_,
                )
                futs[fut] = (mc, batch[0])
        for fut in concurrent.futures.as_completed(futs):
            mc, g0 = futs[fut]
            try:
                all_event_maps.append(fut.result())
            except Exception as e:
                print(f"[WARNING] Classification error for {mc} (batch @{g0}): {e}")
                traceback.print_exc()

    merged: Dict[Tuple, set] = defaultdict(set)
    for event_map in all_event_maps:
        for key, events in event_map.items():
            merged[key].update(events)

    sig_df = sig_df.copy()
    sig_records = sig_df[["sample", "gene", "phasing", "junction"]].to_dict("records")
    event_strs = []
    for rec in sig_records:
        key = (rec["sample"], rec["gene"], rec["phasing"], rec["junction"])
        evs = merged.get(key, set())
        event_strs.append(",".join(sorted(evs)) if evs else "none")
    sig_df["event_type"] = event_strs
    return sig_df


# ---------------------------------------------------------------------------
# QC figures
# ---------------------------------------------------------------------------

# QC figure definitions: each entry is
#   (coverage_col, file_suffix, [bar_metrics])
# The file will be: {prefix}_qc_{file_suffix}_{subset}.pdf
_QC_FIGURES = [
    ("junction_coverage_approx", "junction_coverage_approx", ["junction_PSI_approx"]),
    ("junction_coverage",        "junction_coverage",        ["junction_PSI", "junction_full_IR_ratio"]),
    ("5ss_coverage",             "5ss_coverage",             ["5ss_IR_ratio", "junction_IPA_ratio"]),
    ("3ss_coverage",             "3ss_coverage",             ["3ss_IR_ratio"]),
]
# Coverage columns that belong to splice-site (rather than junction) metrics
_SS_COVERAGE_COLS = frozenset(("5ss_coverage", "3ss_coverage"))


def make_qc_figure(
    final_df: pd.DataFrame,
    gtf_junctions: Dict[str, Dict],
    cov_col: str,
    companions: List[str],
    file_suffix: str,
    jxn_subset: str,
    coverage_threshold: int,
    out_pdf: str,
    has_ipa: bool,
    fit_prefix: str = "alpha_",
    not_fittable: Tuple[str, ...] = ("low_n", "error"),
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        print("[WARNING] matplotlib not available; skipping QC figures.")
        return

    color = "#d95d5b" if jxn_subset == "canonical" else "#4c8fca"
    cmap  = mcolors.LinearSegmentedColormap.from_list("c", ["white", color])

    companions = [m for m in companions
                  if has_ipa or m != "junction_IPA_ratio"]

    df = final_df  # already bulk-only, pre-sliced by caller
    genes = [g for g in gtf_junctions if g in df["gene"].unique()]
    if not genes:
        print(f"[WARNING] No overlapping genes for {file_suffix} ({jxn_subset}). Skipping.")
        return
    samples = sorted(df["sample"].unique())
    jxn_key = "canonical_junctions" if jxn_subset == "canonical" else "all_junctions"
    genes = [g for g in genes if len(gtf_junctions[g][jxn_key]) > 0]
    if not genes:
        print(f"[WARNING] No genes with {jxn_subset} junctions for {file_suffix}. Skipping.")
        return

    n_g = len(genes); n_s = len(samples)
    cov_mat = np.full((n_g, n_s), np.nan)
    cnt_mat = np.full((n_g, n_s), np.nan)
    fit_mat = np.full((n_g, len(companions)), np.nan)
    jxn_n   = []
    sample_idx = {s: i for i, s in enumerate(samples)}
    gene_idx   = {g: i for i, g in enumerate(genes)}

    for gene in genes:
        gi      = gene_idx[gene]
        jxn_set = gtf_junctions[gene][jxn_key]
        jxn_n.append(len(jxn_set))
        gdf     = df[df["gene"] == gene]
        for ci, cm in enumerate(companions):
            acol = f"{fit_prefix}{cm}"
            if acol not in gdf.columns:
                fit_mat[gi, ci] = 0; continue
            sub = gdf[gdf["junction"].isin(jxn_set)].drop_duplicates("junction")
            if sub.empty:
                fit_mat[gi, ci] = 0
            else:
                n_fitted = sub[acol].apply(lambda x: x not in not_fittable and x is not None).sum()
                fit_mat[gi, ci] = float(n_fitted) / len(jxn_set) if jxn_set else np.nan
        for sample in samples:
            si = sample_idx[sample]
            sub = gdf[(gdf["sample"] == sample) & gdf["junction"].isin(jxn_set)]
            if not jxn_set or sub.empty:
                cov_mat[gi, si] = 0.0; cnt_mat[gi, si] = 0.0; continue
            cov = pd.to_numeric(sub[cov_col], errors="coerce")
            n   = int((cov >= coverage_threshold).sum())
            cov_mat[gi, si] = float(n) / len(jxn_set)
            cnt_mat[gi, si] = float(n)

    gene_order   = np.lexsort((-np.nanmean(cnt_mat, axis=1),  -np.nanmean(cov_mat, axis=1)))
    sample_order = np.lexsort((-np.nanmean(cnt_mat, axis=0),  -np.nanmean(cov_mat, axis=0)))
    cov_mat_s = cov_mat[np.ix_(gene_order, sample_order)]
    fit_mat_s = fit_mat[gene_order]
    genes_s   = [genes[i]  for i in gene_order]
    jxn_n_s   = [jxn_n[i]  for i in gene_order]
    ylabels   = [f"{g}  (n={n})" for g, n in zip(genes_s, jxn_n_s)]

    heat_w = min(max(2.0, n_s * 0.05), 10.0)
    bar_w  = 1.5
    n_bars = len(companions)
    fig_w  = heat_w + bar_w * n_bars + 0.8
    fig_h  = max(2.0, n_g * 0.18 + 1.2)
    width_ratios = [heat_w] + [bar_w] * n_bars
    fig, axes = plt.subplots(1, 1 + n_bars, figsize=(fig_w, fig_h),
                              gridspec_kw={"width_ratios": width_ratios, "wspace": 0.05})
    if 1 + n_bars == 1: axes = [axes]

    ax0 = axes[0]
    feat_word = "splice sites" if cov_col in _SS_COVERAGE_COLS else "junctions"
    im = ax0.imshow(cov_mat_s, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ax0.set_title(f"Proportion of {feat_word} with {cov_col} ≥ {coverage_threshold}", fontsize=8, pad=3)
    ax0.set_xticks([]); ax0.set_yticks(range(n_g))
    ax0.set_yticklabels(ylabels, fontsize=6)
    fig.colorbar(im, ax=ax0, fraction=0.03, pad=0.01)

    y_pos = np.arange(n_g)
    for ci, cm in enumerate(companions):
        ax = axes[ci + 1]
        ax.barh(y_pos, fit_mat_s[:, ci], 0.6, color=color, alpha=0.85)
        ax.set_xlim(0, 1)
        ax.set_title(f"Proportion of {feat_word} with\n{cm} modeled", fontsize=7, pad=3)
        ax.set_yticks(y_pos); ax.set_yticklabels([])
        ax.tick_params(left=False); ax.spines["left"].set_visible(False)
        ax.invert_yaxis()
    for ax in axes:
        ax.set_ylim(n_g - 0.5, -0.5)

    fig.suptitle(f"{cov_col} — {jxn_subset.capitalize()} junctions", fontsize=9, y=1.01)
    t_qc_fig = time.time()
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"  QC figure → {out_pdf} ({time.time()-t_qc_fig:.2f}s)")


# ---------------------------------------------------------------------------
# Outlier heatmap + box plots
# ---------------------------------------------------------------------------

def make_outlier_heatmap(
    sig_df: pd.DataFrame,
    metric_col: str,
    effect_col: str,
    stat_label: str,
    threshold_desc: str,
    out_pdf: str,
) -> None:
    """sig_df must already be filtered to outlier rows for this metric
    (the caller applies the padj/delta or |modZ| threshold before calling).
    effect_col holds the signed effect value (delta_{metric} for
    beta_binomial, modz_{metric} for modified_zscore); stat_label/
    threshold_desc are just display text ("delta"/"modZ" and the threshold
    description shown in the title)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        return

    effect_v = pd.to_numeric(sig_df[effect_col], errors="coerce")
    df = sig_df.copy()
    df["_abs_delta"] = effect_v.abs()
    if df.empty: return

    agg = (df.groupby(["gene", "sample"])["_abs_delta"]
             .max().reset_index().rename(columns={"_abs_delta": "max_delta"}))
    genes   = sorted(agg["gene"].unique())
    samples = sorted(agg["sample"].unique())
    mat = np.full((len(genes), len(samples)), np.nan)
    g_idx = {g: i for i, g in enumerate(genes)}
    s_idx = {s: i for i, s in enumerate(samples)}
    for _, row in agg.iterrows():
        mat[g_idx[row["gene"]], s_idx[row["sample"]]] = row["max_delta"]

    gene_order   = np.lexsort((-np.nanmean(mat, axis=1), -(~np.isnan(mat)).sum(axis=1)))
    sample_order = np.lexsort((-np.nanmean(mat, axis=0), -(~np.isnan(mat)).sum(axis=0)))
    mat_s   = mat[np.ix_(gene_order, sample_order)]
    genes_s = [genes[i] for i in gene_order]
    n_g = len(genes_s); n_s = len(samples)

    cmap = mcolors.LinearSegmentedColormap.from_list("wd", ["#ffffff", "#912321"])
    cmap.set_bad(color="#f2f3f4")
    heat_w = min(max(2.0, n_s * 0.05), 12.0)
    fig_h  = max(2.0, n_g * 0.18 + 1.2)
    fig, ax = plt.subplots(figsize=(heat_w, fig_h))
    im = ax.imshow(mat_s, aspect="auto", cmap=cmap,
                   vmin=0, vmax=np.nanmax(mat_s), interpolation="nearest")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01, label=f"|{stat_label} {metric_col}|")
    ax.set_yticks(range(n_g)); ax.set_yticklabels(genes_s, fontsize=6)
    ax.set_xticks([])
    ax.set_title(f"Outlier heatmap: {metric_col}  ({threshold_desc})",
                 fontsize=9)
    fig.tight_layout()
    t_heat = time.time()
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"  Outlier heatmap → {out_pdf} ({time.time()-t_heat:.2f}s)")


def make_hit_boxplots(
    fmt_updated: pd.DataFrame,
    outlier_map: dict,
    metric_col: str,
    rescaled_col: str,
    outdir: str,
    prefix_name: str,
    tmp_dir: str,
    threshold_desc: str,
    stat_label: str,
    effect_threshold: float,
    is_ir_ipa: bool = False,
    has_jxn_type_filter: bool = False,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        return

    gene_jxn_map = outlier_map.get(metric_col, {})
    if not gene_jxn_map: return

    ind_dir = os.path.join(tmp_dir, metric_col)
    os.makedirs(ind_dir, exist_ok=True)
    out_pdf = os.path.join(outdir, f"{prefix_name}_boxplots_{metric_col}.pdf")
    n_hits = sum(len(jxns) for jxns in gene_jxn_map.values())
    n_written = 0

    t_box = time.time()
    with PdfPages(out_pdf) as pdf:
        # --- Title page ---
        fig_t, ax_t = plt.subplots(figsize=(4.5, 4.0))
        ax_t.axis("off")
        if is_ir_ipa:
            effect_str = f">= {effect_threshold}"
        else:
            effect_str = f"|{effect_threshold}|"
        jt_line = "\njunction_type: canonical, annotated" if has_jxn_type_filter else ""
        title_text = (
            f"Metric: {metric_col}\n"
            f"{threshold_desc}\n"
            f"{stat_label} threshold: {effect_str}"
            f"{jt_line}"
        )
        ax_t.text(0.5, 0.5, title_text, transform=ax_t.transAxes,
                  fontsize=12, va="center", ha="center",
                  bbox=dict(boxstyle="round,pad=0.6", facecolor="#ffffff", edgecolor="#0068a9"))
        fig_t.tight_layout()
        pdf.savefig(fig_t, bbox_inches="tight")
        plt.close(fig_t)

        for gene, jxn_map in sorted(gene_jxn_map.items()):
            for jxn, hit_samples in sorted(jxn_map.items()):
                try:
                    sub = fmt_updated[(fmt_updated["gene"] == gene) & (fmt_updated["junction"] == jxn)]
                    if sub.empty: continue
                    vals = pd.to_numeric(sub[rescaled_col], errors="coerce")
                    sub = sub.assign(_val=vals).dropna(subset=["_val"])
                    if sub.empty: continue

                    fig, ax = plt.subplots(figsize=(4.0, 3.2))
                    phasings = [p for p in ("bulk", "hap1", "hap2") if p in sub["phasing"].unique()]
                    data = [sub[sub["phasing"] == p]["_val"].to_numpy() for p in phasings]
                    bp = ax.boxplot(data, labels=phasings, showfliers=False, widths=0.5)
                    for i, p in enumerate(phasings):
                        p_vals = sub[sub["phasing"] == p]
                        is_hit = p_vals["sample"].isin(hit_samples)
                        x_jitter = np.random.normal(i + 1, 0.04, size=len(p_vals))
                        colors = ["#c0392b" if h else "#7f8c8d" for h in is_hit]
                        ax.scatter(x_jitter, p_vals["_val"], c=colors, s=14, alpha=0.8, zorder=3)
                    ax.set_title(f"{gene}: {jxn}", fontsize=9)
                    ax.set_ylabel(rescaled_col, fontsize=8)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)
                    n_written += 1
                except Exception as e:
                    # One bad hit shouldn't truncate the whole PDF -- log which
                    # (gene, junction) failed and why, close any half-built
                    # figure so it doesn't leak, and move on to the next hit.
                    print(f"[WARNING] Box plot failed for {gene}: {jxn} ({metric_col}): {e}")
                    traceback.print_exc()
                    try:
                        plt.close(fig)
                    except Exception:
                        pass
    print(f"  Box plots ({metric_col}, {n_written}/{n_hits}) → {out_pdf} ({time.time()-t_box:.2f}s)")


# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

def main() -> None:
    print("\n" + "*"*80)
    print("  Cohort Junction Outlier Detection")
    print("*"*80 + "\n")

    args        = parse_args()
    approx_only = args.approx

    # ---- Method + threshold specs ----
    # method determines both which statistic gets computed (beta-binomial
    # fit+test, or median/MAD + modified z-score) and how outliers are
    # identified from it.
    if args.bb_thresholds is not None:
        method = "beta_binomial"
        threshold_specs: List = []
        for tok in args.bb_thresholds:
            try:
                p_str, d_str = tok.split(":")
                threshold_specs.append((float(p_str), float(d_str)))
            except Exception:
                raise ValueError(f"Invalid threshold format '{tok}'. Expected padj:delta, e.g. 0.01:0.1")
        if not threshold_specs:
            print("[INFO] --bb-thresholds given with no values; will compute beta-binomial "
                  "statistics and write per-gene TSVs / QC figures only (no outlier identification).")
    else:
        method = "modified_zscore"
        z_vals = args.z_thresholds if args.z_thresholds else [3.5]
        threshold_specs = [float(z) for z in z_vals]

    print(f"Method: {method}")

    prefix        = args.outprefix.rstrip("/")
    outdir        = os.path.dirname(os.path.abspath(prefix))
    prefix_name   = os.path.basename(prefix)
    # Per-gene results columns differ by method (alpha/beta/expected/... for
    # beta_binomial vs. median/mad/modz for modified_zscore), so the method
    # is baked into this directory's name -- otherwise, re-running against
    # the same _raw/ input with the other method would silently overwrite
    # results computed under the first method, even though _raw/ itself
    # (written once by cohort_junction_analysis.py) supports being re-scored
    # under either method as many times as you like.
    base_results  = os.path.join(outdir, f"{prefix_name}_results_{method}")
    # Same reasoning as base_results above -- which columns a QC figure's
    # "proportion fit" bar reads (alpha_ vs. median_) depends on method.
    qc_dir        = os.path.join(outdir, f"{prefix_name}_qc_{method}")
    tmp_dir       = os.path.join(outdir, f"{prefix_name}_tmp")

    def _results_dir():
        os.makedirs(base_results, exist_ok=True); return base_results

    def _threshold_dir(subdir_name: str) -> str:
        d = os.path.join(outdir, f"{prefix_name}_{subdir_name}")
        os.makedirs(d, exist_ok=True)
        return d

    def _thr_subdir_name(thr) -> str:
        if method == "beta_binomial":
            padj_threshold, effect_threshold = thr
            return f"padj{padj_threshold}_delta{effect_threshold}"
        return f"z{thr}"

    gene_info = load_bed(args.bed)
    manifest  = load_manifest(args.manifest)

    manifest_valid = [(g, p) for g, p in manifest if p is not None]
    missing_from_bed = set(g for g, _ in manifest_valid) - set(gene_info)
    if missing_from_bed:
        print(f"[WARNING] {len(missing_from_bed)} gene(s) in manifest not in BED, skipping: {sorted(missing_from_bed)}")
    manifest_valid = [(g, p) for g, p in manifest_valid if g in gene_info]

    if args.test_n_genes is not None:
        manifest_valid = manifest_valid[:args.test_n_genes]
        print(f"[INFO] --test-n-genes {args.test_n_genes}: processing first {len(manifest_valid)} gene(s) only.")

    strand_map = {g: info[2] for g, info in gene_info.items()}

    # Parse GTF upfront if provided
    gtf_junctions: Optional[Dict[str, Dict]] = None
    if args.gtf:
        print(f"\nParsing GTF for annotated junctions (upfront) ...")
        t_gtf = time.time()
        gene_names = [g for g, _ in manifest_valid]
        gtf_junctions = parse_gtf_junctions(args.gtf, gene_names)
        n_matched = sum(1 for g in gene_names if g in gtf_junctions)
        print(f"  → matched {n_matched}/{len(gene_names)} genes ({time.time()-t_gtf:.2f}s)")
    else:
        print("[INFO] --gtf not provided; junction_type column will be omitted.")

    n_genes = len(manifest_valid)
    _metrics = ["PSI_approx"]
    if not approx_only:
        _metrics += ["PSI", "full_IR_ratio"]
        if not args.no_ss_IR:
            _metrics += ["5ss_IR_ratio", "3ss_IR_ratio"]
        if args.has_ipa:
            _metrics.append("IPA_ratio")
    print(f"\nWill process {n_genes} gene(s)")
    print(f"Metrics: {', '.join(_metrics)}")
    print(f"Threads per gene: {args.threads}\n")

    if n_genes == 0:
        print("[WARNING] Manifest has no genes with results -- this group's cohort_junction_analysis "
              "run either skipped entirely (see that rule's --note output) or every gene failed.")

    def _write_tsv(df, path): df.to_csv(path, sep="\t", index=False)

    def _empty_outlier_outputs(reason: str) -> None:
        """Writes an empty (header-only) outliers.tsv / outliers_filtered.tsv for
        every requested threshold, so Snakemake's declared output file exists
        even when there was nothing to analyze."""
        print(f"\n[WARNING] {reason}")
        empty_cols = _OUTPUT_COLS + ["event_type"]
        for thr in threshold_specs:
            thr_dir  = _threshold_dir(_thr_subdir_name(thr))
            empty_df = pd.DataFrame(columns=empty_cols)
            empty_df.to_csv(os.path.join(thr_dir, f"{prefix_name}_outliers.tsv"), sep="\t", index=False)
            empty_df.to_csv(os.path.join(thr_dir, f"{prefix_name}_outliers_filtered.tsv"), sep="\t", index=False)
        print("  Wrote empty outlier file(s) so downstream outputs still exist.")

    if n_genes == 0:
        _empty_outlier_outputs("No genes to process.")
        print("\nDone (nothing to do).")
        return

    def _finalize_results(
        all_results: List[pd.DataFrame],
        computed_metrics: List[str],
    ) -> Optional[pd.DataFrame]:
        """Concatenates every gene's per-metric results, assigns junction_type
        (if a GTF was provided), and -- for beta_binomial only -- applies FDR
        (Benjamini-Hochberg) correction across the whole cohort's p-values.
        modified_zscore has no p-values, so there's nothing to correct there;
        |modZ| is used directly as the outlier statistic."""
        if not all_results: return None

        total_rows = sum(len(r) for r in all_results)
        print(f"\n{'=' * 70}")
        if method == "beta_binomial":
            print("  Correcting p-values and identifying outliers...")
        else:
            print("  Assembling cohort-wide results...")
        print(f"{'=' * 70}")
        t_fdr = time.time()

        final_df = pd.concat(all_results, ignore_index=True)
        final_df = final_df.sort_values(
            ["gene", "junction", "sample", "phasing"], ignore_index=True
        )

        # Assign junction_type if GTF was provided
        if gtf_junctions is not None:
            final_df = assign_junction_types(final_df, gtf_junctions)

        fmt_df = final_df.copy()

        if method == "beta_binomial":
            print(f"\n  Applying FDR-BH ({total_rows:,} rows across {len(computed_metrics)} metric(s)) ...")
            for mc in computed_metrics:
                p_col    = f"p_value_{mc}"
                padj_col = f"padj_{mc}"
                if p_col not in fmt_df.columns: continue
                p_vals = pd.to_numeric(fmt_df[p_col], errors="coerce")
                is_ss  = mc in ("5ss_IR_ratio", "3ss_IR_ratio")
                if is_ss:
                    ss_pos_col = "5ss" if "5ss" in mc else "3ss"
                    dedup_key = list(zip(fmt_df["sample"], fmt_df["gene"],
                                         fmt_df["phasing"], fmt_df[ss_pos_col]))
                    seen: Dict[tuple, int] = {}
                    dedup_idx = []
                    for i, k in enumerate(dedup_key):
                        if k not in seen:
                            seen[k] = i; dedup_idx.append(i)
                    dedup_p = p_vals.iloc[dedup_idx]
                    valid_mask = dedup_p.notna()
                    padj_dedup = np.full(len(dedup_p), np.nan)
                    if valid_mask.sum() > 0:
                        _, pv, _, _ = multipletests(dedup_p[valid_mask].to_numpy(), method="fdr_bh")
                        padj_dedup[valid_mask.to_numpy()] = pv
                    key_to_padj = {k: padj_dedup[j] for j, k in enumerate(
                        [dedup_key[i] for i in dedup_idx])}
                    p_str = fmt_df[p_col]
                    fmt_df[padj_col] = [
                        key_to_padj.get(k, p_str.iloc[i])
                        for i, k in enumerate(dedup_key)
                    ]
                else:
                    p_str  = fmt_df[p_col]
                    valid  = p_vals.notna()
                    padj_vals = p_str.copy().astype(object)
                    padj_arr  = np.full(valid.sum(), np.nan)
                    if valid.sum() > 0:
                        _, pv, _, _ = multipletests(p_vals[valid].to_numpy(), method="fdr_bh")
                        padj_arr = pv
                    padj_vals[valid] = padj_arr
                    # Rows with sentinel strings (low_n, low_coverage, etc.) should mirror p_value
                    is_sentinel = ~valid & p_str.notna()
                    padj_vals[is_sentinel] = p_str[is_sentinel]
                    fmt_df[padj_col] = padj_vals

        fmt_df = select_output_columns(fmt_df)
        print(f"       Done ({time.time()-t_fdr:.2f}s)")
        return fmt_df

    results_subdir = _results_dir()
    all_results: List[pd.DataFrame] = []

    for gene, path in manifest_valid:
        t_gene = time.time()
        try:
            combined = load_gene_raw_metrics(path)
            res = run_gene_stats_pipeline(
                gene, combined, approx_only,
                args.coverage_threshold, args.PSI_rescale_factor, args.n_threshold,
                args.threads, args.has_ipa, args.no_ss_IR, method,
            )
        except Exception as e:
            print(f"[ERROR] Gene {gene}: {e}"); traceback.print_exc()
            res = None
        if res is not None and len(res):
            all_results.append(res)
        print(f"  {gene} complete ({time.time() - t_gene:.0f}s)")

    if not all_results:
        _empty_outlier_outputs("No gene produced usable results (every gene errored, or "
                                "produced empty output, during statistical testing).")
        print("\nDone (nothing to do).")
        return

    sample_df   = all_results[0]
    test_prefix = "p_value_" if method == "beta_binomial" else "modz_"
    computed_metrics = [c.replace(test_prefix, "") for c in sample_df.columns
                        if c.startswith(test_prefix)]

    # FDR correction (beta_binomial only; done once, shared across all thresholds)
    final_df = _finalize_results(all_results, computed_metrics)
    if final_df is None:
        _empty_outlier_outputs("Result assembly produced no rows.")
        print("\nDone (nothing to do).")
        return

    # ---- Write per-gene TSVs once (threshold-independent) ----
    fmt_base = select_output_columns(final_df)
    n_genes_write = fmt_base["gene"].nunique()
    print(f"  Writing {n_genes_write} per-gene results files ...")
    t_w = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = [
            ex.submit(_write_tsv, gdf,
                      os.path.join(results_subdir, f"{g}.tsv"))
            for g, gdf in fmt_base.groupby("gene")
        ]
        for fut in concurrent.futures.as_completed(futs):
            fut.result()
    print(f"       File writes done ({time.time()-t_w:.2f}s)")

    _IR_IPA_METRICS = frozenset((
        "5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio", "junction_IPA_ratio"
    ))

    # ---- Per-threshold outlier identification, classification, and output ----
    for thr in threshold_specs:
        t_thr = time.time()

        if method == "beta_binomial":
            padj_threshold, effect_threshold = thr
            stat_label = "delta"
            thr_desc   = f"padj <= {padj_threshold}, |delta| >= {effect_threshold}"

            def _effect_col(mc): return f"delta_{mc}"

            def _metric_cols_ok(df, mc):
                return f"padj_{mc}" in df.columns and _effect_col(mc) in df.columns

            def _outlier_mask(df, mc):
                if not _metric_cols_ok(df, mc):
                    return pd.Series(False, index=df.index)
                padj_v  = pd.to_numeric(df[f"padj_{mc}"], errors="coerce")
                delta_v = pd.to_numeric(df[_effect_col(mc)], errors="coerce")
                return padj_v.le(padj_threshold) & delta_v.abs().ge(effect_threshold)
        else:
            z_threshold      = thr
            effect_threshold = z_threshold
            stat_label = "modZ"
            thr_desc   = f"|modZ| >= {z_threshold}"

            def _effect_col(mc): return f"modz_{mc}"

            def _metric_cols_ok(df, mc):
                return _effect_col(mc) in df.columns

            def _outlier_mask(df, mc):
                if not _metric_cols_ok(df, mc):
                    return pd.Series(False, index=df.index)
                zv = pd.to_numeric(df[_effect_col(mc)], errors="coerce")
                return zv.abs().ge(z_threshold)

        def _pos_mask(df, mc):
            ec = _effect_col(mc)
            if ec not in df.columns:
                return pd.Series(False, index=df.index)
            ev = pd.to_numeric(df[ec], errors="coerce")
            return ev.ge(effect_threshold)

        thr_dir = _threshold_dir(_thr_subdir_name(thr))
        print(f"\n{'='*70}")
        print(f"  Threshold: {thr_desc}")
        print(f"  Output: {thr_dir}")
        print(f"{'='*70}")

        # ---- Outlier mask (any metric) ----
        outlier_mask = pd.Series(False, index=final_df.index)
        for mc in computed_metrics:
            outlier_mask |= _outlier_mask(final_df, mc)

        sig_df = final_df[outlier_mask].copy()

        print(f"  Identifying outliers with {thr_desc} ...")
        for mc in computed_metrics:
            if not _metric_cols_ok(final_df, mc): continue
            mask_mc = _outlier_mask(final_df, mc)
            n_rows  = int(mask_mc.sum())
            n_jxns  = int(final_df[mask_mc]["junction"].nunique())
            if mc == "5ss_IR_ratio" and "5ss" in final_df.columns:
                n_ss = int(final_df[mask_mc]["5ss"].nunique())
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions, {n_ss} unique 5ss)")
            elif mc == "3ss_IR_ratio" and "3ss" in final_df.columns:
                n_ss = int(final_df[mask_mc]["3ss"].nunique())
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions, {n_ss} unique 3ss)")
            else:
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions)")

        # ---- Unreliable haplotype outlier filter per metric ----
        # A haplotype (sample, junction) pair is unreliable when EITHER:
        #   1. bulk is NOT sandwiched: NOT (hap1 <= bulk <= hap2 OR hap1 >= bulk >= hap2)
        #   2. one haplotype dominates asymmetrically:
        #      max(|hap1-bulk|, |hap2-bulk|) / min(|hap1-bulk|, |hap2-bulk|) > 10
        # If either condition holds, the pair is unreliable and outlier_{mc} is set to False.
        unreliable_hap_outliers: Dict[str, set] = {}

        if not sig_df.empty:
            key_col = list(zip(sig_df["sample"], sig_df["junction"]))
            sig_df["_key"] = key_col

            for mc in computed_metrics:
                rescaled_c = f"rescaled_{mc}"
                if rescaled_c not in final_df.columns:
                    unreliable_hap_outliers[mc] = set()
                    continue
                phasing_df = final_df[final_df["phasing"].isin(["bulk", "hap1", "hap2"])][
                    ["sample", "junction", "phasing", rescaled_c]
                ].copy()
                phasing_df[rescaled_c] = pd.to_numeric(phasing_df[rescaled_c], errors="coerce")
                pivot = phasing_df.pivot_table(
                    index=["sample", "junction"], columns="phasing",
                    values=rescaled_c, aggfunc="first"
                )
                missing = [c for c in ("bulk", "hap1", "hap2") if c not in pivot.columns]
                if missing:
                    unreliable_hap_outliers[mc] = set()
                    continue
                pivot = pivot.dropna(subset=["bulk", "hap1", "hap2"])
                if pivot.empty:
                    unreliable_hap_outliers[mc] = set()
                    continue

                b  = pivot["bulk"]
                h1 = pivot["hap1"]
                h2 = pivot["hap2"]

                # Condition 1: bulk is sandwiched between hap1 and hap2
                sandwiched = ((h1 <= b) & (b <= h2)) | ((h2 <= b) & (b <= h1))

                # Condition 2: ratio of larger to smaller deviation <= 10
                d1   = (b - h1).abs()
                d2   = (b - h2).abs()
                dmax = np.maximum(d1, d2)
                dmin = np.minimum(d1, d2)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(dmin > 0, dmax / dmin, np.inf)
                symmetric = ratio <= 10

                unreliable_hap_outliers[mc] = set(pivot.index[~sandwiched | ~symmetric])

            sig_df = sig_df.drop(columns=["_key"])
            print(f"  Removing unreliable haplotype-only outliers ...")
            for mc in computed_metrics:
                unreliable = unreliable_hap_outliers.get(mc, set())
                if not unreliable or not _metric_cols_ok(sig_df, mc): continue
                hap_rows_mc = sig_df[sig_df["phasing"].isin(["hap1", "hap2"])]
                passes_mc = _outlier_mask(hap_rows_mc, mc)
                keys = list(zip(
                    hap_rows_mc.loc[passes_mc, "sample"],
                    hap_rows_mc.loc[passes_mc, "junction"],
                ))
                n_tot = sum(1 for k in keys if k in unreliable)
                print(f"       {mc}: removed {n_tot:,} unreliable haplotype outlier rows")

            n_total_jxns = sig_df["junction"].nunique() if not sig_df.empty else 0
            print(f"  {len(sig_df):,} total outlier rows ({n_total_jxns} unique junctions)")
        else:
            unreliable_hap_outliers = {mc: set() for mc in computed_metrics}

        # ---- outlier_{metric} boolean columns ----
        sig_df = sig_df.copy()
        for mc in computed_metrics:
            ocol = f"outlier_{mc}"
            if not _metric_cols_ok(sig_df, mc):
                sig_df[ocol] = False; continue
            passes     = _outlier_mask(sig_df, mc)
            unreliable = unreliable_hap_outliers.get(mc, set())
            if unreliable:
                is_hap        = sig_df["phasing"].isin(["hap1", "hap2"])
                keys          = list(zip(sig_df["sample"], sig_df["junction"]))
                is_unreliable = pd.Series([k in unreliable for k in keys], index=sig_df.index)
                sig_df[ocol]  = passes & ~(is_hap & is_unreliable)
            else:
                sig_df[ocol] = passes

        # Per-metric breakdown of final outlier counts
        _has_jxn_type = "junction_type" in sig_df.columns
        for mc in computed_metrics:
            ocol = f"outlier_{mc}"
            if ocol not in sig_df.columns: continue
            mask_mc = sig_df[ocol].astype(bool)
            n_rows  = int(mask_mc.sum())
            n_jxns  = int(sig_df[mask_mc]["junction"].nunique())
            is_5ss  = mc == "5ss_IR_ratio" and "5ss" in sig_df.columns
            is_3ss  = mc == "3ss_IR_ratio" and "3ss" in sig_df.columns
            ss_col  = "5ss" if is_5ss else ("3ss" if is_3ss else None)
            ss_lbl  = "5ss" if is_5ss else ("3ss" if is_3ss else None)
            if mc in _IR_IPA_METRICS and _effect_col(mc) in sig_df.columns:
                pos_mask = mask_mc & _pos_mask(sig_df, mc)
                n_pos    = int(sig_df[pos_mask]["junction"].nunique())
                ss_total = f", {int(sig_df[mask_mc][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                ss_pos   = f", {int(sig_df[pos_mask][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                if _has_jxn_type:
                    can_ann = sig_df["junction_type"].isin(["canonical", "annotated"])
                    n_can   = int(sig_df[pos_mask & can_ann]["junction"].nunique())
                    ss_can  = f", {int(sig_df[pos_mask & can_ann][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                    print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions{ss_total} → {n_pos} unique junctions{ss_pos} with {stat_label} > 0 → {n_can} unique junctions{ss_can} canonical or annotated)")
                else:
                    print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions{ss_total} → {n_pos} unique junctions{ss_pos} with {stat_label} > 0)")
            elif ss_col:
                n_ss = int(sig_df[mask_mc][ss_col].nunique())
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions, {n_ss} unique {ss_lbl})")
            else:
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions)")

        # ---- Compute n_sample_outlier for outlier TSV (threshold-specific) ----
        fmt_updated = final_df.copy()
        for mc in computed_metrics:
            n_col = f"n_sample_outlier_{mc}"
            ocol  = f"outlier_{mc}"
            if sig_df.empty or ocol not in sig_df.columns:
                fmt_updated[n_col] = 0
                continue
            counts = (
                sig_df[sig_df[ocol].astype(bool)]
                .groupby(["gene", "junction"])["sample"]
                .nunique().rename(n_col).reset_index()
            )
            fmt_updated = fmt_updated.merge(counts, on=["gene", "junction"], how="left")
            fmt_updated[n_col] = fmt_updated[n_col].fillna(0).astype(int)

        # ---- Event classification ----
        n_sig = len(sig_df)
        print(f"  Outliers identified ({time.time() - t_thr:.2f}s)")
        print(f"  Classifying events ...")
        t_cls = time.time()
        if n_sig > 0:
            sig_df = classify_all_events(
                sig_df, final_df, strand_map,
                effect_threshold, args.threads, args.has_ipa,
                computed_metrics, _effect_col, unreliable_hap_outliers,
            )
        else:
            sig_df["event_type"] = "none"
        print(f"       done ({time.time()-t_cls:.2f}s)  "
              f"{sig_df['event_type'].value_counts().to_dict() if n_sig else {}}")

        # ---- Sort outliers ----
        if n_sig > 0:
            effect_cols = [c for c in sig_df.columns
                          if c.startswith("delta_") or c.startswith("modz_")]
            if effect_cols:
                effect_num = sig_df[effect_cols].apply(lambda col: pd.to_numeric(col, errors="coerce"))
                sig_df["_max_abs_delta"] = effect_num.abs().max(axis=1)

                # Per-row: does it have a named (non-'other') event, and what's the
                # max abs effect value of the metric(s) that produced it?
                all_metric_event_sets = {
                    mc: frozenset(evs) for mc, evs in _METRIC_EVENTS.items()
                }
                def _named_event_delta(row):
                    et  = row.get("event_type", "other") or "other"
                    evs = set(e.strip() for e in et.split(",")) - {"other", ""}
                    if not evs:
                        return 0, 0.0
                    best = 0.0
                    for mc, mc_evs in all_metric_event_sets.items():
                        if evs & mc_evs:
                            dc = _effect_col(mc)
                            if dc in row:
                                v = pd.to_numeric(row[dc], errors="coerce")
                                if pd.notna(v):
                                    best = max(best, abs(v))
                    return 1, best

                if "event_type" in sig_df.columns:
                    named_info = sig_df.apply(_named_event_delta, axis=1, result_type="expand")
                    sig_df["_has_named_event"]   = named_info[0]
                    sig_df["_named_event_delta"]  = named_info[1]
                else:
                    sig_df["_has_named_event"]   = 0
                    sig_df["_named_event_delta"]  = 0.0

                gene_agg = (
                    sig_df.groupby(["sample", "gene"])
                    .agg(
                        _gene_has_named  =("_has_named_event",   "max"),
                        _gene_named_delta=("_named_event_delta", "max"),
                        _gene_max_delta  =("_max_abs_delta",     "max"),
                    )
                    .reset_index()
                )
                gene_agg["_sort_key"] = list(zip(
                    -gene_agg["_gene_has_named"],
                    -gene_agg["_gene_named_delta"],
                    -gene_agg["_gene_max_delta"],
                ))
                gene_agg["_gene_rank"] = (
                    gene_agg.groupby("sample")["_sort_key"]
                    .rank(ascending=True, method="min")
                )
                sig_df = sig_df.merge(
                    gene_agg[["sample", "gene", "_gene_rank"]],
                    on=["sample", "gene"], how="left"
                )
                sig_df = sig_df.sort_values(
                    ["sample", "_gene_rank", "junction"],
                    ascending=[True, True, True],
                ).drop(columns=["_max_abs_delta", "_has_named_event", "_named_event_delta"])
                sig_df = sig_df.rename(columns={"_gene_rank": "gene_rank"})

        # ---- Build outlier TSV column order ----
        # _OUTPUT_COLS + n_sample_outlier_{metric} + outlier_{metric} after each padj_{metric}
        # (beta_binomial) or modz_{metric} (modified_zscore) + event_type at end.
        # Merge n_sample_outlier columns from fmt_updated into sig_df
        n_sample_cols = [f"n_sample_outlier_{mc}" for mc in computed_metrics]
        n_sample_cols_present = [c for c in n_sample_cols if c in fmt_updated.columns]
        if n_sample_cols_present:
            sig_df = sig_df.merge(
                fmt_updated[["sample", "gene", "junction", "phasing"] + n_sample_cols_present]
                .drop_duplicates(subset=["sample", "gene", "junction", "phasing"]),
                on=["sample", "gene", "junction", "phasing"], how="left"
            )

        outlier_anchor_col = "padj_" if method == "beta_binomial" else "modz_"
        outlier_tsv_cols = []
        for col in _OUTPUT_COLS:
            outlier_tsv_cols.append(col)
            for mc in computed_metrics:
                if col == f"{outlier_anchor_col}{mc}":
                    n_sc = f"n_sample_outlier_{mc}"
                    if n_sc in sig_df.columns:
                        outlier_tsv_cols.append(n_sc)
                    outlier_tsv_cols.append(f"outlier_{mc}")
        outlier_tsv_cols.append("event_type")
        outlier_tsv_cols = [c for c in outlier_tsv_cols if c in sig_df.columns]

        # ---- Parallel output: TSVs + heatmaps + box plots ----
        t_out = time.time()
        out_jobs = []

        # TSV jobs — need the full prepared slice
        _outliers_data = sig_df[outlier_tsv_cols].copy()
        if "event_type" in sig_df.columns:
            filt = sig_df[sig_df["event_type"] != "none"]
            _outliers_filt_data = filt[outlier_tsv_cols].copy() if not filt.empty else pd.DataFrame(columns=outlier_tsv_cols)
        else:
            filt = sig_df.iloc[0:0]
            _outliers_filt_data = pd.DataFrame(columns=outlier_tsv_cols)

        def _write_outliers(data=_outliers_data):
            t0 = time.time()
            out = os.path.join(thr_dir, f"{prefix_name}_outliers.tsv")
            data.to_csv(out, sep="\t", index=False)
            print(f"  Outliers → {out} ({time.time()-t0:.2f}s)")

        def _write_outliers_filtered(data=_outliers_filt_data):
            t0 = time.time()
            out_filt = os.path.join(thr_dir, f"{prefix_name}_outliers_filtered.tsv")
            data.to_csv(out_filt, sep="\t", index=False)
            print(f"  Outliers (filtered) → {out_filt} ({time.time()-t0:.2f}s)")


        # ---- Build outlier_map for box plots ----
        # Built from filt (outliers_filtered rows: event_type != "none", outlier_{mc} True),
        # so boxplots exactly match outliers_filtered.tsv.
        outlier_map: Dict[str, Dict[str, Dict[str, set]]] = {}
        for mc in computed_metrics:
            ocol = f"outlier_{mc}"
            gene_jxn_map: Dict[str, Dict[str, set]] = defaultdict(lambda: defaultdict(set))
            if not filt.empty and ocol in filt.columns:
                passing = filt[filt[ocol].astype(bool)]
                for _, row in passing.iterrows():
                    gene_jxn_map[row["gene"]][row["junction"]].add(row["sample"])
            outlier_map[mc] = gene_jxn_map

        out_jobs.append(_write_outliers)
        out_jobs.append(_write_outliers_filtered)

        # Heatmap + box plot jobs per metric — pre-slice to only needed columns/rows
        if n_sig > 0:
            for mc in computed_metrics:
                dc = _effect_col(mc)
                if dc not in sig_df.columns:
                    continue

                # Heatmap: only gene, sample, effect col — rows passing threshold
                heat_mask = _outlier_mask(sig_df, mc)
                heat_df   = sig_df.loc[heat_mask, ["gene", "sample", dc]].copy()

                def _make_heatmap(df=heat_df, mc=mc, dc=dc):
                    make_outlier_heatmap(
                        df, mc, dc, stat_label, thr_desc,
                        os.path.join(thr_dir, f"{prefix_name}_outlier_heatmap_{mc}.pdf"),
                    )
                out_jobs.append(_make_heatmap)

                rc = f"rescaled_{mc}"
                if rc in fmt_updated.columns:
                    # Box plots: only gene, junction, phasing, sample, rescaled col
                    # filtered to genes/junctions in outlier_map[mc]
                    gj_map = outlier_map.get(mc, {})
                    if gj_map:
                        bp_genes = set(gj_map.keys())
                        bp_jxns  = set(j for jxns in gj_map.values() for j in jxns)
                        bp_df = fmt_updated.loc[
                            fmt_updated["gene"].isin(bp_genes) &
                            fmt_updated["junction"].isin(bp_jxns),
                            ["gene", "junction", "phasing", "sample", rc]
                        ].copy()
                        def _make_boxplot(df=bp_df, mc=mc, rc=rc):
                            make_hit_boxplots(
                                df, outlier_map, mc, rc, thr_dir, prefix_name, tmp_dir,
                                threshold_desc=thr_desc,
                                stat_label=stat_label,
                                effect_threshold=effect_threshold,
                                is_ir_ipa=(mc in _IR_IPA_METRICS),
                                has_jxn_type_filter=("junction_type" in fmt_updated.columns
                                                     and mc in _IR_IPA_METRICS),
                            )
                        out_jobs.append(_make_boxplot)

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = [ex.submit(fn) for fn in out_jobs]
            for fut in concurrent.futures.as_completed(futs):
                try:
                    fut.result()
                except Exception as e:
                    print(f"[WARNING] Output job failed: {e}"); traceback.print_exc()

        print(f"  All outputs written ({time.time()-t_out:.2f}s)")

        if os.path.exists(tmp_dir):
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)

        # ---- Outlier summary for this threshold ----
        try:
            _METRIC_DISPLAY = {
                "junction_PSI":            "PSI",
                "junction_PSI_approx":     "PSI_approx",
                "5ss_IR_ratio":            "5ss_IR",
                "3ss_IR_ratio":            "3ss_IR",
                "junction_full_IR_ratio":  "full_IR",
                "junction_IPA_ratio":      "IPA",
            }
            _METRIC_ORDER   = list(_METRIC_DISPLAY.keys())
            computed_mets   = set(computed_metrics)
            metrics_present = [m for m in _METRIC_ORDER if m in computed_mets]

            if "sample" in sig_df.columns and "gene" in sig_df.columns:
                rows_any = sig_df[["sample", "gene"]].drop_duplicates()
                any_ser  = rows_any.groupby("sample")["gene"].nunique().rename("_any")
                metric_sers = []
                for mc in metrics_present:
                    ocol = f"outlier_{mc}"
                    if ocol not in sig_df.columns: continue
                    ms = (sig_df[sig_df[ocol].astype(bool)][["sample", "gene"]]
                          .drop_duplicates()
                          .groupby("sample")["gene"].nunique().rename(mc))
                    metric_sers.append(ms)
                summary = pd.concat([any_ser] + metric_sers, axis=1)
                summary = summary.fillna(0)
                for col in summary.columns:
                    summary[col] = pd.to_numeric(summary[col], errors="coerce").fillna(0).astype(int)
                summary = summary.sort_values("_any", ascending=False).reset_index()
                header_metrics = ["Genes"] + [_METRIC_DISPLAY[m] for m in metrics_present
                                               if m in summary.columns]
                col_w  = max(10, max(len(h) for h in header_metrics) + 2)
                header = f"  {'Sample':<40}" + "".join(f"{h:>{col_w}}" for h in header_metrics)
                sep    = f"  {'-'*40}" + ("-"*col_w) * len(header_metrics)
                print("\n" + "="*len(sep.rstrip()))
                print(f"  Outlier summary — {thr_desc}")
                print(f"  (samples ranked by total genes with outlier)")
                print("="*len(sep.rstrip()))
                print(header); print(sep)
                for _, row in summary.iterrows():
                    vals = f"  {row['sample']:<40}{row['_any']:>{col_w}}"
                    for m in metrics_present:
                        if m in summary.columns:
                            vals += f"{row.get(m, 0):>{col_w}}"
                    print(vals)
                print("="*len(sep.rstrip()))
        except Exception as e:
            print(f"[WARNING] Outlier summary failed: {e}"); traceback.print_exc()

        print(f"\n  Finished threshold: {thr_desc} ({time.time()-t_thr:.2f}s)")

    # end threshold loop

    # ---- QC figures (uses already-parsed gtf_junctions) ----
    if args.gtf and gtf_junctions is not None:
        gene_names = [g for g, _ in manifest_valid]
        print(f"\n{'='*56}")
        print(f"  QC")
        print(f"{'='*56}")
        t_qc = time.time()
        print(f"\n  {'Gene':<20} {'Canonical jxns':>16} {'Annotated jxns':>16}")
        print(f"  {'-'*20} {'-'*16} {'-'*16}")
        for g in gene_names:
            if g in gtf_junctions:
                n_can = len(gtf_junctions[g]["canonical_junctions"])
                n_ann = len(gtf_junctions[g]["all_junctions"])
                print(f"  {g:<20} {n_can:>16} {n_ann:>16}")
            else:
                print(f"  {g:<20} {'NOT FOUND':>16} {'NOT FOUND':>16}")

        results_subdir = _results_dir()

        # Use final_df directly — bulk rows only, pre-sliced per figure
        os.makedirs(qc_dir, exist_ok=True)
        bulk_df = final_df[final_df["phasing"] == "bulk"]

        fit_prefix   = "alpha_" if method == "beta_binomial" else "median_"
        not_fittable = ("low_n", "error") if method == "beta_binomial" else ("low_n", "no_variance")

        qc_jobs = []
        for cov_col, file_suffix, bar_metrics in _QC_FIGURES:
            if approx_only and cov_col != "junction_coverage_approx":
                continue
            if args.no_ss_IR and cov_col in ("5ss_coverage", "3ss_coverage"):
                continue
            companions = [m for m in bar_metrics
                          if args.has_ipa or m != "junction_IPA_ratio"]
            if not companions:
                continue
            needed_cols = ["gene", "sample", "junction", cov_col]
            for cm in companions:
                ac = f"{fit_prefix}{cm}"
                if ac in bulk_df.columns:
                    needed_cols.append(ac)
            needed_cols = list(dict.fromkeys(c for c in needed_cols if c in bulk_df.columns))
            qc_slice = bulk_df[needed_cols].copy()

            for subset in ("canonical", "annotated"):
                out_pdf = os.path.join(qc_dir, f"{prefix_name}_qc_{file_suffix}_{subset}.pdf")
                def _qc_job(df=qc_slice, cc=cov_col, cm=companions, fs=file_suffix,
                            ss=subset, op=out_pdf):
                    make_qc_figure(
                        df, gtf_junctions,
                        cc, cm, fs, ss, args.coverage_threshold, op,
                        args.has_ipa, fit_prefix, not_fittable,
                    )
                qc_jobs.append(_qc_job)

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
                futs = [ex.submit(fn) for fn in qc_jobs]
                for fut in concurrent.futures.as_completed(futs):
                    try:
                        fut.result()
                    except Exception as e:
                        print(f"[WARNING] QC figure failed: {e}"); traceback.print_exc()

        print(f"\n  Finished QC ({time.time()-t_qc:.2f}s)")

    print("\nDone.")


if __name__ == "__main__":
    main()