
from __future__ import annotations

import os
import sys
import argparse
import warnings
import traceback
import time
import math
import glob as _glob
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Fits distributions, runs statistical tests, applies FDR correction, "
                     "and writes a combined scored TSV for cohort junction analysis (step 8B)."
    )
    p.add_argument("--manifest",                   required=True,
                   help="gene -> result-file manifest TSV from cohort_junction_analysis.py.")
    p.add_argument("--bed",                        required=True)
    p.add_argument("--outprefix",                  required=True)
    p.add_argument("--approx",                     action="store_true")
    p.add_argument("--has-ipa",                    action="store_true",
                   help="Set if --genome was provided to cohort_junction_analysis.py "
                        "(enables testing junction_IPA_ratio as a metric).")
    p.add_argument("--coverage-threshold",         type=int,   default=20)
    p.add_argument("--phasing-threshold",          type=float, default=0.5,
                   help="Minimum fraction of bulk coverage that must be phased "
                        "(hap1+hap2 denom / bulk denom) for haplotype metrics to be scored.")
    p.add_argument("--PSI-rescale-factor",         type=float, default=1e-3)
    p.add_argument("--n-threshold",                type=int,   default=30)
    p.add_argument("--method",                     required=True,
                   choices=["beta_binomial", "modified_zscore"],
                   help="Statistical method to use for scoring.")
    p.add_argument("--no-ss-IR",                   action="store_true")
    p.add_argument("--gtf",                        default=None,
                   help="GTF/GTF.gz. If provided, adds junction_type column.")
    p.add_argument("--threads",                    type=int,   default=1)
    p.add_argument("--test-n-genes",               type=int,   default=None)
    return p.parse_args()


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
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (ValueError, TypeError):
            pass
    return df


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
    mat_block, tol, n_threshold = args
    return [_fit_one_beta(mat_block[i], tol, n_threshold)
            for i in range(len(mat_block))]


def fit_beta_dist_chunk(mat, feat_names, tol, n_threshold, threads: int = 1):
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
            results = []
            for block in ex.map(_fit_beta_rows, chunks):
                results.extend(block)

    return pd.DataFrame(results, index=feat_names,
                        columns=["n", "alpha", "beta_param", "expected"])


def _betabinom_test_rows(args):
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
            p = np.concatenate(list(ex.map(_betabinom_test_rows, chunks)))

    df[p_col] = p
    return df


_MODZ_CONST = 0.6745
_MEANAD_CONST = 0.7979


def _fit_one_modz(x: np.ndarray, tol: float, n_threshold: int):
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
    mat_block, tol, n_threshold = args
    return [_fit_one_modz(mat_block[i], tol, n_threshold)
            for i in range(len(mat_block))]


def fit_modz_dist_chunk(mat, feat_names, tol, n_threshold, threads: int = 1):
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


def _run_one_metric(
    combined_df:        pd.DataFrame,
    metric_col:         str,
    rescaled_col:       str,
    usage_col:          str,
    coverage_col:       str,
    id_col:             str,
    coverage_threshold: int,
    phasing_threshold:  float,
    PSI_rescale_factor: float,
    n_threshold:        int,
    threads:            int,
    method:             str,
) -> pd.DataFrame:
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

    n_col   = f"n_{metric_col}"
    p1_col  = f"p1_{metric_col}"
    p99_col = f"p99_{metric_col}"
    delta_col = f"delta_{metric_col}"
    if method == "beta_binomial":
        alpha_col    = f"alpha_{metric_col}"
        beta_col     = f"beta_{metric_col}"
        expected_col = f"expected_{metric_col}"
        p_col        = f"p_value_{metric_col}"
        fit_indicator_col = alpha_col
        empty_cols = (alpha_col, beta_col, expected_col, p1_col, p99_col, delta_col, p_col)
    else:
        median_col = f"median_{metric_col}"
        mad_col    = f"mad_{metric_col}"
        modz_col   = f"modz_{metric_col}"
        fit_indicator_col = median_col
        empty_cols = (median_col, mad_col, modz_col, p1_col, p99_col, delta_col)

    if bulk_wide.empty:
        combined_df[n_col] = 0
        for col in empty_cols:
            combined_df[col] = "low_n"
        return combined_df

    mat        = bulk_wide.to_numpy(dtype=np.float32)
    feat_names = bulk_wide.index.tolist()

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
    is_error     = fit_vals == "error"
    is_no_var    = fit_vals == "no_variance"
    has_fit      = ~is_low_n & ~is_error & ~is_no_var
    coverage_arr = combined_df[coverage_col].to_numpy(dtype=float)
    has_cov      = coverage_arr >= coverage_threshold

    is_hap = combined_df["phasing"].isin(["hap1", "hap2"])
    bulk_denom = (
        combined_df[combined_df["phasing"] == "bulk"]
        [["sample", fit_id_col, coverage_col]]
        .rename(columns={coverage_col: "_bulk_denom"})
    )
    hap_denom_sum = (
        combined_df[is_hap]
        .groupby(["sample", fit_id_col])[coverage_col]
        .sum()
        .reset_index()
        .rename(columns={coverage_col: "_hap_denom_sum"})
    )
    phased_check = bulk_denom.merge(hap_denom_sum, on=["sample", fit_id_col], how="left")
    phased_check["_low_phased"] = (
        phased_check["_bulk_denom"].fillna(0) * phasing_threshold
        > phased_check["_hap_denom_sum"].fillna(0)
    )
    lp_map = phased_check.set_index(["sample", fit_id_col])["_low_phased"].to_dict()
    _lp_keys = list(zip(combined_df["sample"], combined_df[fit_id_col]))
    is_low_phased = pd.Series(
        [lp_map.get(k, False) for k in _lp_keys],
        index=combined_df.index,
        dtype=bool,
    )
    is_low_phased_hap = is_hap & is_low_phased & has_cov

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
    delta_v[is_low_n]                    = "low_n"
    delta_v[is_error]                    = "error"
    delta_v[is_no_var]                   = "no_variance"
    delta_v[has_fit & ~has_cov]          = "low_coverage"
    delta_v[has_fit & is_low_phased_hap] = "low_phased_coverage"
    combined_df[delta_col] = delta_v

    if method == "beta_binomial":
        combined_df[p_col] = pd.Series(np.nan, index=combined_df.index, dtype=object)
        combined_df.loc[is_low_n,                    p_col] = "low_n"
        combined_df.loc[is_error,                    p_col] = "error"
        combined_df.loc[has_fit & ~has_cov,          p_col] = "low_coverage"
        combined_df.loc[has_fit & is_low_phased_hap, p_col] = "low_phased_coverage"

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
        combined_df[modz_col] = pd.Series(np.nan, index=combined_df.index, dtype=object)
        combined_df.loc[is_low_n,                    modz_col] = "low_n"
        combined_df.loc[is_no_var,                    modz_col] = "no_variance"
        combined_df.loc[has_fit & ~has_cov,          modz_col] = "low_coverage"
        combined_df.loc[has_fit & is_low_phased_hap, modz_col] = "low_phased_coverage"

        testable = has_fit & has_cov & ~is_low_phased_hap
        if testable.any():
            rescaled_v2 = pd.to_numeric(combined_df.loc[testable, rescaled_col], errors="coerce")
            median_v   = pd.to_numeric(combined_df.loc[testable, median_col],   errors="coerce")
            mad_v      = pd.to_numeric(combined_df.loc[testable, mad_col],      errors="coerce")
            modz_v     = _MODZ_CONST * (rescaled_v2 - median_v) / mad_v
            combined_df.loc[testable, modz_col] = modz_v.astype(object)

    if is_ss_metric:
        combined_df = combined_df.drop(columns=["_ss_id"])

    return combined_df


def run_all_metrics(
    combined_df:        pd.DataFrame,
    coverage_threshold: int,
    phasing_threshold:  float,
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
            id_col, coverage_threshold, phasing_threshold, PSI_rescale_factor,
            n_threshold, threads, method,
        )
    return combined_df


def run_gene_stats_pipeline(
    gene:                str,
    combined:            pd.DataFrame,
    approx_only:         bool,
    coverage_threshold:  int,
    phasing_threshold:   float,
    PSI_rescale_factor:  float,
    n_threshold:         int,
    threads:             int,
    has_ipa:             bool,
    no_ss_ir:            bool,
    method:              str,
) -> pd.DataFrame:
    print(f"  Gene: {gene}  Fitting + scoring ({method}) ...")
    t0 = time.time()

    if approx_only:
        combined = _run_one_metric(
            combined, "junction_PSI_approx", "rescaled_junction_PSI_approx",
            "junction_usage", "junction_coverage_approx", "junction",
            coverage_threshold, phasing_threshold, PSI_rescale_factor, n_threshold,
            threads, method,
        )
    else:
        combined = run_all_metrics(combined, coverage_threshold, phasing_threshold,
                                   PSI_rescale_factor, n_threshold, threads, has_ipa,
                                   method, no_ss_ir)

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


def _per_metric_output_cols(metric_col: str) -> List[str]:
    return [
        f"n_{metric_col}",
        f"alpha_{metric_col}", f"beta_{metric_col}", f"expected_{metric_col}",
        f"p1_{metric_col}", f"p99_{metric_col}",
        f"delta_{metric_col}", f"p_value_{metric_col}", f"padj_{metric_col}",
        f"median_{metric_col}", f"mad_{metric_col}", f"modz_{metric_col}",
    ]


_OUTPUT_COLS = [
    "sample", "gene", "gene_rank", "region", "phasing", "junction", "5ss", "3ss",
    "junction_type",
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
    "junction_full_IR_count",
    "junction_full_IR_ratio", "rescaled_junction_full_IR_ratio",
    *_per_metric_output_cols("junction_full_IR_ratio"),
    "junction_IPA_count",
    "junction_IPA_ratio", "rescaled_junction_IPA_ratio",
    *_per_metric_output_cols("junction_IPA_ratio"),
]


def select_output_columns(df: pd.DataFrame, *_) -> pd.DataFrame:
    present = [c for c in _OUTPUT_COLS if c in df.columns]
    return df[present]


def parse_gtf_junctions(
    gtf_path: str,
    gene_names: List[str],
) -> Dict[str, Dict]:
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
            start = int(parts[3]) - 1
            end   = int(parts[4])
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
    df = df.copy()
    jtype = pd.Series("novel", index=df.index, dtype=object)

    for gene, gdf in df.groupby("gene"):
        if gene not in gtf_junctions:
            continue
        canonical_set = gtf_junctions[gene]["canonical_junctions"]
        annotated_set = gtf_junctions[gene]["all_junctions"]
        jxns = df.loc[gdf.index, "junction"]
        is_annotated = jxns.isin(annotated_set) & ~jxns.isin(canonical_set)
        is_canonical = jxns.isin(canonical_set)
        jtype.loc[gdf.index[is_annotated]] = "annotated"
        jtype.loc[gdf.index[is_canonical]] = "canonical"

    df["junction_type"] = jtype
    return df


def main() -> None:
    print("\n" + "*"*80)
    print("  Cohort Junction Fitting & Scoring")
    print("*"*80 + "\n")

    args        = parse_args()
    approx_only = args.approx
    method      = args.method

    print(f"Method: {method}")

    prefix        = args.outprefix.rstrip("/")
    outdir        = os.path.dirname(os.path.abspath(prefix))
    prefix_name   = os.path.basename(prefix)
    base_results  = os.path.join(outdir, f"{prefix_name}_results_{method}")
    sentinel_path = os.path.join(outdir, f"{prefix_name}_scoring.done")

    def _results_dir():
        os.makedirs(base_results, exist_ok=True); return base_results

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
        # Write empty outputs so downstream steps still have their inputs
        results_subdir = _results_dir()
        all_scored_path = os.path.join(results_subdir, "all_scored_results.tsv")
        pd.DataFrame(columns=_OUTPUT_COLS).to_csv(all_scored_path, sep="\t", index=False)
        with open(os.path.join(results_subdir, "method.txt"), "w") as fh:
            fh.write(method + "\n")
        with open(os.path.join(results_subdir, "computed_metrics.txt"), "w") as fh:
            pass
        with open(sentinel_path, "w") as fh:
            fh.write(f"done\nmethod={method}\nall_scored={all_scored_path}\n")
        print("\nDone (nothing to do).")
        return

    def _finalize_results(
        all_results: List[pd.DataFrame],
        computed_metrics: List[str],
    ) -> Optional[pd.DataFrame]:
        if not all_results: return None

        total_rows = sum(len(r) for r in all_results)
        print(f"\n{'=' * 70}")
        if method == "beta_binomial":
            print("  Correcting p-values and assembling results...")
        else:
            print("  Assembling cohort-wide results...")
        print(f"{'=' * 70}")
        t_fdr = time.time()

        final_df = pd.concat(all_results, ignore_index=True)
        final_df = final_df.sort_values(
            ["gene", "junction", "sample", "phasing"], ignore_index=True
        )

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
                    is_sentinel = ~valid & p_str.notna()
                    padj_vals[is_sentinel] = p_str[is_sentinel]
                    fmt_df[padj_col] = padj_vals

        fmt_df = select_output_columns(fmt_df)
        print(f"       Done ({time.time()-t_fdr:.2f}s)")
        return fmt_df

    all_results: List[pd.DataFrame] = []

    for gene, path in manifest_valid:
        t_gene = time.time()
        try:
            combined = load_gene_raw_metrics(path)
            res = run_gene_stats_pipeline(
                gene, combined, approx_only,
                args.coverage_threshold, args.phasing_threshold,
                args.PSI_rescale_factor, args.n_threshold,
                args.threads, args.has_ipa, args.no_ss_IR, method,
            )
        except Exception as e:
            print(f"[ERROR] Gene {gene}: {e}"); traceback.print_exc()
            res = None
        if res is not None and len(res):
            all_results.append(res)
        print(f"  {gene} complete ({time.time() - t_gene:.0f}s)")

    if not all_results:
        print("[WARNING] No gene produced usable results (every gene errored, or "
              "produced empty output, during statistical testing).")
        results_subdir = _results_dir()
        all_scored_path = os.path.join(results_subdir, "all_scored_results.tsv")
        pd.DataFrame(columns=_OUTPUT_COLS).to_csv(all_scored_path, sep="\t", index=False)
        with open(os.path.join(results_subdir, "method.txt"), "w") as fh:
            fh.write(method + "\n")
        with open(os.path.join(results_subdir, "computed_metrics.txt"), "w") as fh:
            pass
        with open(sentinel_path, "w") as fh:
            fh.write(f"done\nmethod={method}\nall_scored={all_scored_path}\n")
        print("\nDone (nothing to do).")
        return

    sample_df    = all_results[0]
    test_prefix  = "p_value_" if method == "beta_binomial" else "modz_"
    computed_metrics = [c.replace(test_prefix, "") for c in sample_df.columns
                        if c.startswith(test_prefix)]

    final_df = _finalize_results(all_results, computed_metrics)
    if final_df is None:
        print("[WARNING] Result assembly produced no rows.")
        results_subdir = _results_dir()
        all_scored_path = os.path.join(results_subdir, "all_scored_results.tsv")
        pd.DataFrame(columns=_OUTPUT_COLS).to_csv(all_scored_path, sep="\t", index=False)
        with open(os.path.join(results_subdir, "method.txt"), "w") as fh:
            fh.write(method + "\n")
        with open(os.path.join(results_subdir, "computed_metrics.txt"), "w") as fh:
            pass
        with open(sentinel_path, "w") as fh:
            fh.write(f"done\nmethod={method}\nall_scored={all_scored_path}\n")
        print("\nDone (nothing to do).")
        return

    results_subdir = _results_dir()

    # Write per-gene TSVs
    def _write_tsv(df, path): df.to_csv(path, sep="\t", index=False)

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

    # Write combined all_scored_results.tsv
    all_scored_path = os.path.join(results_subdir, "all_scored_results.tsv")
    print(f"  Writing combined scored results → {all_scored_path} ...")
    t_all = time.time()
    final_df.to_csv(all_scored_path, sep="\t", index=False)
    print(f"       Done ({time.time()-t_all:.2f}s, {len(final_df):,} rows)")

    # Write method.txt
    with open(os.path.join(results_subdir, "method.txt"), "w") as fh:
        fh.write(method + "\n")

    # Write computed_metrics.txt
    with open(os.path.join(results_subdir, "computed_metrics.txt"), "w") as fh:
        for mc in computed_metrics:
            fh.write(mc + "\n")

    # Write sentinel
    with open(sentinel_path, "w") as fh:
        fh.write(f"done\nmethod={method}\nall_scored={all_scored_path}\n")

    print(f"\nSentinel written → {sentinel_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
