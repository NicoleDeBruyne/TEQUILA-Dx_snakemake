# small bits that used to be duplicated across scripts

import io
from typing import List

import numpy as np
import pandas as pd

HIT_COL = "hit"


def is_hit(values):
    return values.astype(str).str.strip().str.upper() == "TRUE"


# only the rows flagged as hits
def read_hit_rows(path, **read_csv_kwargs):
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if HIT_COL in raw.columns:
        raw = raw[is_hit(raw[HIT_COL])]
    buf = io.StringIO(raw.to_csv(sep="\t", index=False))
    read_csv_kwargs.setdefault("sep", "\t")
    return pd.read_csv(buf, **read_csv_kwargs)


# PSI, but NaN where the junction or its splice sites don't have enough coverage
def compute_rescaled_psi_column(jxn_counts, ss1, ss2, PSI_rescale_factor, gtex_coverage_threshold):
    jxn_counts = pd.to_numeric(pd.Series(jxn_counts), errors='coerce').fillna(0).astype(int)

    ss1_sums = jxn_counts.groupby(ss1).sum()
    ss2_sums = jxn_counts.groupby(ss2).sum()

    ss1_coverage = pd.Series(ss1).map(ss1_sums).fillna(0).to_numpy()
    ss2_coverage = pd.Series(ss2).map(ss2_sums).fillna(0).to_numpy()

    jxn_counts = jxn_counts.to_numpy(dtype=np.float64)
    jxn_coverage = ss1_coverage + ss2_coverage - jxn_counts

    with np.errstate(divide='ignore', invalid='ignore'):
        psi = np.where(jxn_coverage == 0, np.nan, jxn_counts / jxn_coverage)

    rescaled_psi = psi * (1 - 2 * PSI_rescale_factor) + PSI_rescale_factor
    rescaled_psi[jxn_coverage < gtex_coverage_threshold] = np.nan

    return rescaled_psi


def _per_metric_output_cols(metric_col: str) -> List[str]:
    return [
        f"n_{metric_col}",
        f"alpha_{metric_col}", f"beta_{metric_col}", f"expected_{metric_col}",
        f"p1_{metric_col}", f"p99_{metric_col}",
        f"delta_{metric_col}", f"p_value_{metric_col}", f"padj_{metric_col}",
        f"median_{metric_col}", f"mad_{metric_col}", f"modz_{metric_col}",
    ]


# column order for the cohort junction outputs
OUTPUT_COLS = [
    "sample", "gene", "gene_rank", "region", "phasing", "junction", "5ss", "3ss",
    "junction_type",
    "junction_usage",
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
    present = [c for c in OUTPUT_COLS if c in df.columns]
    return df[present]
