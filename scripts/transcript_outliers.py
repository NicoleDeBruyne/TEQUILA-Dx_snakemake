# isoform-usage outlier calls + sample QC, shared by 9G and the report so they can't drift apart

import numpy as np
import pandas as pd

QC_PASS = "pass"
QC_FAIL = "fail"


def _num(df, col):
    return pd.to_numeric(df[col], errors="coerce")


# two flavours: isoform outliers (usage shift + padj/z) and landscape outliers (TVD)
def call_outliers(iso, land, method, padj=0.05, z=3.5, delta=0.1, tvd=0.1):
    big_shift = _num(iso, "delta_usage").abs() >= delta
    big_tvd = _num(land, "nearest_tvd") >= tvd
    if method == "beta_binomial":
        iso_out = iso[big_shift & (_num(iso, "padj") <= padj)].sort_values(["sample", "padj", "gene"])
        land_out = land[big_tvd & (_num(land, "padj") <= padj)].sort_values(["sample", "padj", "gene"])
    else:
        iso_out = iso[big_shift & (_num(iso, "modz").abs() >= z)]
        iso_out = iso_out.assign(_k=-_num(iso_out, "modz").abs()).sort_values(["sample", "_k", "gene"]).drop(columns="_k")
        land_out = land[big_tvd & (_num(land, "modz") >= z)]
        land_out = land_out.assign(_k=-_num(land_out, "modz")).sort_values(["sample", "_k", "gene"]).drop(columns="_k")
    return iso_out, land_out


def outlier_gene_counts(iso_out, land_out, samples):
    pairs = pd.concat([iso_out[["sample", "gene"]], land_out[["sample", "gene"]]], ignore_index=True)
    counts = pairs.drop_duplicates().groupby("sample")["gene"].nunique()
    return counts.reindex([str(s) for s in samples], fill_value=0)


# median + k*MAD on outlier counts, anything above looks like a noisy sample
def qc_threshold(counts, mad_multiplier, mad_floor=1.0):
    values = np.asarray(counts, dtype=float)
    med = float(np.median(values)) if len(values) else 0.0
    mad = float(np.median(np.abs(values - med))) if len(values) else 0.0
    return med + mad_multiplier * max(mad, mad_floor), med, mad


def load_sample_qc(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    return {str(s): str(q) for s, q in zip(df["sample"], df["qc"])}
