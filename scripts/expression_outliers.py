"""
scripts/expression_outliers.py
Shared low-expression outlier score, computed identically for all four gene
quantification methods (by count, by coverage, by splice-site assignment,
and AMALGAM) so they're all on the same scale and report outliers the same
way. Deliberately NOT a parametric (e.g. negative-binomial) model: with the
sample sizes this pipeline actually sees per (bed, sample_type) group --
"kind of all over the place", per design discussion, sometimes well under
what a per-gene dispersion fit needs to be trustworthy -- a full parametric
fit tends to produce miscalibrated p-values rather than genuinely more
sensitive ones. This is a robust, nonparametric alternative: a leave-one-out
modified (median/MAD) z-score in log-space, with each gene's MAD shrunk
toward a cohort-wide MAD-vs-expression trend so that genes with too few
informative samples (or a spuriously tiny raw MAD) don't produce absurd
z-scores from ordinary noise.

One-sided by design: only low expression is scored as a possible outlier
(z <= -z_threshold), per design discussion -- unusually HIGH expression is
not flagged. The signed z-score is still returned in full (not clipped),
so nothing is thrown away for a future two-sided use.

Algorithm, per (bed, sample_type) group's CPTM-like matrix (genes x samples):
  1. Log2-transform (+pseudocount) -- expression values are multiplicative/
     right-skewed; z-scores or percentiles computed on the raw linear scale
     get dominated by the high-expression tail and compress the low-
     expression direction this score is meant to detect.
  2. Cohort-wide MAD-vs-expression trend: every gene's own (median, MAD)
     across ALL samples, binned by expression level, trend value per gene =
     the median raw MAD within its own bin. A simple, dependency-free local
     smoother (no distributional assumptions, no extra library) -- similar
     in spirit to DESeq2's dispersion trend fit, without fitting a
     parametric curve. Computed once from the full cohort (stable, since it
     pools every gene), not recomputed per leave-one-out sample below.
  3. Leave-one-out center/spread per (gene, sample): when scoring sample s
     for gene g, the median and MAD are computed from every OTHER sample in
     the group, not including s -- otherwise a single extreme value can
     visibly shift its own reference distribution, especially at small n
     (self-masking).
  4. Shrinkage: each gene's leave-one-out MAD is blended with its trend
     value, weighted by how many samples back the gene's own estimate --
     more samples -> trust the gene's own MAD more; fewer -> lean on the
     trend (same idea as empirical-Bayes dispersion shrinkage in DESeq2/
     edgeR, without the GLM). A floor (--min-mad) additionally guards
     against a near-zero MAD (e.g. a gene with nearly identical values in
     every sample) producing an artificially huge z-score from a trivial
     fluctuation.
  5. z = (this sample's log2 value - leave-one-out median) / shrunk MAD.
"""

import numpy as np
import pandas as pd


def compute_outlier_scores(cptm_df, pseudocount=1.0, shrinkage_k=10.0, min_mad=0.1,
                            n_bins=None, z_threshold=3.0):
    """
    cptm_df: genes x samples DataFrame of a normalized (CPTM-like) value.
        May contain NaN (e.g. a gene missing for some samples) -- NaN
        entries are excluded from every median/MAD calculation and their
        own z-score comes out as NaN.

    Returns (z_df, is_outlier_df), both the same shape/index/columns as
    cptm_df:
      z_df: signed robust z-score in log2 space. More negative = lower
        expression relative to the rest of the group. NOT clipped to one
        side -- the sign is informative even though only the low side is
        flagged as an outlier by is_outlier_df.
      is_outlier_df: boolean, True where z <= -abs(z_threshold) (one-sided,
        low-expression only -- see module docstring). This is a ranking/
        thresholding heuristic, not a calibrated p-value -- deliberately no
        multiple-testing correction is applied on top of it (see design
        discussion: correcting a noisy small-n statistic doesn't buy
        calibration back).
    """
    L = np.log2(cptm_df.astype(float) + pseudocount)
    n_genes, n_samples = L.shape

    if n_genes == 0 or n_samples < 2:
        # Not enough samples for a leave-one-out comparison at all.
        nan_df = pd.DataFrame(np.nan, index=cptm_df.index, columns=cptm_df.columns)
        return nan_df, nan_df.notna()  # all-False (nan_df.notna() is all False here)

    if n_bins is None:
        n_bins = max(1, min(20, n_genes // 5))

    # ---- Step 1: cohort-wide MAD-vs-expression trend (all samples, no LOO) ----
    gene_median_all = L.median(axis=1, skipna=True)
    gene_mad_all = 1.4826 * L.sub(gene_median_all, axis=0).abs().median(axis=1, skipna=True)

    # Bin genes by overall expression level; each gene's trend value is the
    # median raw MAD within its own bin -- equal-sized bins by rank, not by
    # expression value, so every bin has enough genes to give a stable
    # median even when expression is unevenly distributed across the panel.
    rank_order = gene_median_all.rank(method="first", na_option="bottom").astype(int) - 1
    bin_of_gene = pd.Series(np.minimum(rank_order.to_numpy() * n_bins // max(n_genes, 1), n_bins - 1),
                             index=L.index)
    trend_by_bin = gene_mad_all.groupby(bin_of_gene).median()
    trend_mad = bin_of_gene.map(trend_by_bin)

    # Per-gene sample count backs the shrinkage weight -- how much to trust
    # this gene's own MAD vs. the cohort-wide trend.
    n_per_gene = L.notna().sum(axis=1)
    weight = n_per_gene / (n_per_gene + shrinkage_k)

    # ---- Step 2: leave-one-out center/spread per (gene, sample), shrunk MAD ----
    z_df = pd.DataFrame(index=L.index, columns=L.columns, dtype=float)
    for sample in L.columns:
        others = L.drop(columns=sample)
        med_loo = others.median(axis=1, skipna=True)
        mad_loo = 1.4826 * others.sub(med_loo, axis=0).abs().median(axis=1, skipna=True)
        shrunk_mad = (weight * mad_loo + (1 - weight) * trend_mad).clip(lower=min_mad)
        z_df[sample] = (L[sample] - med_loo) / shrunk_mad

    is_outlier_df = z_df <= -abs(z_threshold)
    return z_df, is_outlier_df


def outliers_long_format(z_df, is_outlier_df):
    """Long-format (sample, gene, zscore) rows for every TRUE cell in
    is_outlier_df, sorted most-negative-z (most outlier-like) first --
    the shape merge_hits.py-style *_outliers.tsv files use elsewhere in
    this pipeline."""
    stacked = z_df.stack()
    stacked = stacked.dropna()
    stacked.index.names = ["gene", "sample"]
    mask = is_outlier_df.stack().reindex(stacked.index).fillna(False)
    out = stacked[mask].reset_index()
    out.columns = ["gene", "sample", "zscore"]
    out = out[["sample", "gene", "zscore"]].sort_values("zscore", ascending=True).reset_index(drop=True)
    return out
