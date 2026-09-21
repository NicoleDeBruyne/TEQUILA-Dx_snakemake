
import numpy as np
import pandas as pd


def compute_outlier_scores(cptm_df, pseudocount=1.0, shrinkage_k=10.0, min_mad=0.1,
                            n_bins=None, z_threshold=3.0):
    L = np.log2(cptm_df.astype(float) + pseudocount)
    n_genes, n_samples = L.shape

    if n_genes == 0 or n_samples < 2:
        nan_df = pd.DataFrame(np.nan, index=cptm_df.index, columns=cptm_df.columns)
        return nan_df, nan_df.notna()

    if n_bins is None:
        n_bins = max(1, min(20, n_genes // 5))

    gene_median_all = L.median(axis=1, skipna=True)
    gene_mad_all = 1.4826 * L.sub(gene_median_all, axis=0).abs().median(axis=1, skipna=True)

    rank_order = gene_median_all.rank(method="first", na_option="bottom").astype(int) - 1
    bin_of_gene = pd.Series(np.minimum(rank_order.to_numpy() * n_bins // max(n_genes, 1), n_bins - 1),
                             index=L.index)
    trend_by_bin = gene_mad_all.groupby(bin_of_gene).median()
    trend_mad = bin_of_gene.map(trend_by_bin)

    n_per_gene = L.notna().sum(axis=1)
    weight = n_per_gene / (n_per_gene + shrinkage_k)

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
    stacked = z_df.stack()
    stacked = stacked.dropna()
    stacked.index.names = ["gene", "sample"]
    mask = is_outlier_df.stack().reindex(stacked.index).fillna(False)
    out = stacked[mask].reset_index()
    out.columns = ["gene", "sample", "zscore"]
    out = out[["sample", "gene", "zscore"]].sort_values("zscore", ascending=True).reset_index(drop=True)
    return out
