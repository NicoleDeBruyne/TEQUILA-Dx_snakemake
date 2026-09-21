
import numpy as np


def add_motr_args(parser):
    parser.add_argument('--motr-max-zero-fraction', type=float, default=0.25,
        help="If a given SAMPLE has a zero count for at least this fraction of the targeted-"
             "panel genes in the matrix, that sample is excluded from MOTR entirely (its size "
             "factor is set to NaN, so every gene's MOTR value for that sample comes out NaN) "
             "rather than risk an unstable estimate from too little evidence. "
             "Default: 0.25 (25%% zero)")


def compute_size_factors(raw_df, max_zero_fraction=0.25):
    nonzero = raw_df.gt(0)
    n_samples = raw_df.shape[1]
    n_genes = raw_df.shape[0]

    log_counts_where_nonzero = np.log(raw_df.where(nonzero))
    gene_has_any_nonzero = nonzero.any(axis=1)
    log_geomeans = log_counts_where_nonzero.sum(axis=1) / n_samples
    log_geomeans = log_geomeans.where(gene_has_any_nonzero)

    log_ratios = log_counts_where_nonzero.sub(log_geomeans, axis=0)

    size_factors = np.exp(log_ratios.median(axis=0, skipna=True))

    zero_fraction = (~nonzero).sum(axis=0) / n_genes
    excluded = zero_fraction >= max_zero_fraction
    if excluded.any():
        for sample in raw_df.columns[excluded]:
            print(f"WARNING: sample '{sample}' had a zero count for "
                  f"{zero_fraction[sample]:.1%} of the {n_genes} targeted-panel gene(s) "
                  f"(>= {max_zero_fraction:.0%} threshold) -- excluding this sample from MOTR "
                  f"entirely (its MOTR values will be blank/NaN for every gene).")
        size_factors = size_factors.where(~excluded, np.nan)

    return size_factors
