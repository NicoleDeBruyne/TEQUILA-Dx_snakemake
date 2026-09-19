#!/usr/bin/env python3
"""
scripts/motr.py
Shared MOTR ("median of target ratios") size-factor computation, used by
all quantification methods that produce a targeted-panel MOTR matrix
(scripts/quantify_gene_by_assignment.py, scripts/quantify_gene_expression.py,
scripts/normalize_amalgam_matrix.py) so the algorithm and CLI args stay in
one place instead of being copy-pasted per method.

MOTR is DESeq2's own "poscounts" median-of-ratios size-factor normalization
(their fix for the classic method's all-or-nothing zero-count handling),
restricted to genes on a run's BED panel rather than the whole
transcriptome, since that's the only gene set this pipeline ever reports
on. See compute_size_factors()'s docstring for the algorithm and its
per-sample exclusion rule.
"""

import numpy as np


def add_motr_args(parser):
    """Shared CLI option for MOTR's per-sample exclusion threshold -- same
    default/semantics across every quantification method. See
    compute_size_factors()'s docstring."""
    parser.add_argument('--motr-max-zero-fraction', type=float, default=0.25,
        help="If a given SAMPLE has a zero count for at least this fraction of the targeted-"
             "panel genes in the matrix, that sample is excluded from MOTR entirely (its size "
             "factor is set to NaN, so every gene's MOTR value for that sample comes out NaN) "
             "rather than risk an unstable estimate from too little evidence. "
             "Default: 0.25 (25%% zero)")


def compute_size_factors(raw_df, max_zero_fraction=0.25):
    """DESeq2's "poscounts" median-of-ratios size factor (their own fix for
    the classic method's all-or-nothing zero-count handling -- see
    estimateSizeFactorsForMatrix(type="poscounts") in the DESeq2 R source),
    restricted to whatever gene set raw_df already contains (the caller
    decides that scope -- here, the targeted BED panel only, hence "MOTR":
    median of TARGET ratios).

    Why not the classic method (drop any gene with a zero ANYWHERE) or a
    flat +1 pseudocount: the classic rule can drop most/all genes on a
    smaller panel with many samples (one sample's true zero for a gene
    kills that gene as a reference for EVERY sample). A flat pseudocount
    avoids that, but a low-count gene's ratio gets disproportionately
    shrunk toward 1 by adding the same constant a high-count gene barely
    notices -- and since every gene's ratio feeds into one shared
    per-sample size factor (via the median), that shrinkage biases the
    size factor for genes that had nothing to do with it. Poscounts avoids
    both: nothing is ever added to a real count, and only genes that are
    zero in EVERY sample (no evidence at all) are excluded.

    Algorithm:
      1. Per gene, geometric mean using only that gene's NONZERO samples,
         in log-space: sum(log(nonzero counts)) / (total number of
         samples) -- note the denominator is every sample, not just the
         nonzero ones, matching DESeq2's own geoMeanNZ(). A gene is
         excluded only if it's zero in every sample (nothing to reference).
      2. Per sample, per gene: log(this sample's raw count) - log(that
         gene's geometric mean) -- using the REAL count, not a pseudocount.
         A gene that's zero in this particular sample simply contributes no
         ratio for this sample (skipped, not treated as a fabricated
         value) -- it may still be usable for other samples where it's
         nonzero.
      3. Each sample's size factor = the median of its own contributed
         ratios (from step 2). Taking the median (not the mean) is what
         makes this robust to a handful of genes that genuinely differ in
         expression between samples, unlike a simple total-count
         normalization.
      4. Normalized value = raw count / size factor.

    Per-sample exclusion: if a sample has a zero count for at least
    max_zero_fraction of the (targeted-panel) genes in raw_df, that ONE
    sample is excluded from MOTR entirely -- its size factor is set to
    NaN, so raw_df.div(size_factors, axis=1) naturally NaNs out every
    gene's MOTR value for that sample -- with a warning, rather than
    either fabricating a number from too little evidence or (as the
    classic method would) letting one sparse sample sink every gene for
    the whole cohort.

    Returns a pandas Series of size factors, indexed by raw_df's columns
    (samples). A sample excluded per the above has NaN in this Series.
    """
    nonzero = raw_df.gt(0)
    n_samples = raw_df.shape[1]
    n_genes = raw_df.shape[0]

    # Step 1: per-gene geometric mean over nonzero samples only, denominator
    # is every sample (DESeq2's geoMeanNZ convention) -- excludes only
    # all-zero genes (sum of an empty log-set is 0, so an all-zero gene's
    # "geometric mean" comes out as log(1)=0 in this formula; explicitly
    # mark those as unusable rather than silently including them as if
    # they were legitimately flat).
    log_counts_where_nonzero = np.log(raw_df.where(nonzero))
    gene_has_any_nonzero = nonzero.any(axis=1)
    log_geomeans = log_counts_where_nonzero.sum(axis=1) / n_samples
    log_geomeans = log_geomeans.where(gene_has_any_nonzero)  # NaN for all-zero genes -> excluded below

    # Step 2/3: per sample, ratio only over genes usable for THIS sample
    # (nonzero here AND has a valid geomean), then median.
    log_ratios = log_counts_where_nonzero.sub(log_geomeans, axis=0)  # NaN wherever this sample is 0, or gene is all-zero

    size_factors = np.exp(log_ratios.median(axis=0, skipna=True))

    # Per-sample exclusion: fraction of zero counts is measured against
    # ALL genes in raw_df (the full targeted panel), not just the ones
    # that ended up contributing a ratio -- this is a statement about how
    # sparse the sample itself is, independent of which other genes
    # happened to be all-zero cohort-wide.
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
