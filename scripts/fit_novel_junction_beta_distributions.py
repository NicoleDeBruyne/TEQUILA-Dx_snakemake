#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2026.09.17

"""
scripts/fit_novel_junction_beta_distributions.py

Fits a beta distribution for the subset of ONE sample's own splice junctions
that are absent from a GTEx tissue's raw junction count matrix entirely
("novel" junctions) -- i.e. exactly the junctions
fit_gtex_beta_distributions.py's genome-wide table (rule _5B1) does not
cover, because that table only spans junctions actually present in the raw
GTEx matrix.

This still has to run per (sample, tissue), because which junctions are
"novel" is only known once you have a specific sample's own splice-junction
calls -- it can't be precomputed ahead of time the way _5B1's table can.

The fit itself, though, does NOT depend on the sample's own read counts --
only on how much coverage the novel junction's two splice sites (ss1/ss2)
have elsewhere in the GTEx reference (from other, real junctions sharing
one of those sites). This mirrors exactly what
perform_splice_junction_beta_binomial_tests.py used to do inline, per
region, via its "ADD MISSING JUNCTIONS" step (padding a novel junction into
the GTEx frame with 0 of its own reads, then letting the normal
groupby-based coverage calculation pick up whatever the shared splice sites
already have) -- computed here via the same compute_rescaled_psi_column()
helper fit_gtex_beta_distributions.py uses (kept as a compact float64
array, never accumulated as columns on the padded frame -- see that
function's docstring for why), just run once per sample across all of its
own regions rather than repeated inside every region's loop.

Output has the same schema as fit_gtex_beta_distributions.py's table (see
that module's docstring), except every row here has in_gtex_matrix=False.
"""

import argparse
import os
from math import ceil

import numpy as np
import pandas as pd
from scipy.stats import beta


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit beta distributions for one sample's splice junctions that are entirely "
                    "absent from a GTEx tissue's raw junction count matrix.")
    parser.add_argument('--jxn-info-file', required=True,
        help="Path to this sample's own splice junction counts TSV (rule _5A's output).")
    parser.add_argument('--gtexfile', required=True,
        help="Path to the raw per-GTEx-sample splice junction count matrix for this tissue.")
    parser.add_argument('--gtex-beta-fits', required=True,
        help="Path to fit_gtex_beta_distributions.py's (_5B1) output for this tissue -- used only "
             "to determine which junctions are already covered there.")
    parser.add_argument('--outfile', required=True, type=str)
    parser.add_argument('--gtex-coverage-threshold', type=int, default=20)
    parser.add_argument('--PSI-rescale-factor', type=float, default=1e-3)
    parser.add_argument('--gtex-n-threshold', type=int, default=100)
    return parser.parse_args()


def compute_rescaled_psi_column(jxn_counts, ss1, ss2, PSI_rescale_factor, gtex_coverage_threshold):
    """Identical to fit_gtex_beta_distributions.py's function of the same name -- see that
    module's docstring for why this replaces calculate_coverage()+calculate_PSI() (which used
    to add 5 new columns to the shared, genome-wide-sized padded dataframe for EVERY GTEx
    sample -- the same OOM cause _5B1 had, since this script also builds a padded frame the
    size of the whole tissue's raw matrix, just to get correct splice-site coverage for a
    handful of novel junctions).

    Returns a float64 numpy array, one entry per row of the padded frame (both the real GTEx
    junctions and the novel ones), NaN where PSI is undefined or coverage is below threshold.
    """
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


def fit_beta_dist(x, tol, n_threshold):
    """Identical to fit_gtex_beta_distributions.py / perform_splice_junction_beta_binomial_tests.py's
    function of the same name."""
    x = pd.to_numeric(x, errors='coerce')
    x = x[~np.isnan(x)]
    n = len(x)

    if n == 0:
        p1_value, p99_value = "low_n", "low_n"
    else:
        sorted_x = np.sort(x)
        if n <= 10:
            p1_value, p99_value = sorted_x[0], sorted_x[-1]
        else:
            k = ceil(n * 0.01)
            p1_value, p99_value = sorted_x[k], sorted_x[n - 1 - k]

    if n < n_threshold:
        return (n, "low_n", "low_n", "low_n", p1_value, p99_value)
    if x.var() < tol:
        return (n, x.mean() / tol, (1 - x.mean()) / tol, x.mean(), p1_value, p99_value)
    try:
        alpha_value, beta_value = beta.fit(x, floc=0, fscale=1)[0:2]
        expected_PSI = alpha_value / (alpha_value + beta_value)
        return (n, alpha_value, beta_value, expected_PSI, p1_value, p99_value)
    except Exception:
        return (n, "error", "error", "error", p1_value, p99_value)


def main():
    print("\n\n\n******************************************************************************************")
    print("Fitting beta distributions for this sample's novel (not-in-GTEx-matrix) junctions...")
    print("******************************************************************************************\n")

    args = parse_args()

    for path in (args.jxn_info_file, args.gtexfile, args.gtex_beta_fits):
        if not os.path.exists(path):
            print(f"\nERROR: {path} not found.")
            return

    os.makedirs(os.path.dirname(args.outfile) or ".", exist_ok=True)

    empty_cols = ['num_gtex_samples_with_good_coverage', 'alpha', 'beta',
                  'expected_PSI', 'p1_PSI', 'p99_PSI', 'in_gtex_matrix']

    jxn_info_df = pd.read_csv(args.jxn_info_file, sep='\t', keep_default_na=False, header=0,
                               dtype={'junction': str})

    gtex_fits = pd.read_csv(args.gtex_beta_fits, sep='\t', comment='#', index_col=0)
    known_junctions = set(gtex_fits.index)

    sample_junctions = set(jxn_info_df['junction'].unique())
    novel_junctions = sorted(sample_junctions - known_junctions)

    if not novel_junctions:
        print("\nNo novel junctions (all of this sample's junctions are already in the GTEx "
              "reference matrix) -- writing an empty output.")
        pd.DataFrame(columns=empty_cols).to_csv(args.outfile, sep='\t', index_label='junction')
        return

    print(f"\nFound {len(novel_junctions)} junction(s) in this sample absent from the GTEx "
          f"reference matrix -- fitting these against GTEx splice-site coverage.")

    print(f"\nReading GTEx file {args.gtexfile}...")
    gtex_df = pd.read_csv(args.gtexfile, sep='\t', index_col=0)

    gtex_samples = list(gtex_df.columns)
    gtex_df.columns = [c + '_jxn_alignment_count' for c in gtex_samples]

    idx_chrom_coord = gtex_df.index.str.split(':')
    gtex_chrom = idx_chrom_coord.map(lambda p: p[0])
    idx_coords = idx_chrom_coord.map(lambda p: p[1]).str.split('-')
    gtex_start = idx_coords.map(lambda p: int(p[0]))
    gtex_end = idx_coords.map(lambda p: int(p[1]))

    gtex_df.index = gtex_chrom.astype(str) + '_' + gtex_start.astype(str) + '_' + gtex_end.astype(str)
    gtex_ss1 = gtex_chrom.values + '_' + (gtex_start.values).astype(str)
    gtex_ss2 = gtex_chrom.values + '_' + (gtex_end.values).astype(str)

    ##################################################################
    # Pad the reference frame with the novel junctions, each contributing
    # 0 of its own reads to every GTEx sample -- same padding convention
    # perform_splice_junction_beta_binomial_tests.py used to apply inline,
    # per region, in its "ADD MISSING JUNCTIONS" step. Splice sites shared
    # with real junctions elsewhere in the reference still contribute their
    # own coverage via the groupby inside compute_rescaled_psi_column().
    #
    # ss1/ss2 and the padded read counts are kept as plain numpy arrays --
    # never written back as columns onto a genome-wide-sized dataframe --
    # for the same reason fit_gtex_beta_distributions.py avoids that (see
    # compute_rescaled_psi_column()'s docstring).
    ##################################################################

    novel_ss1 = np.array([j.rsplit('_', 1)[0] for j in novel_junctions])
    novel_ss2 = np.array([j.rsplit('_', 2)[0] + '_' + j.rsplit('_', 1)[1] for j in novel_junctions])

    padded_ss1 = np.concatenate([gtex_ss1, novel_ss1])
    padded_ss2 = np.concatenate([gtex_ss2, novel_ss2])
    n_gtex_rows = len(gtex_df)
    n_novel_rows = len(novel_junctions)

    print(f"Computing coverage and PSI for {n_gtex_rows} GTEx junctions + {n_novel_rows} novel "
          f"junction(s) across {len(gtex_samples)} samples...")

    # Only the novel junctions' rows are ever kept beyond this loop --
    # build the padded-PSI column per sample, then immediately slice out
    # just the tail (novel) rows into the small matrix that gets fit.
    novel_psi_matrix = np.empty((n_novel_rows, len(gtex_samples)), dtype=np.float64)

    for i, sample in enumerate(gtex_samples):
        padded_counts = np.concatenate([
            gtex_df[sample + '_jxn_alignment_count'].to_numpy(),
            np.zeros(n_novel_rows, dtype=np.float64),
        ])
        padded_psi = compute_rescaled_psi_column(
            padded_counts, padded_ss1, padded_ss2,
            args.PSI_rescale_factor, args.gtex_coverage_threshold
        )
        novel_psi_matrix[:, i] = padded_psi[n_gtex_rows:]

    del gtex_df

    results = []
    for junction, row in zip(novel_junctions, novel_psi_matrix):
        results.append((junction,) + fit_beta_dist(row, args.PSI_rescale_factor, args.gtex_n_threshold))

    out_df = pd.DataFrame(
        results,
        columns=['junction', 'num_gtex_samples_with_good_coverage', 'alpha', 'beta',
                 'expected_PSI', 'p1_PSI', 'p99_PSI']
    )
    out_df['in_gtex_matrix'] = False
    out_df = out_df.set_index('junction').sort_index()
    out_df.to_csv(args.outfile, sep='\t')

    print(f"\nWrote beta fits for {len(out_df)} novel junction(s) to {args.outfile}")
    print("\nFinished fitting novel junction beta distributions.\n")


if __name__ == '__main__':
    main()
