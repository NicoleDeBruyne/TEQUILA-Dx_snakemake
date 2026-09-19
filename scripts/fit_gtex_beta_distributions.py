#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2026.09.17

"""
scripts/fit_gtex_beta_distributions.py

Fits a beta distribution to every junction present in one tissue's raw GTEx
splice-junction count matrix (gtex_{tissue}_jxn_counts.txt), genome-wide --
NOT scoped to any sample's gene panel or region.

This is rule _5B1_fit_gtex_beta_distributions' script. It exists because the
beta fit only ever depends on the GTEx reference data plus three threshold
config values (gtex_coverage_threshold, gtex_n_threshold, PSI_rescale_factor)
-- never on any individual sample's own read counts -- so refitting it once
per SAMPLE (as the pipeline used to, inside
perform_splice_junction_beta_binomial_tests.py) recomputed the exact same
answer, from scratch, for every sample that queried a given tissue. This
script computes it once per (tissue, threshold combination) and writes a
lookup table that every sample's rule _5C job reads instead of re-fitting.

Junctions present in a SAMPLE but absent from the raw GTEx matrix entirely
("novel" junctions) are NOT covered by this table -- see
fit_novel_junction_beta_distributions.py (rule _5B2), which handles exactly
that per-sample, per-tissue gap, since the set of novel junctions can't be
known ahead of any sample's own data.

Fitting (scipy.stats.beta.fit()'s iterative MLE, called once per junction)
is the actual bottleneck here, not the coverage/PSI computation -- a
genome-wide GTEx matrix can have far more junctions than any single BED
gene panel will ever query. So this script only FITS junctions overlapping
at least one BED panel used anywhere in the run (the union across all
samples -- see --bed-files); coverage/PSI are still computed genome-wide
first, exactly as before, so a kept junction's numbers are identical to
what a full genome-wide fit would have produced. Every sample's own
splice-junction counts (rule _5A) are already restricted to that sample's
own BED panel, so the union of all samples' panels is guaranteed to cover
every junction any _5C job in this run could look up.

Output columns (index = junction, "chrom_start_end"):
    num_gtex_samples_with_good_coverage, alpha, beta, expected_PSI, p1_PSI, p99_PSI, in_gtex_matrix
"in_gtex_matrix" is always True here (every row in this table came from the
raw GTEx matrix) -- it exists so the downstream script (_5C) can tell this
table's rows apart from fit_novel_junction_beta_distributions.py's rows
after concatenating both, without re-deriving membership itself.
"""

import argparse
import concurrent.futures
import os

import numpy as np
import pandas as pd
from math import ceil


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fit a beta distribution to every junction in a GTEx tissue's raw junction "
                    "count matrix, genome-wide, once -- reused across every sample that queries "
                    "this tissue with the same thresholds.")
    parser.add_argument('--gtexfile', required=True,
        help="Path to the raw per-GTEx-sample splice junction count matrix for one tissue.")
    parser.add_argument('--bed-files', required=True, nargs='+',
        help="One or more BED panel files (the union across every sample in the run). Only "
             "junctions overlapping at least one of these get the expensive beta.fit() call; "
             "coverage/PSI are still computed genome-wide first, so this only skips fitting "
             "junctions no sample in this run could ever query.")
    parser.add_argument('--outfile', required=True, type=str)
    parser.add_argument('--gtex-coverage-threshold', type=int, default=20)
    parser.add_argument('--PSI-rescale-factor', type=float, default=1e-3)
    parser.add_argument('--gtex-n-threshold', type=int, default=100)
    parser.add_argument('--threads', type=int, default=1)
    return parser.parse_args()


def compute_rescaled_psi_column(jxn_counts, ss1, ss2, PSI_rescale_factor, gtex_coverage_threshold):
    """Compute one GTEx sample's rescaled PSI value for every junction, as a
    plain float64 numpy array -- nothing is written back onto a shared
    dataframe.

    This replaces calling calculate_coverage()+calculate_PSI() (which used to
    add 5 new columns -- ss1_coverage, ss2_coverage, jxn_coverage, sample_PSI,
    rescaled_sample_PSI -- to the shared genome-wide dataframe, for EVERY
    GTEx sample) with a function that returns only the one column actually
    needed for fitting (rescaled_sample_PSI), as a compact float64 array
    instead of pandas Series (and never as the object-dtype array the
    original produced by mixing floats with the "n/a" string sentinel).
    Genome-wide, with hundreds of GTEx samples, that column accumulation was
    the direct cause of _5B1 OOMing in production -- see the conversation
    that led to this rewrite. The arithmetic itself (coverage via ss1/ss2
    groupby-sums, PSI = jxn/coverage, rescale, then NaN below the coverage
    threshold) is unchanged from calculate_coverage/calculate_PSI.

    Returns a float64 numpy array (NaN where PSI is undefined or coverage
    is below threshold) -- never the string "n/a"; that sentinel now only
    ever appears in the final CSV, written by pandas' NaN handling on
    output, not carried through the computation.
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


def load_merged_bed_regions(bed_paths):
    """Load one or more BED files and merge their intervals per chromosome into a sorted,
    non-overlapping union -- standard "sort, then merge touching/overlapping intervals" sweep.
    Only the first 3 columns (chrom, start, end) are used; BED coordinates are treated as plain
    integers, since the overlap test below (see junctions_overlap_regions) uses a closed-interval
    comparison that's deliberately a little permissive about the BED half-open-vs-junction
    1-based conventions, rather than risk excluding a junction at a panel's edge over an
    off-by-one.

    Returns {chrom: [(start, end), ...]} with each chromosome's list sorted by start and merged.
    """
    by_chrom = {}
    for path in bed_paths:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                chrom = parts[0].strip()
                start, end = int(parts[1]), int(parts[2])
                by_chrom.setdefault(chrom, []).append((start, end))

    merged = {}
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        out = []
        for start, end in intervals:
            if out and start <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged[chrom] = out
    return merged


def junctions_overlap_regions(chrom, start, end, merged_regions):
    """Vectorized test of whether each junction (chrom[i], start[i], end[i]) overlaps ANY
    interval in merged_regions (see load_merged_bed_regions), on the matching chromosome.

    Per chromosome, uses np.searchsorted against the merged intervals' start coordinates:
    since the intervals are sorted and non-overlapping, the only interval that could possibly
    overlap a query is the one immediately at-or-before the query's end -- if that one doesn't
    overlap, no earlier interval (which starts and ends even earlier) can either. This keeps the
    cost at roughly O(junctions * log(panel intervals)) rather than O(junctions * panel
    intervals), which matters since a genome-wide junction count can be in the hundreds of
    thousands while a gene panel is usually much smaller.

    Returns a boolean numpy array, one entry per junction.
    """
    chrom = np.asarray(chrom)
    start = np.asarray(start, dtype=np.int64)
    end = np.asarray(end, dtype=np.int64)
    keep = np.zeros(len(chrom), dtype=bool)

    for c, intervals in merged_regions.items():
        mask = chrom == c
        if not mask.any():
            continue
        interval_starts = np.array([iv[0] for iv in intervals], dtype=np.int64)
        interval_ends = np.array([iv[1] for iv in intervals], dtype=np.int64)

        idx = np.where(mask)[0]
        js, je = start[idx], end[idx]

        cand = np.searchsorted(interval_starts, je, side='right') - 1
        valid = cand >= 0
        overlap = np.zeros(len(idx), dtype=bool)
        overlap[valid] = js[valid] <= interval_ends[cand[valid]]
        keep[idx[valid]] = overlap[valid]

    return keep


def fit_beta_dist(x, tol, n_threshold):
    """Identical to perform_splice_junction_beta_binomial_tests.py's function of the same name.
    See that module's docstring for the p1/p99 convention."""
    from scipy.stats import beta

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


def _fit_chunk(rows, tol, n_threshold):
    """Fit a chunk of (junction, psi_row_values) pairs. Runs in a worker
    process -- avoids pandas .apply(axis=1) overhead by iterating over plain
    numpy rows instead of constructing a pd.Series per call."""
    out = []
    for junction, row in rows:
        out.append((junction,) + fit_beta_dist(row, tol, n_threshold))
    return out


def main():
    print("\n\n\n******************************************************************************************")
    print("Fitting beta distributions on GTEx reference junctions (genome-wide, once per tissue)...")
    print("******************************************************************************************\n")

    args = parse_args()

    if not os.path.exists(args.gtexfile):
        print(f"\nERROR: {args.gtexfile} not found.")
        return

    os.makedirs(os.path.dirname(args.outfile) or ".", exist_ok=True)

    print(f"\nReading GTEx file {args.gtexfile}...")
    gtex_df = pd.read_csv(args.gtexfile, sep='\t', index_col=0)

    gtex_samples = list(gtex_df.columns)
    gtex_df.columns = [c + '_jxn_alignment_count' for c in gtex_samples]

    # Normalize the GTEx index ("chr:start-end:strand") to this pipeline's
    # "chr_start_end" convention -- same normalization
    # perform_splice_junction_beta_binomial_tests.py used to do per-region.
    idx_chrom_coord = gtex_df.index.str.split(':')
    gtex_chrom = idx_chrom_coord.map(lambda p: p[0])
    idx_coords = idx_chrom_coord.map(lambda p: p[1]).str.split('-')
    gtex_start = idx_coords.map(lambda p: int(p[0]))
    gtex_end = idx_coords.map(lambda p: int(p[1]))

    gtex_df.index = gtex_chrom.astype(str) + '_' + gtex_start.astype(str) + '_' + gtex_end.astype(str)
    ss1 = (gtex_chrom.values + '_' + (gtex_start.values).astype(str))
    ss2 = (gtex_chrom.values + '_' + (gtex_end.values).astype(str))
    junctions = gtex_df.index.to_numpy()
    junction_chrom = gtex_chrom.to_numpy()
    junction_start = gtex_start.to_numpy()
    junction_end = gtex_end.to_numpy()

    print(f"Computing coverage and PSI for {len(gtex_df)} GTEx junctions across "
          f"{len(gtex_samples)} samples...")

    # Build the rescaled-PSI matrix one GTEx sample (column) at a time,
    # directly into a pre-allocated float64 array -- never accumulating
    # per-sample coverage/PSI columns on gtex_df itself (see
    # compute_rescaled_psi_column()'s docstring for why that mattered).
    psi_matrix = np.empty((len(gtex_df), len(gtex_samples)), dtype=np.float64)

    for i, sample in enumerate(gtex_samples):
        psi_matrix[:, i] = compute_rescaled_psi_column(
            gtex_df[sample + '_jxn_alignment_count'].to_numpy(),
            ss1, ss2,
            args.PSI_rescale_factor, args.gtex_coverage_threshold
        )

    del gtex_df

    print(f"\nLoading BED panel(s) to restrict fitting to junctions any sample in this run "
          f"could query: {args.bed_files}")
    merged_regions = load_merged_bed_regions(args.bed_files)
    in_panel = junctions_overlap_regions(junction_chrom, junction_start, junction_end, merged_regions)

    print(f"{in_panel.sum()} of {len(junctions)} GTEx junctions overlap the run's BED panel(s) "
          f"-- only these will be fit; the rest are skipped (they can never be queried by any "
          f"sample in this run).")

    junctions = junctions[in_panel]
    psi_matrix = psi_matrix[in_panel]

    print(f"\nFitting beta distributions for {len(junctions)} junctions using {args.threads} thread(s). "
          f"This may take a while...")

    rows = list(zip(junctions, psi_matrix))

    if args.threads <= 1 or len(rows) < 2:
        results = _fit_chunk(rows, args.PSI_rescale_factor, args.gtex_n_threshold)
    else:
        chunk_size = max(1, ceil(len(rows) / args.threads))
        chunks = [rows[i:i + chunk_size] for i in range(0, len(rows), chunk_size)]
        results = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = [
                executor.submit(_fit_chunk, chunk, args.PSI_rescale_factor, args.gtex_n_threshold)
                for chunk in chunks
            ]
            for future in concurrent.futures.as_completed(futures):
                results.extend(future.result())

    out_df = pd.DataFrame(
        results,
        columns=['junction', 'num_gtex_samples_with_good_coverage', 'alpha', 'beta',
                 'expected_PSI', 'p1_PSI', 'p99_PSI']
    )
    out_df['in_gtex_matrix'] = True
    out_df = out_df.set_index('junction').sort_index()

    with open(args.outfile, 'w') as f:
        f.write(
            f"# gtex_coverage_threshold={args.gtex_coverage_threshold}\t"
            f"gtex_n_threshold={args.gtex_n_threshold}\t"
            f"PSI_rescale_factor={args.PSI_rescale_factor}\n"
        )
    out_df.to_csv(args.outfile, sep='\t', mode='a')

    print(f"\nWrote beta fits for {len(out_df)} GTEx junctions to {args.outfile}")
    print("\nFinished fitting GTEx beta distributions.\n")


if __name__ == '__main__':
    main()
