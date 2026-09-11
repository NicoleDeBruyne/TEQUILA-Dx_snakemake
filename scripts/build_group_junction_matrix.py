#!/usr/bin/env python3

# Combines per-sample junction count matrices (each: columns = [junction, <sample>])
# into a single matrix with one column per sample, for use as one --matrix-query
# entry in validate_sample_type.py (which expects one matrix per query group, with
# a column per sample belonging to that group).
#
# Reading is the whole cost here (each file is small, but there can be
# hundreds of samples in a group) -- purely I/O-bound, no CPU-heavy work per
# file, so this uses a ThreadPoolExecutor (not processes): pandas.read_csv's
# C parser releases the GIL for the bulk of its work, and threads avoid the
# cost of pickling/copying file paths and small DataFrames across process
# boundaries for what's fundamentally a "wait on disk" workload. The shared
# `combined` matrix is only ever written from the main thread, after every
# worker's read has completed -- workers return data, they never touch
# `combined` themselves, so there's no cross-thread write race.

import argparse
import concurrent.futures

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-sample junction count matrices into one group-level matrix.")
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Per-sample junction count matrix TSVs (columns: junction, <sample>).")
    parser.add_argument("--sample-names", nargs="+", required=True,
        help="Clean sample name for each --infiles entry (same order), used as the output column name.")
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--threads", type=int, default=1,
        help="Number of files to read in parallel (I/O-bound; thread-based). Default: 1")
    return parser.parse_args()


def _junction_index(path):
    """Read just the first (junction ID) column of a per-sample matrix, as an
    Index -- avoids loading the count column during the union pass."""
    col = pd.read_csv(path, sep="\t", usecols=[0])
    return pd.Index(col.iloc[:, 0])


def _read_sample_column(path, sample):
    """Read one sample's full matrix, indexed by junction ID -- returns
    (sample, indexed_series) so the caller can fill `combined` itself."""
    df = pd.read_csv(path, sep="\t")
    df = df.set_index(df.columns[0])   # first column is the junction ID, whatever its header is
    return sample, df.iloc[:, 0]


def main():
    args = parse_args()
    if len(args.infiles) != len(args.sample_names):
        raise ValueError("--infiles and --sample-names must have the same number of entries")

    # Pass 1: build the full union of junction IDs across all samples without
    # ever holding more than one file's junction column in memory at a time
    # per thread. Reads happen in parallel; the union reduction itself stays
    # single-threaded (it's cheap relative to the file reads, and doing it
    # incrementally in the main thread as results arrive avoids needing a
    # lock around a shared Index).
    union_index = None
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
        for idx in ex.map(_junction_index, args.infiles):
            union_index = idx if union_index is None else union_index.union(idx)
    union_index.name = "junction"

    # Pass 2: read every sample's full column in parallel, then fill the
    # preallocated final matrix one column at a time in the main thread --
    # `combined` is never touched from inside a worker, so concurrent reads
    # can't race on it.
    combined = pd.DataFrame(0.0, index=union_index, columns=args.sample_names)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
        futures = [ex.submit(_read_sample_column, path, sample)
                   for path, sample in zip(args.infiles, args.sample_names)]
        for f in concurrent.futures.as_completed(futures):
            sample, col = f.result()
            combined.loc[col.index, sample] = col.values

    combined.to_csv(args.outfile, sep="\t")
    print(f"Wrote combined junction count matrix for {len(args.sample_names)} sample(s) to {args.outfile}")


if __name__ == "__main__":
    main()
