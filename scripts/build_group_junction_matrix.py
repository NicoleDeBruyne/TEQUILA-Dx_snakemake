#!/usr/bin/env python3

# Combines per-sample junction count matrices (each: columns = [junction, <sample>])
# into a single matrix with one column per sample, for use as one --matrix-query
# entry in validate_sample_type.py (which expects one matrix per query group, with
# a column per sample belonging to that group).

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-sample junction count matrices into one group-level matrix.")
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Per-sample junction count matrix TSVs (columns: junction, <sample>).")
    parser.add_argument("--sample-names", nargs="+", required=True,
        help="Clean sample name for each --infiles entry (same order), used as the output column name.")
    parser.add_argument("--outfile", required=True)
    return parser.parse_args()


def _junction_index(path):
    """Read just the first (junction ID) column of a per-sample matrix, as an
    Index -- avoids loading the count column during the union pass."""
    col = pd.read_csv(path, sep="\t", usecols=[0])
    return pd.Index(col.iloc[:, 0])


def main():
    args = parse_args()
    if len(args.infiles) != len(args.sample_names):
        raise ValueError("--infiles and --sample-names must have the same number of entries")

    # Pass 1: build the full union of junction IDs across all samples without
    # ever holding more than one file's junction column in memory at a time.
    # This avoids the old approach's repeated whole-DataFrame copies (every
    # incremental `combined.join(df, how="outer")` allocated a brand-new,
    # ever-growing DataFrame).
    union_index = None
    for path in args.infiles:
        idx = _junction_index(path)
        union_index = idx if union_index is None else union_index.union(idx)
    union_index.name = "junction"

    # Pass 2: preallocate the final matrix at its true final size, then fill
    # in one sample's column at a time, reading (and discarding) one file at
    # a time -- same one-file-in-memory-at-a-time discipline as before, but
    # without re-copying the growing combined matrix on every iteration.
    combined = pd.DataFrame(0.0, index=union_index, columns=args.sample_names)
    for path, sample in zip(args.infiles, args.sample_names):
        df = pd.read_csv(path, sep="\t")
        df = df.set_index(df.columns[0])   # first column is the junction ID, whatever its header is
        combined.loc[df.index, sample] = df.iloc[:, 0].values

    combined.to_csv(args.outfile, sep="\t")
    print(f"Wrote combined junction count matrix for {len(args.sample_names)} sample(s) to {args.outfile}")


if __name__ == "__main__":
    main()
