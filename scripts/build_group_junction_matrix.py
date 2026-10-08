

import argparse
import concurrent.futures

import pandas as pd

from junction_bed_filter import read_bed, build_bed_dict, touches_bed as _touches_bed


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-sample junction count matrices into one group-level matrix.")
    parser.add_argument("--infiles", nargs="+", required=True,
        help="Per-sample junction count matrix TSVs (columns: junction, <sample>).")
    parser.add_argument("--sample-names", nargs="+", required=True,
        help="Clean sample name for each --infiles entry (same order), used as the output column name.")
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--bed", default=None,
        help="Panel BED. If given, only junctions with at least one splice site inside it are kept "
             "(all that the sample-type validation uses), instead of every junction genome-wide.")
    parser.add_argument("--threads", type=int, default=1,
        help="Number of files to read in parallel (I/O-bound; thread-based). Default: 1")
    return parser.parse_args()


def _read_sample_column(path, sample, bed_dict=None):
    df = pd.read_csv(path, sep="\t")
    df = df.set_index(df.columns[0])
    col = df.iloc[:, 0]
    if bed_dict is not None:
        col = col[_touches_bed(col.index, bed_dict)]
    return sample, col


def main():
    args = parse_args()
    if len(args.infiles) != len(args.sample_names):
        raise ValueError("--infiles and --sample-names must have the same number of entries")

    bed_dict = build_bed_dict(read_bed(args.bed)) if args.bed else None
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
        cols = dict(ex.map(lambda ps: _read_sample_column(ps[0], ps[1], bed_dict),
                           zip(args.infiles, args.sample_names)))

    union_index = None
    for sample in args.sample_names:
        idx = pd.Index(cols[sample].index)
        union_index = idx if union_index is None else union_index.union(idx)
    union_index.name = "junction"

    combined = pd.DataFrame(0.0, index=union_index, columns=args.sample_names)
    for sample in args.sample_names:
        col = cols.pop(sample)
        combined.loc[col.index, sample] = col.values

    combined.to_csv(args.outfile, sep="\t")
    print(f"Wrote combined junction count matrix for {len(args.sample_names)} sample(s) to {args.outfile}")


if __name__ == "__main__":
    main()
