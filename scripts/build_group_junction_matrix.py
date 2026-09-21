

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
    col = pd.read_csv(path, sep="\t", usecols=[0])
    return pd.Index(col.iloc[:, 0])


def _read_sample_column(path, sample):
    df = pd.read_csv(path, sep="\t")
    df = df.set_index(df.columns[0])
    return sample, df.iloc[:, 0]


def main():
    args = parse_args()
    if len(args.infiles) != len(args.sample_names):
        raise ValueError("--infiles and --sample-names must have the same number of entries")

    union_index = None
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
        for idx in ex.map(_junction_index, args.infiles):
            union_index = idx if union_index is None else union_index.union(idx)
    union_index.name = "junction"

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
