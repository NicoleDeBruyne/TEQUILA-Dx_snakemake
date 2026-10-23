import argparse
import concurrent.futures

import pandas as pd
import numpy as np


def read_bed(path: str) -> pd.DataFrame:
    bed = pd.read_csv(
        path, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2], names=["chrom", "start", "end"]
    )
    bed["start"] = bed["start"].astype(int)
    bed["end"]   = bed["end"].astype(int)
    return bed


def build_bed_dict(bed: pd.DataFrame) -> dict:
    bed_dict = {}
    for chrom, sub in bed.groupby("chrom", sort=False):
        starts = np.sort(sub["start"].to_numpy())
        ends   = sub["end"].to_numpy()[np.argsort(sub["start"].to_numpy())]
        bed_dict[chrom] = (starts, ends)
    return bed_dict


def _in_bed(bed_dict: dict, chrom_arr: np.ndarray, pos_arr: np.ndarray) -> np.ndarray:
    result = np.zeros(len(pos_arr), dtype=bool)
    for chrom, (starts, ends) in bed_dict.items():
        mask = (chrom_arr == chrom)
        if not mask.any():
            continue
        p = pos_arr[mask]
        idx = np.searchsorted(starts, p, side="right") - 1
        valid = idx >= 0
        hit = np.zeros(len(p), dtype=bool)
        if valid.any():
            hit[valid] = p[valid] < ends[idx[valid]]
        result[mask] = hit
    return result


def touches_bed(index: pd.Index, bed_dict: dict) -> np.ndarray:
    keep = np.zeros(len(index), dtype=bool)
    split = index.to_series(index=range(len(index))).astype(str).str.rsplit("_", n=2, expand=True)
    if split.shape[1] < 3:
        return keep
    ss1 = pd.to_numeric(split[1], errors="coerce")
    ss2 = pd.to_numeric(split[2], errors="coerce")
    ok = (split[0].notna() & ss1.notna() & ss2.notna()).to_numpy()
    chrom_arr = split[0].to_numpy()[ok]
    hit = (_in_bed(bed_dict, chrom_arr, ss1.to_numpy()[ok].astype(np.int64))
           | _in_bed(bed_dict, chrom_arr, ss2.to_numpy()[ok].astype(np.int64)))
    keep[np.flatnonzero(ok)[hit]] = True
    return keep


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
        col = col[touches_bed(col.index, bed_dict)]
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
