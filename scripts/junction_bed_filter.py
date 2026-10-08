"""
BED helpers shared by build_group_junction_matrix.py (_9C1) and validate_sample_type.py (_9C2):
keep only junctions with at least one splice site inside the panel BED.
"""
import numpy as np
import pandas as pd


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
    """True for junctions with at least one splice site inside the BED. PSI of a junction with
    both sites in the BED only uses junctions sharing one of its sites, all of which pass this
    test, so dropping the rest first leaves every BED junction's PSI unchanged."""
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
