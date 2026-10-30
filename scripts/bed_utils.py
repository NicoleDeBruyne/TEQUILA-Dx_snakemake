# BED helpers used by a few of the junction scripts (were copy-pasted in each one before)

from typing import Dict, Tuple

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


# sorted starts/ends per chrom so overlap checks can use searchsorted
def build_bed_dict(bed: pd.DataFrame) -> dict:
    bed_dict = {}
    for chrom, sub in bed.groupby("chrom", sort=False):
        starts = np.sort(sub["start"].to_numpy())
        ends   = sub["end"].to_numpy()[np.argsort(sub["start"].to_numpy())]
        bed_dict[chrom] = (starts, ends)
    return bed_dict


def in_bed(bed_dict: dict, chrom_arr: np.ndarray, pos_arr: np.ndarray) -> np.ndarray:
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


# a junction counts if either end lands in the BED
def touches_bed(index: pd.Index, bed_dict: dict) -> np.ndarray:
    keep = np.zeros(len(index), dtype=bool)
    split = index.to_series(index=range(len(index))).astype(str).str.rsplit("_", n=2, expand=True)
    if split.shape[1] < 3:
        return keep
    ss1 = pd.to_numeric(split[1], errors="coerce")
    ss2 = pd.to_numeric(split[2], errors="coerce")
    ok = (split[0].notna() & ss1.notna() & ss2.notna()).to_numpy()
    chrom_arr = split[0].to_numpy()[ok]
    hit = (in_bed(bed_dict, chrom_arr, ss1.to_numpy()[ok].astype(np.int64))
           | in_bed(bed_dict, chrom_arr, ss2.to_numpy()[ok].astype(np.int64)))
    keep[np.flatnonzero(ok)[hit]] = True
    return keep


# gene -> (chrom, start, end), from the 4th BED column
def load_bed(path: str) -> Dict[str, Tuple[str, str, str]]:
    gene_info: Dict[str, Tuple[str, str, str]] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                print(f"[WARNING] BED line has <6 columns, skipping: {line!r}")
                continue
            chrom, start, end, gene, _, strand = parts[:6]
            gene   = gene.strip()
            strand = strand.strip()
            region = f"{chrom.strip()}:{start.strip()}-{end.strip()}"
            if gene in gene_info:
                print(f"[WARNING] Gene '{gene}' appears more than once in BED. Using last entry.")
            gene_info[gene] = (chrom.strip(), region, strand)
    print(f"BED file loaded: {len(gene_info)} gene(s)")
    return gene_info
