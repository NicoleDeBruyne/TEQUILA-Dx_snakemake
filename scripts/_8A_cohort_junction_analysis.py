# 8A - per-gene junction tables for a sample-type group. junction reads plus splice-site, intron retention
# and intronic polyA (IPA) coverage, for bulk and both haplotypes. one TSV per gene

from __future__ import annotations

import os
import argparse
import warnings
import traceback
import time
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
import concurrent.futures
from pandas.errors import PerformanceWarning
from bed_utils import load_bed
from gtf_utils import parse_gtf_junctions

warnings.filterwarnings("ignore", category=PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


_REF_OPS  = frozenset((0, 2, 3, 7, 8))
_SKIP_OP  = 3

# what each ratio gets divided by
_METRIC_DENOMINATOR: Dict[str, str] = {
    "junction_PSI_approx":    "junction_coverage_approx",
    "junction_PSI":           "junction_coverage",
    "5ss_IR_ratio":           "5ss_coverage",
    "3ss_IR_ratio":           "3ss_coverage",
    "junction_full_IR_ratio": "junction_coverage",
    "junction_IPA_ratio":     "5ss_coverage",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Computes per-gene, per-sample splice junction coverage/usage "
                     "metrics for cohort-level outlier analysis (core analysis stage)."
    )
    p.add_argument("--mapping-file",               required=True)
    p.add_argument("--bed",                        required=True)
    p.add_argument("--outdir",                     required=True,
                   help="Directory to write one raw metrics TSV per gene into.")
    p.add_argument("--manifest",                    required=True,
                   help="Path to write the gene -> result-file manifest TSV to. "
                        "Contains one row per gene in --bed.")
    p.add_argument("--note",                        required=True,
                   help="Path to write a short human-readable note to, explaining "
                        "whether the analysis ran or was skipped (and why).")
    p.add_argument("--min-samples",                 type=int,   default=10,
                   help="Minimum number of unique samples in --mapping-file required "
                        "to run the analysis at all. Groups below this are skipped "
                        "entirely (manifest is written with every gene set to 'None', "
                        "and --note explains why) -- fitting a per-junction cohort "
                        "distribution from a handful of samples isn't meaningful.")
    p.add_argument("--approx",                     action="store_true")
    p.add_argument("--PSI-rescale-factor",         type=float, default=1e-3)
    p.add_argument("--min-jxn-reads",              type=int,   default=20)
    p.add_argument("--include-monoexonic",         action="store_true")
    p.add_argument("--genome",                     default=None)
    p.add_argument("--alu-bed",                    default=None)
    p.add_argument("--gtf",                        default=None,
                   help="GTF/GTF.gz. If provided, adds a junction_type column (canonical / annotated / novel).")
    p.add_argument("--threads",                    type=int,   default=1)
    p.add_argument("--test-n-genes",               type=int,   default=None)
    return p.parse_args()


def load_and_validate_mapping(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()
    _CANON = ["gene", "sample", "bulk", "hap1", "hap2"]
    rename: Dict[str, str] = {}
    used: set = set()
    for col in df.columns:
        cl = col.lower().replace("-", "").replace("_", "")
        for canonical in _CANON:
            if canonical in cl and canonical not in used:
                rename[col] = canonical
                used.add(canonical)
                break
    df = df.rename(columns=rename)
    missing_cols = [c for c in ("gene", "sample", "bulk") if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Could not identify required columns {missing_cols}")
    for col in ("hap1", "hap2"):
        if col not in df.columns:
            df[col] = np.nan
    df = df.apply(lambda col: col.str.strip() if col.dtype == object else col)
    df.replace({"NA": np.nan, "": np.nan, "na": np.nan, "N/A": np.nan,
                "None": np.nan, "none": np.nan}, inplace=True)
    print(f"Mapping file loaded: {df['gene'].nunique()} gene(s), "
          f"{df['sample'].nunique()} sample(s), {len(df)} row(s)")
    return df


def _parse_region(region: str) -> Tuple[str, int, int]:
    chrom, se = region.split(":")
    start, end = map(int, se.split("-"))
    return chrom, start, end


def collect_read_data(
    bam_path: str,
    region: str,
    include_monoexonic: bool = False,
    collect_softclips: bool = False,
    strand: str = "+",
    genome_seq: Optional[str] = None,
    gene_region_start: int = 0,
    keep_reads: bool = True,
) -> Tuple[Dict[Tuple[int, int], int], int,
           List[Tuple[List[Tuple[int, int]], List[int], Optional[Tuple]]],
           ]:
    reads: List = []
    counts: Dict[str, object] = {}
    for read in iter_read_data(bam_path, region, include_monoexonic, collect_softclips, strand,
                               genome_seq, gene_region_start, counts):
        if keep_reads:
            reads.append(read)
    return counts["jxn_raw"], counts["gene_cov"], reads


def iter_read_data(
    bam_path: str,
    region: str,
    include_monoexonic: bool = False,
    collect_softclips: bool = False,
    strand: str = "+",
    genome_seq: Optional[str] = None,
    gene_region_start: int = 0,
    counts: Optional[Dict[str, object]] = None,
):
    counts = {} if counts is None else counts
    chrom_r, region_start, region_end = _parse_region(region)
    jxn_raw:  Dict[Tuple[int, int], int] = defaultdict(int)
    gene_cov  = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(region=region):
            if aln.is_secondary or aln.cigartuples is None:
                continue
            gene_cov   += 1
            pos          = aln.reference_start
            has_splice   = False
            has_softclip = False
            blocks:      List[Tuple[int, int]] = []
            block_start: Optional[int]         = None
            read_jxns:   List[Tuple[int, int]] = []

            for op, length in aln.cigartuples:
                if op == _SKIP_OP:
                    has_splice = True
                    if block_start is not None:
                        blocks.append((block_start, pos))
                        block_start = None
                    j_start = pos + 1
                    j_end   = pos + length
                    if j_start < region_end and j_end > region_start:
                        jxn_raw[(j_start, j_end)] += 1
                        read_jxns.append((j_start, j_end))
                    pos += length
                elif op == 4:
                    has_softclip = True
                elif op in _REF_OPS:
                    if block_start is None:
                        block_start = pos
                    pos += length

            if block_start is not None:
                blocks.append((block_start, pos))

            if not has_splice:
                if not include_monoexonic:
                    gene_cov -= 1
                    continue
                yield (blocks, [], None)
                continue

            sc_tuple = None
            if collect_softclips and has_softclip and genome_seq is not None:
                qs   = aln.query_sequence
                qa_s = aln.query_alignment_start
                qa_e = aln.query_alignment_end
                rs0  = aln.reference_start
                re0  = aln.reference_end
                if qs:
                    if strand == "+":
                        sc   = qs[qa_e:] if qa_e < aln.query_length else ""
                        g    = genome_seq[re0 - gene_region_start - 20:
                                          re0 - gene_region_start]
                        pos3 = re0 - 1
                        lead = False
                    else:
                        sc   = qs[:qa_s] if qa_s > 0 else ""
                        g    = genome_seq[rs0 - gene_region_start:
                                          rs0 - gene_region_start + 20]
                        pos3 = rs0
                        lead = True
                    sc_tuple = (sc, g, pos3, lead)

            yield (blocks, read_jxns, sc_tuple)

    counts["jxn_raw"] = jxn_raw
    counts["gene_cov"] = gene_cov


# approx mode: count unique read footprints instead of going through every read
def footprint_unique_counts(footprints, jxn_coords):
    index = {c: i for i, c in enumerate(jxn_coords)}
    by_start, by_end = defaultdict(list), defaultdict(list)
    for i, (a, b) in enumerate(jxn_coords):
        by_start[a].append(i)
        by_end[b].append(i)
    per_chain = defaultdict(int)
    for fp in footprints:
        if len(fp) > 1:
            per_chain[tuple((fp[k][1] + 1, fp[k + 1][0]) for k in range(len(fp) - 1))] += 1
    usage = np.zeros(len(jxn_coords), dtype=np.int64)
    approx = np.zeros(len(jxn_coords), dtype=np.int64)
    for chain, n_fp in per_chain.items():
        used = {index[c] for c in chain if c in index}
        if used:
            usage[list(used)] += n_fp
        sites = set()
        for a, b in chain:
            sites.update(by_start.get(a, ()))
            sites.update(by_end.get(b, ()))
        if sites:
            approx[list(sites)] += n_fp
    return usage, approx


# splice-site / IR / IPA counting, done in batches of reads
class CoverageAccumulator:

    BATCH = 4096

    def __init__(self, jxn_coords_for_cov, strand, alu, chrom):
        self.jxn_coords = list(jxn_coords_for_cov)
        n = self.n = len(self.jxn_coords)
        self.strand, self.alu, self.chrom = strand, alu, chrom
        self.ss1 = np.array([a for a, b in self.jxn_coords], dtype=np.int64)
        self.ss2 = np.array([b for a, b in self.jxn_coords], dtype=np.int64)
        self.o1 = np.argsort(self.ss1, kind="stable"); self.s1 = self.ss1[self.o1]
        self.o2 = np.argsort(self.ss2, kind="stable"); self.s2 = self.ss2[self.o2]
        self.five_ss_arr = self.ss1 if strand == "+" else self.ss2
        z = lambda: np.zeros(n, dtype=np.int64)
        self.ss1_cov, self.ss2_cov, self.jxn_cov = z(), z(), z()
        self.ss1_ir, self.ss2_ir, self.full_ir, self.ipa = z(), z(), z(), z()
        self.u_ss1_cov, self.u_ss2_cov, self.u_jxn_cov = z(), z(), z()
        self.u_ss1_ir, self.u_ss2_ir, self.u_full_ir, self.u_ipa = z(), z(), z(), z()
        self._buf: List = []

    # a polyA-looking soft clip means oligo-dT priming. Alu polyA counts too
    def is_polya(self, sc_tuple) -> bool:
        sc, g, pos3, leading = sc_tuple
        if _is_oligo_dt_priming(sc, g, leading):
            return True
        if self.alu is not None:
            frag = sc[-10:] if leading else sc[:10]
            if len(frag) == 10 and (frag.count("A")/10 >= 0.8 or frag.count("T")/10 >= 0.8):
                if _in_alu(self.chrom, pos3, self.alu):
                    return True
        return False

    def add(self, blocks, splice_jxns, weight=1, n_polya=0):
        if self.n == 0 or not blocks:
            return
        self._buf.append((blocks, splice_jxns, weight, n_polya))
        if len(self._buf) >= self.BATCH:
            self._flush()

    @staticmethod
    def _pairs(order, sorted_vals, lo_vals, hi_vals):
        lo = np.searchsorted(sorted_vals, lo_vals, side="left")
        hi = np.searchsorted(sorted_vals, hi_vals, side="right")
        cnt = np.maximum(hi - lo, 0)
        total = int(cnt.sum())
        if total == 0:
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
        qi = np.repeat(np.arange(len(lo)), cnt)
        offs = np.arange(total) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        return qi, order[np.repeat(lo, cnt) + offs]

    def _add_counts(self, keys, w, target, utarget):
        if len(keys) == 0:
            return
        j = keys % self.n
        target += np.bincount(j, weights=w[keys // self.n], minlength=self.n).astype(np.int64)
        utarget += np.bincount(j, minlength=self.n)

    def _count(self, read_idx, jidx, w, target, utarget):
        if len(jidx) == 0:
            return np.empty(0, dtype=np.int64)
        keys = np.unique(read_idx * self.n + jidx)
        self._add_counts(keys, w, target, utarget)
        return keys

    def _flush(self):
        buf, self._buf = self._buf, []
        if not buf:
            return
        n = self.n
        w = np.fromiter((b[2] for b in buf), dtype=np.int64, count=len(buf))
        nb = np.fromiter((len(b[0]) for b in buf), dtype=np.int64, count=len(buf))
        block_read = np.repeat(np.arange(len(buf)), nb)
        flat = [blk for b in buf for blk in b[0]]
        bs = np.fromiter((x[0] for x in flat), dtype=np.int64, count=len(flat))
        be = np.fromiter((x[1] for x in flat), dtype=np.int64, count=len(flat))
        r_start = np.fromiter((b[0][0][0] for b in buf), dtype=np.int64, count=len(buf))
        r_end = np.fromiter((b[0][-1][1] for b in buf), dtype=np.int64, count=len(buf))
        ss1, ss2 = self.ss1, self.ss2

        qi, j = self._pairs(self.o1, self.s1, bs + 4, be + 1)
        r = block_read[qi]; keep = r_end[r] > ss1[j] + 1
        k1 = self._count(r[keep], j[keep], w, self.ss1_cov, self.u_ss1_cov)
        qi, j = self._pairs(self.o2, self.s2, bs, be - 3)
        r = block_read[qi]; keep = r_start[r] <= ss2[j] - 3
        k2 = self._count(r[keep], j[keep], w, self.ss2_cov, self.u_ss2_cov)
        if len(k1) or len(k2):
            self._add_counts(np.union1d(k1, k2), w, self.jxn_cov, self.u_jxn_cov)
        qi, j = self._pairs(self.o1, self.s1, bs + 4, be - 2)
        self._count(block_read[qi], j, w, self.ss1_ir, self.u_ss1_ir)
        qi, j = self._pairs(self.o2, self.s2, bs + 3, be - 3)
        self._count(block_read[qi], j, w, self.ss2_ir, self.u_ss2_ir)
        qi, j = self._pairs(self.o1, self.s1, bs + 4, be - 3)
        keep = ss2[j] <= be[qi] - 3
        self._count(block_read[qi][keep], j[keep], w, self.full_ir, self.u_full_ir)

        five_keys = k1 if self.strand == "+" else k2
        if not len(five_keys):
            return
        five_read = five_keys // n
        five_j = five_keys % n
        for ri, (blocks, splice_jxns, _weight, n_polya) in enumerate(buf):
            if not n_polya:
                continue
            a = np.searchsorted(five_read, ri, side="left")
            b = np.searchsorted(five_read, ri, side="right")
            if a == b:
                continue
            five = five_j[a:b]
            ref_end = blocks[-1][1]
            cand = five[(ss1[five] <= ref_end) & (ref_end <= ss2[five])]
            if not len(cand):
                continue
            splice_set = set(jj[0] if self.strand == "+" else jj[1] for jj in splice_jxns)
            for ji in np.sort(cand):
                five_ss = int(self.five_ss_arr[ji])
                if not any(s >= five_ss for s in splice_set):
                    self.ipa[ji] += n_polya
                    self.u_ipa[ji] += 1

    def result(self) -> Dict[Tuple[int, int], List[int]]:
        self._flush()
        out = {}
        for i, jc in enumerate(self.jxn_coords):
            out[jc] = [
                int(self.ss1_cov[i]), int(self.ss2_cov[i]), int(self.jxn_cov[i]),
                int(self.ss1_ir[i]),  int(self.ss2_ir[i]),  int(self.full_ir[i]),
                int(self.ipa[i]),
                int(self.u_ss1_cov[i]), int(self.u_ss2_cov[i]), int(self.u_jxn_cov[i]),
                int(self.u_ss1_ir[i]),  int(self.u_ss2_ir[i]),  int(self.u_full_ir[i]),
                int(self.u_ipa[i]),
            ]
        return out


def compute_coverage_metrics(
    reads:              List,
    jxn_coords_for_cov: List[Tuple[int, int]],
    strand:             str,
    genome_seq:         Optional[str],
    gene_region_start:  int,
    alu:                Optional[dict],
    chrom:              str,
) -> Dict[Tuple[int, int], List[int]]:
    if not jxn_coords_for_cov:
        return {}
    acc = CoverageAccumulator(jxn_coords_for_cov, strand, alu, chrom)
    for blocks, splice_jxns, sc_tuple in reads:
        acc.add(blocks, splice_jxns, 1, int(sc_tuple is not None and acc.is_polya(sc_tuple)))
    return acc.result()


# ignore it if the genome itself has an A/T run there, that's internal priming not a real end
def _is_oligo_dt_priming(sc: str, g: str, leading: bool) -> bool:
    frag = sc[-10:] if leading else sc[:10]
    if len(frag) < 10:
        return False
    if not (frag.count("A") / 10 >= 0.8 or frag.count("T") / 10 >= 0.8):
        return False
    f5  = g[:5]  if leading else g[-5:]
    f10 = g[:10] if leading else g[-10:]
    if f5 in ("AAAAA", "TTTTT") or f10.count("A") / 10 >= 0.8 or f10.count("T") / 10 >= 0.8:
        return False
    return True


def _in_alu(chrom: str, pos0: int, alu: dict) -> bool:
    if alu is None or chrom not in alu:
        return False
    ivs = alu[chrom]; lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if ivs[mid][0] <= pos0: lo = mid + 1
        else:                   hi = mid - 1
    return hi >= 0 and ivs[hi][0] <= pos0 < ivs[hi][1]


def load_alu_intervals(alu_bed: str) -> dict:
    import gzip
    open_fn = gzip.open if alu_bed.endswith(".gz") else open
    ivs: Dict[str, list] = defaultdict(list)
    with open_fn(alu_bed, "rt" if alu_bed.endswith(".gz") else "r") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip().split("\t")
            try:
                chrom = f[0]; s, e = int(f[1]), int(f[2])
            except (ValueError, IndexError):
                continue
            ivs[chrom].append((s, e))
    for c in ivs: ivs[c].sort()
    return dict(ivs)


# approx coverage = reads using either splice site
def compute_junction_coverage_approx(
    junction_counts: Dict[str, int],
    all_jxns: Optional[List[str]] = None,
) -> Dict[str, int]:
    ss1_usage: Dict[str, int] = defaultdict(int)
    ss2_usage: Dict[str, int] = defaultdict(int)
    parsed: Dict[str, Tuple[str, str, str]] = {}
    for jxn, count in junction_counts.items():
        parts = jxn.split("_")
        chrom = "_".join(parts[:-2])
        ss1, ss2 = parts[-2], parts[-1]
        parsed[jxn] = (chrom, ss1, ss2)
        ss1_usage[f"{chrom}_{ss1}"] += count
        ss2_usage[f"{chrom}_{ss2}"] += count

    target_jxns = all_jxns if all_jxns is not None else list(junction_counts.keys())
    result: Dict[str, int] = {}
    for jxn in target_jxns:
        if jxn in parsed:
            chrom, ss1, ss2 = parsed[jxn]
        else:
            parts = jxn.split("_")
            chrom = "_".join(parts[:-2])
            ss1, ss2 = parts[-2], parts[-1]
        own_count = junction_counts.get(jxn, 0)
        result[jxn] = (ss1_usage[f"{chrom}_{ss1}"]
                       + ss2_usage[f"{chrom}_{ss2}"]
                       - own_count)
    return result


def process_sample(
    sample_name: str,
    bulk_bam:    Optional[str],
    hap1_bam:    Optional[str],
    hap2_bam:    Optional[str],
    region: str,
    gene: str,
    all_jxns: List[str],
    approx_only: bool,
    PSI_rescale_factor: float,
    strand: str = "+",
    include_monoexonic: bool = False,
    genome_path: Optional[str] = None,
    alu: Optional[dict] = None,
    jxn_raw_bulk: Optional[Dict[Tuple[int, int], int]] = None,
) -> Tuple[pd.DataFrame, float]:
    t0 = time.time()

    for label, bam in (("bulk", bulk_bam), ("hap1", hap1_bam), ("hap2", hap2_bam)):
        if bam and isinstance(bam, str) and not os.path.exists(bam):
            print(f"[WARNING] {label} BAM not found for sample '{sample_name}': {bam}.")
            if label == "hap1":   hap1_bam = None
            elif label == "hap2": hap2_bam = None
            else:                 bulk_bam = None

    if bool(hap1_bam) ^ bool(hap2_bam):
        present = "hap1" if hap1_bam else "hap2"
        print(f"[WARNING] Sample '{sample_name}' has {present} but not the other. Reverting both to NA.")
        hap1_bam = None; hap2_bam = None

    chrom = region.split(":")[0]
    empty = pd.DataFrame()

    if approx_only:
        # approx mode: just junction counts + approximate coverage
        raw: List[Tuple[str, int, dict, set]] = []
        for label, bam in (("bulk", bulk_bam), ("hap1", hap1_bam), ("hap2", hap2_bam)):
            if not bam or not isinstance(bam, str): continue
            read_counts: Dict[str, object] = {}
            footprints = set()
            for blocks, splice_jxns, _sc in iter_read_data(bam, region, include_monoexonic, counts=read_counts):
                footprints.add(tuple(blocks))
            gene_cov = read_counts["gene_cov"]
            if gene_cov == 0: continue
            jxn_counts = {f"{chrom}_{a}_{b}": c for (a, b), c in read_counts["jxn_raw"].items()}
            raw.append((label, gene_cov, jxn_counts, footprints))
        if not raw: return empty, time.time() - t0

        if not all_jxns:
            uni: set = set()
            for _, _, jc, _fp in raw: uni.update(jc)
            all_jxns = sorted(uni)

        jxn_index = {j: i for i, j in enumerate(all_jxns)}
        n_jxns = len(all_jxns)
        coords = [(int(j.split("_")[-2]), int(j.split("_")[-1])) for j in all_jxns]
        chunks = []
        for label, gene_cov, jxn_counts, footprints in raw:
            u_unique, c_unique = footprint_unique_counts(footprints, coords)
            cov_dict  = compute_junction_coverage_approx(jxn_counts)
            u = np.zeros(n_jxns, dtype=np.int32)
            c = np.zeros(n_jxns, dtype=np.int32)
            for jxn, cnt in jxn_counts.items():
                if jxn in jxn_index: u[jxn_index[jxn]] = cnt
            for jxn, cnt in cov_dict.items():
                if jxn in jxn_index: c[jxn_index[jxn]] = cnt
            chunks.append((label, u, c, u_unique, c_unique))

        rows = []
        for label, u, c, u_unique, c_unique in chunks:
            cf = c.astype(np.float32); uf = u.astype(np.float32)
            with np.errstate(divide="ignore", invalid="ignore"):
                psi = np.where(cf > 0, uf / cf, np.nan)
            rsc = np.where(np.isfinite(psi), psi*(1-2*PSI_rescale_factor)+PSI_rescale_factor, np.nan)
            rows.append(pd.DataFrame({
                "junction":                     all_jxns,
                "junction_usage":               u,
                "junction_usage_unique":        u_unique,
                "junction_coverage_approx":     c,
                "junction_coverage_approx_unique": c_unique,
                "junction_PSI_approx":          psi,
                "rescaled_junction_PSI_approx": rsc,
                "phasing":                      label,
            }))
        df = pd.concat(rows, ignore_index=True)
        df["sample"] = sample_name; df["region"] = region; df["gene"] = gene
        return df, time.time() - t0

    jxn_index = {j: i for i, j in enumerate(all_jxns)}
    n_jxns = len(all_jxns)

    genome_seq: Optional[str] = None
    gene_start = int(region.split(":")[1].split("-")[0])
    gene_end   = int(region.split(":")[1].split("-")[1])
    if genome_path:
        try:
            fa = pysam.FastaFile(genome_path)
            genome_seq = fa.fetch(chrom, gene_start, gene_end)
            fa.close()
        except Exception:
            genome_seq = None

    do_ipa   = genome_seq is not None

    jxn_coords_for_cov: List[Tuple[int, int]] = []
    for jxn in all_jxns:
        parts = jxn.split("_")
        jxn_coords_for_cov.append((int(parts[-2]), int(parts[-1])))

    # full mode: stream every BAM through the accumulator
    chunks_wide: List = []
    for label, bam in (("bulk", bulk_bam), ("hap1", hap1_bam), ("hap2", hap2_bam)):
        if not bam or not isinstance(bam, str): continue
        acc = CoverageAccumulator(jxn_coords_for_cov, strand, alu, chrom)
        counts: Dict[str, object] = {}
        fp_ids: Dict[Tuple, int] = {}
        for blocks, splice_jxns, sc_tuple in iter_read_data(
                bam, region, include_monoexonic,
                collect_softclips=do_ipa, strand=strand,
                genome_seq=genome_seq, gene_region_start=gene_start, counts=counts):
            fpk = tuple(blocks)
            entry = fp_ids.get(fpk)
            if entry is None:
                entry = fp_ids[fpk] = [splice_jxns, 0, 0]
            entry[1] += 1
            if sc_tuple is not None and acc.is_polya(sc_tuple):
                entry[2] += 1
        gene_cov = counts["gene_cov"]
        if gene_cov == 0: continue
        jxn_raw_lbl = jxn_raw_bulk if (label == "bulk" and jxn_raw_bulk is not None) else counts["jxn_raw"]


        for fpk, (fp_jxns, n_reads, n_polya) in fp_ids.items():
            acc.add(list(fpk), fp_jxns, n_reads, n_polya)
        cov_metrics = acc.result() if jxn_coords_for_cov else {}
        usage_u, approx_u = footprint_unique_counts(fp_ids.keys(), jxn_coords_for_cov)
        del fp_ids

        ss1_usage_agg: Dict[int, int] = defaultdict(int)
        ss2_usage_agg: Dict[int, int] = defaultdict(int)
        for (j_start, j_end), cnt in jxn_raw_lbl.items():
            ss1_usage_agg[j_start] += cnt
            ss2_usage_agg[j_end]   += cnt

        cov_result: Dict[str, Dict] = {}
        for i, jxn in enumerate(all_jxns):
            ss1, ss2 = jxn_coords_for_cov[i]
            five_ss_pos  = ss1 if strand == "+" else ss2
            three_ss_pos = ss2 if strand == "+" else ss1
            m = cov_metrics.get((ss1, ss2), [0]*14)
            cov_result[jxn] = {
                "junction_coverage":       m[2],
                "ss1_usage":               ss1_usage_agg.get(ss1, 0),
                "ss1_coverage":            m[0],
                "ss1_ir":                  m[3],
                "ss2_usage":               ss2_usage_agg.get(ss2, 0),
                "ss2_coverage":            m[1],
                "ss2_ir":                  m[4],
                "five_ss":                 f"{chrom}_{five_ss_pos}",
                "three_ss":                f"{chrom}_{three_ss_pos}",
                "junction_full_ir":        m[5],
                "junction_ipa":            m[6],
                "ss1_coverage_unique":     m[7],
                "ss2_coverage_unique":     m[8],
                "junction_coverage_unique": m[9],
                "ss1_ir_unique":           m[10],
                "ss2_ir_unique":           m[11],
                "junction_full_ir_unique": m[12],
                "junction_ipa_unique":     m[13],
                "junction_usage_unique":   int(usage_u[i]) if len(usage_u) else 0,
                "junction_coverage_approx_unique": int(approx_u[i]) if len(approx_u) else 0,
            }

        jxn_raw_str = {f"{chrom}_{s}_{e}": c for (s, e), c in jxn_raw_lbl.items()}
        approx_cov  = compute_junction_coverage_approx(jxn_raw_str, all_jxns)
        chunks_wide.append((label, cov_result, jxn_raw_str, approx_cov))
    if not chunks_wide: return empty, time.time() - t0


    def _arr(label, cov_key):
        for lbl, cov, jrs, ac in chunks_wide:
            if lbl == label:
                return np.array([cov.get(j, {}).get(cov_key, 0) for j in all_jxns], dtype=np.int32)
        return np.zeros(n_jxns, dtype=np.int32)

    def _arr_approx(label):
        for lbl, _cov, _jrs, ac in chunks_wide:
            if lbl == label:
                return np.array([ac.get(j, 0) for j in all_jxns], dtype=np.int32)
        return np.zeros(n_jxns, dtype=np.int32)

    # ratios (PSI, IR, IPA) and their rescaled versions
    denom_arrays: Dict[str, Dict[str, np.ndarray]] = {}
    for denom_col in ("junction_coverage", "5ss_coverage", "3ss_coverage"):
        denom_arrays[denom_col] = {
            "bulk": _arr("bulk", denom_col),
            "hap1": _arr("hap1", denom_col),
            "hap2": _arr("hap2", denom_col),
        }
    denom_arrays["junction_coverage_approx"] = {
        "bulk": _arr_approx("bulk"),
        "hap1": _arr_approx("hap1"),
        "hap2": _arr_approx("hap2"),
    }

    def _rescale(arr_f):
        return np.where(np.isfinite(arr_f),
                        arr_f * (1 - 2*PSI_rescale_factor) + PSI_rescale_factor, np.nan)

    rows = []
    for label, cov_result, jxn_raw_str, approx_cov in chunks_wide:
        ju   = np.array([jxn_raw_str.get(j, 0) for j in all_jxns], dtype=np.int32)
        jc_approx = np.array([approx_cov.get(j, 0) for j in all_jxns], dtype=np.int32)
        jc   = np.array([cov_result.get(j, {}).get("junction_coverage", 0) for j in all_jxns], dtype=np.int32)
        s1u  = np.array([cov_result.get(j, {}).get("ss1_usage",    0) for j in all_jxns], dtype=np.int32)
        s1c  = np.array([cov_result.get(j, {}).get("ss1_coverage", 0) for j in all_jxns], dtype=np.int32)
        s1ir = np.array([cov_result.get(j, {}).get("ss1_ir",       0) for j in all_jxns], dtype=np.int32)
        s2u  = np.array([cov_result.get(j, {}).get("ss2_usage",    0) for j in all_jxns], dtype=np.int32)
        s2c  = np.array([cov_result.get(j, {}).get("ss2_coverage", 0) for j in all_jxns], dtype=np.int32)
        s2ir = np.array([cov_result.get(j, {}).get("ss2_ir",       0) for j in all_jxns], dtype=np.int32)
        uq = lambda key: np.array([cov_result.get(j, {}).get(key, 0) for j in all_jxns], dtype=np.int32)
        ju_u, jca_u, jc_u = uq("junction_usage_unique"), uq("junction_coverage_approx_unique"), uq("junction_coverage_unique")
        s1c_u, s2c_u, s1ir_u, s2ir_u = uq("ss1_coverage_unique"), uq("ss2_coverage_unique"), uq("ss1_ir_unique"), uq("ss2_ir_unique")
        fir_u, ipa_u = uq("junction_full_ir_unique"), uq("junction_ipa_unique")
        if strand == "-":
            s1u, s2u   = s2u,  s1u
            s1c, s2c   = s2c,  s1c
            s1ir, s2ir = s2ir, s1ir
            s1c_u, s2c_u = s2c_u, s1c_u
            s1ir_u, s2ir_u = s2ir_u, s1ir_u
        fir  = np.array([cov_result.get(j, {}).get("junction_full_ir", 0) for j in all_jxns], dtype=np.int32)
        ipa  = np.array([cov_result.get(j, {}).get("junction_ipa",     0) for j in all_jxns], dtype=np.int32)
        five_ss_arr  = [cov_result.get(j, {}).get("five_ss",  f"{chrom}_0") for j in all_jxns]
        three_ss_arr = [cov_result.get(j, {}).get("three_ss", f"{chrom}_0") for j in all_jxns]

        jcf  = jc.astype(np.float32)
        s1cf = s1c.astype(np.float32)
        s2cf = s2c.astype(np.float32)
        jc_approx_f = jc_approx.astype(np.float32)

        with np.errstate(divide="ignore", invalid="ignore"):
            psi_approx_v = np.where(jc_approx_f > 0, ju.astype(np.float32) / jc_approx_f, np.nan)
            psi_v    = np.where(jcf  > 0, ju.astype(np.float32) / jcf,  np.nan)
            irr1     = np.where(s1cf > 0, s1ir.astype(np.float32) / s1cf, np.nan)
            irr2     = np.where(s2cf > 0, s2ir.astype(np.float32) / s2cf, np.nan)
            full_irr = np.where(jcf  > 0, fir.astype(np.float32) / jcf,  np.nan)
            ipar     = np.where(s1cf > 0, ipa.astype(np.float32) / s1cf,  np.nan)

        row_dict = {
            "junction":                        all_jxns,
            "5ss":                             five_ss_arr,
            "3ss":                             three_ss_arr,
            "junction_usage":                  ju,
            "junction_usage_unique":           ju_u,
            "junction_coverage_approx":        jc_approx,
            "junction_coverage_approx_unique": jca_u,
            "junction_PSI_approx":             psi_approx_v,
            "rescaled_junction_PSI_approx":    _rescale(psi_approx_v),
            "junction_coverage":               jc,
            "junction_coverage_unique":        jc_u,
            "junction_PSI":                    psi_v,
            "rescaled_junction_PSI":           _rescale(psi_v),
            "5ss_usage":                       s1u,
            "5ss_coverage":                    s1c,
            "5ss_coverage_unique":             s1c_u,
            "5ss_IR_count":                    s1ir,
            "5ss_IR_count_unique":             s1ir_u,
            "5ss_IR_ratio":                    irr1,
            "rescaled_5ss_IR_ratio":           _rescale(irr1),
            "3ss_usage":                       s2u,
            "3ss_coverage":                    s2c,
            "3ss_coverage_unique":             s2c_u,
            "3ss_IR_count":                    s2ir,
            "3ss_IR_count_unique":             s2ir_u,
            "3ss_IR_ratio":                    irr2,
            "rescaled_3ss_IR_ratio":           _rescale(irr2),
            "junction_full_IR_count":          fir,
            "junction_full_IR_count_unique":   fir_u,
            "junction_full_IR_ratio":          full_irr,
            "rescaled_junction_full_IR_ratio": _rescale(full_irr),
            "junction_IPA_count":              ipa,
            "junction_IPA_count_unique":       ipa_u,
            "junction_IPA_ratio":              ipar,
            "rescaled_junction_IPA_ratio":     _rescale(ipar),
            "phasing":                         label,
        }
        rows.append(pd.DataFrame(row_dict))

    df = pd.concat(rows, ignore_index=True)
    df["sample"] = sample_name; df["region"] = region; df["gene"] = gene
    return df, time.time() - t0


# first pass: which junctions exist at all (bulk counts per sample)
def _discover_junctions_worker(
    sample_name: str, bulk_bam: str, region: str,
    include_monoexonic: bool = False,
) -> Tuple[str, Dict[Tuple[int, int], int], int]:
    if not os.path.exists(bulk_bam):
        print(f"[WARNING] Bulk BAM not found for sample '{sample_name}': {bulk_bam}. Skipping.")
        return sample_name, {}, 0
    jxn_raw, gene_cov, _ = collect_read_data(bulk_bam, region, include_monoexonic, keep_reads=False)
    return sample_name, jxn_raw, gene_cov


def _row_bam_size(row: pd.Series) -> int:
    total = 0
    for col in ("bulk", "hap1", "hap2"):
        val = row.get(col, None)
        if pd.notna(val) and isinstance(val, str):
            try:
                total += os.path.getsize(val)
            except OSError:
                pass
    return total


def _sample_args(row: pd.Series):
    bulk = row["bulk"] if pd.notna(row.get("bulk", "")) else None
    hap1 = row["hap1"] if pd.notna(row.get("hap1", "")) else None
    hap2 = row["hap2"] if pd.notna(row.get("hap2", "")) else None
    return row["sample"], bulk, hap1, hap2


# biggest BAMs first so one huge sample doesn't end up running alone at the end
def _rows_by_bam_size(region_df: pd.DataFrame) -> List[pd.Series]:
    rows = [pd.Series(r._asdict()) for r in region_df.itertuples(index=False)]
    return sorted(rows, key=_row_bam_size, reverse=True)


class _GeneJob:

    def __init__(self, gene, region, strand, region_df):
        self.gene, self.region, self.strand, self.region_df = gene, region, strand, region_df
        self.t0 = time.time()
        self.stage = None
        self.pending: set = set()
        self.failed = False
        self.sample_jxn_raw: Dict[str, Dict[Tuple[int, int], int]] = {}
        self.jxn_count_union: Dict[Tuple[int, int], int] = defaultdict(int)
        self.all_jxns: List[str] = []
        self.sample_dfs: List[pd.DataFrame] = []
        self.lines: List[str] = []


# several genes share one process pool, so workers don't sit idle between genes
def run_all_genes(gene_groups, approx_only, PSI_rescale_factor, threads, include_monoexonic,
                  min_jxn_reads, genome_path, alu, on_gene_done) -> None:
    max_inflight = max(2, threads)
    queue = list(gene_groups)
    owner: Dict[concurrent.futures.Future, Tuple[_GeneJob, str]] = {}
    active: List[_GeneJob] = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as pool:

        def submit_psi(job: _GeneJob):
            job.stage = "psi"
            for row in _rows_by_bam_size(job.region_df):
                sname, bulk, hap1, hap2 = _sample_args(row)
                if approx_only:
                    fut = pool.submit(process_sample, sname, bulk, hap1, hap2, job.region, job.gene,
                                      [], True, PSI_rescale_factor, job.strand, include_monoexonic)
                else:
                    fut = pool.submit(process_sample, sname, bulk, hap1, hap2, job.region, job.gene,
                                      job.all_jxns, approx_only, PSI_rescale_factor, job.strand,
                                      include_monoexonic, genome_path, alu,
                                      job.sample_jxn_raw.get(sname))
                owner[fut] = (job, sname); job.pending.add(fut)

        def start(job: _GeneJob):
            print(f"\n{'='*70}\n  Gene: {job.gene}   Region: {job.region}   (started)\n{'='*70}")
            if approx_only:
                submit_psi(job)
                if not job.pending:
                    stage_complete(job)
                return
            job.stage = "discovery"
            bulk_rows = job.region_df[job.region_df["bulk"].notna()]
            for _, row in bulk_rows.iterrows():
                fut = pool.submit(_discover_junctions_worker, row["sample"], row["bulk"],
                                  job.region, include_monoexonic)
                owner[fut] = (job, row["sample"]); job.pending.add(fut)
            if not job.pending:
                finish(job, None)

        def finish(job: _GeneJob, combined: Optional[pd.DataFrame]):
            for f in list(job.pending):
                f.cancel(); owner.pop(f, None)
            job.pending.clear()
            if job in active:
                active.remove(job)
            if combined is not None:
                combined = combined.sort_values(["junction", "sample", "phasing"], ignore_index=True)
            on_gene_done(job.gene, combined, time.time() - job.t0, job.lines)

        def stage_complete(job: _GeneJob):
            if job.stage == "discovery":
                chrom = job.region.split(":")[0]
                job.all_jxns = sorted(f"{chrom}_{s}_{e}" for (s, e), cnt in job.jxn_count_union.items()
                                      if cnt >= min_jxn_reads)
                job.lines.append(f"       {job.gene}: {len(job.jxn_count_union)} junctions found, "
                                 f"{len(job.all_jxns)} kept (with >={min_jxn_reads} reads in >=1 sample)")
                if not job.all_jxns:
                    job.lines.append("  No junctions passed filter. Skipping.")
                    finish(job, None); return
                submit_psi(job)
                if not job.pending:
                    finish(job, None)
                return
            if not job.sample_dfs:
                job.lines.append("  No data. Skipping."); finish(job, None); return
            combined = pd.concat(sorted(job.sample_dfs, key=lambda d: d["sample"].iloc[0]), ignore_index=True)
            if approx_only:
                all_jxns_set = set()
                bulk_max: Dict[str, int] = defaultdict(int)
                for sdf in job.sample_dfs:
                    all_jxns_set.update(sdf["junction"].unique())
                    b = sdf[sdf["phasing"] == "bulk"]
                    for jxn, u in zip(b["junction"], b["junction_usage"]):
                        u = int(u)
                        if u > bulk_max[jxn]: bulk_max[jxn] = u
                kept = {j for j in all_jxns_set if bulk_max.get(j, 0) >= min_jxn_reads}
                if len(all_jxns_set) - len(kept) > 0:
                    combined = combined[combined["junction"].isin(kept)]
                job.lines.append(f"       {job.gene}: {len(all_jxns_set)} junctions found, {len(kept)} kept "
                                 f"(>={min_jxn_reads}), {len(combined)} rows")
            finish(job, combined)

        while queue or active:
            while queue and len(active) < max_inflight:
                gene, region, strand, region_df = queue.pop(0)
                job = _GeneJob(gene, region, strand, region_df)
                active.append(job)
                start(job)
            live = [f for f in owner]
            if not live:
                continue
            done, _ = concurrent.futures.wait(live, return_when=concurrent.futures.FIRST_COMPLETED)
            for fut in done:
                if fut not in owner:
                    continue
                job, sname = owner.pop(fut)
                job.pending.discard(fut)
                if job.failed or job not in active:
                    continue
                try:
                    if job.stage == "discovery":
                        sname, jxn_raw, _ = fut.result()
                        job.sample_jxn_raw[sname] = jxn_raw
                        for coord, cnt in jxn_raw.items():
                            job.jxn_count_union[coord] = max(job.jxn_count_union[coord], cnt)
                    else:
                        sdf, elapsed = fut.result()
                        job.lines.append(f"       {sname}: {elapsed:.2f}s")
                        if len(sdf): job.sample_dfs.append(sdf)
                except Exception as e:
                    if job.stage == "discovery":
                        print(f"  [ERROR] Discovery ({job.gene}, {sname}): {e}"); traceback.print_exc()
                    else:
                        print(f"  [ERROR] {job.gene} / {sname}: {e}"); traceback.print_exc()
                        print(f"  [WARNING] Skipping {job.gene}.")
                        job.failed = True
                        finish(job, None)
                        continue
                if not job.pending:
                    stage_complete(job)


OUTPUT_COLUMNS = [
    "sample", "region", "gene", "phasing", "junction", "5ss", "3ss", "junction_type",
    "junction_usage", "junction_usage_unique",
    "junction_coverage_approx", "junction_coverage_approx_unique", "junction_PSI_approx", "rescaled_junction_PSI_approx",
    "junction_coverage", "junction_coverage_unique", "junction_PSI", "rescaled_junction_PSI",
    "5ss_usage", "5ss_IR_count", "5ss_IR_count_unique", "5ss_coverage", "5ss_coverage_unique",
    "5ss_IR_ratio", "rescaled_5ss_IR_ratio",
    "3ss_usage", "3ss_IR_count", "3ss_IR_count_unique", "3ss_coverage", "3ss_coverage_unique",
    "3ss_IR_ratio", "rescaled_3ss_IR_ratio",
    "junction_full_IR_count", "junction_full_IR_count_unique", "junction_full_IR_ratio",
    "rescaled_junction_full_IR_ratio",
    "junction_IPA_count", "junction_IPA_count_unique", "junction_IPA_ratio", "rescaled_junction_IPA_ratio",
]


def assign_junction_types(
    df: pd.DataFrame,
    gtf_junctions: Dict[str, Dict],
) -> pd.DataFrame:
    df = df.copy()
    jtype = pd.Series("novel", index=df.index, dtype=object)

    for gene, gdf in df.groupby("gene"):
        if gene not in gtf_junctions:
            continue
        canonical_set = gtf_junctions[gene]["canonical_junctions"]
        annotated_set = gtf_junctions[gene]["all_junctions"]
        jxns = df.loc[gdf.index, "junction"]
        is_annotated = jxns.isin(annotated_set) & ~jxns.isin(canonical_set)
        is_canonical = jxns.isin(canonical_set)
        jtype.loc[gdf.index[is_annotated]] = "annotated"
        jtype.loc[gdf.index[is_canonical]] = "canonical"

    df["junction_type"] = jtype
    return df


def finalize_columns(df: pd.DataFrame, gtf_junctions: Optional[Dict[str, Dict]]) -> pd.DataFrame:
    if gtf_junctions is not None:
        df = assign_junction_types(df, gtf_junctions)
    ordered = [c for c in OUTPUT_COLUMNS if c in df.columns]
    return df[ordered + [c for c in df.columns if c not in ordered]]


# gene -> result file, the downstream rules read this
def _write_manifest(gene_info: Dict[str, Tuple[str, str, str]],
                     result_paths: Dict[str, Optional[str]],
                     manifest_path: str) -> None:
    manifest_rows = [
        {"gene": gene, "result_path": result_paths.get(gene) or "None"}
        for gene in gene_info
    ]
    manifest_df = pd.DataFrame(manifest_rows, columns=["gene", "result_path"])
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    n_with_data = sum(1 for r in manifest_rows if r["result_path"] != "None")
    print(f"\nManifest written → {manifest_path} "
          f"({n_with_data}/{len(manifest_rows)} gene(s) have results)")


def main() -> None:
    print("\n" + "*"*80)
    print("  Cohort Junction Analysis (core metrics)")
    print("*"*80 + "\n")

    args        = parse_args()
    approx_only = args.approx

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.manifest)) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.note)) or ".", exist_ok=True)

    gene_info  = load_bed(args.bed)
    mapping_df = load_and_validate_mapping(args.mapping_file)

    n_samples = mapping_df["sample"].nunique()
    print(f"  {n_samples} sample(s) in this group.")

    missing = set(mapping_df["gene"].unique()) - set(gene_info)
    if missing:
        print(f"[WARNING] {len(missing)} gene(s) in mapping not in BED, skipping: {sorted(missing)}")

    gene_groups = []
    for gene, grp in mapping_df.groupby("gene", sort=False):
        if gene not in gene_info: continue
        _, region, strand = gene_info[gene]
        gene_groups.append((gene, region, strand, grp.reset_index(drop=True)))

    if args.test_n_genes is not None:
        gene_groups = gene_groups[:args.test_n_genes]
        print(f"[INFO] --test-n-genes {args.test_n_genes}: processing first {len(gene_groups)} gene(s) only.")

    alu = None
    if args.alu_bed:
        print(f"Loading Alu intervals from {args.alu_bed} ...")
        alu = load_alu_intervals(args.alu_bed)
        print(f"  → {sum(len(v) for v in alu.values())} intervals")

    genome_path = args.genome
    if genome_path:
        print(f"Genome: {genome_path} (IPA detection enabled)")
    else:
        print("[INFO] --genome not provided; IPA detection disabled.")

    n_genes = len(gene_groups)
    print(f"\nWill process {n_genes} gene(s)")
    print(f"Worker processes (shared across genes): {args.threads}\n")

    result_paths: Dict[str, Optional[str]] = {}

    gtf_junctions: Optional[Dict[str, Dict]] = None
    if args.gtf:
        print(f"Parsing GTF for annotated junctions: {args.gtf} ...")
        gtf_junctions = parse_gtf_junctions(args.gtf, [g[0] for g in gene_groups])
        print(f"  -> matched {sum(1 for g in gene_groups if g[0] in gtf_junctions)}/{len(gene_groups)} genes")
    else:
        print("[INFO] --gtf not provided; junction_type column will be omitted.")

    # write each gene out as soon as it finishes
    def on_gene_done(gene, combined, elapsed, lines):
        for line in sorted(lines): print(line)
        if combined is not None and len(combined):
            out_path = os.path.join(args.outdir, f"{gene}.tsv")
            finalize_columns(combined, gtf_junctions).to_csv(out_path, sep="\t", index=False)
            result_paths[gene] = out_path
            print(f"  {gene} → {out_path}")
        else:
            result_paths[gene] = None
            print(f"  {gene}: no output")
        print(f"  {gene} complete ({elapsed:.0f}s)")

    run_all_genes(gene_groups, approx_only, args.PSI_rescale_factor, args.threads,
                  args.include_monoexonic, args.min_jxn_reads, genome_path, alu, on_gene_done)

    _write_manifest(gene_info, result_paths, args.manifest)
    n_with_data = sum(1 for p in result_paths.values() if p)
    with open(args.note, "w") as fh:
        fh.write(f"Ran normally: {n_samples} sample(s), {len(gene_groups)} gene(s) "
                 f"attempted, {n_with_data} gene(s) produced results.\n")
    print("\nDone.")


if __name__ == "__main__":
    main()
