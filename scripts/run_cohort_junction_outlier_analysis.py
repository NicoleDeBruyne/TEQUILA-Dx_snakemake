#!/usr/bin/env python3

from __future__ import annotations

import os
import sys
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
from scipy.stats import betabinom, beta
from statsmodels.stats.multitest import multipletests
from pandas.errors import PerformanceWarning

warnings.filterwarnings("ignore", category=PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*Glyph.*missing from font.*", category=UserWarning)
warnings.filterwarnings("ignore", message=".*Adding colorbar to a different Figure.*", category=UserWarning)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_REF_OPS  = frozenset((0, 2, 3, 7, 8))   # M D N = X
_SKIP_OP  = 3                              # N  — intron

# Maps each metric column to its phasing-check denominator column
_METRIC_DENOMINATOR: Dict[str, str] = {
    "junction_PSI_approx":    "junction_coverage_approx",
    "junction_PSI":           "junction_coverage",
    "5ss_IR_ratio":           "5ss_coverage",
    "3ss_IR_ratio":           "3ss_coverage",
    "junction_full_IR_ratio": "junction_coverage",
    "junction_IPA_ratio":     "5ss_coverage",
}

# Maps each metric to the event types it can produce
_METRIC_EVENTS: Dict[str, List[str]] = {
    "junction_PSI_approx":    ["alt_5ss_approx", "alt_3ss_approx", "exon_skipping_approx", "exon_inclusion_approx"],
    "junction_PSI":           ["alt_5ss", "alt_3ss", "exon_skipping", "exon_inclusion"],
    "5ss_IR_ratio":           ["5ss_IR"],
    "3ss_IR_ratio":           ["3ss_IR"],
    "junction_full_IR_ratio": ["full_IR"],
    "junction_IPA_ratio":     ["IPA"],
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Identifies splice junctions with unusual usage frequencies."
    )
    p.add_argument("--mapping-file",               required=True)
    p.add_argument("--bed",                        required=True)
    p.add_argument("--outprefix",                  required=True)
    p.add_argument("--approx",                     action="store_true")
    p.add_argument("--coverage-threshold",         type=int,   default=20)
    p.add_argument("--PSI-rescale-factor",         type=float, default=1e-3)
    p.add_argument("--n-threshold",                type=int,   default=30)
    p.add_argument("--phasing-threshold",          type=float, default=0.8)
    p.add_argument("--thresholds",
                   nargs="*", default=[],
                   metavar="PADJ:DELTA",
                   help="One or more padj:delta threshold pairs, e.g. 0.05:0.1 0.01:0.1 0.01:0.2. "
                        "Each combination produces its own output subdirectory.")
    p.add_argument("--no-ss-IR",                   action="store_true")
    p.add_argument("--min-jxn-reads",              type=int,   default=20)
    p.add_argument("--include-monoexonic",         action="store_true")
    p.add_argument("--genome",                     default=None)
    p.add_argument("--alu-bed",                    default=None)
    p.add_argument("--gtf",                        default=None,
                   help="GTF/GTF.gz. If provided, adds junction_type column.")
    p.add_argument("--threads",                    type=int,   default=1)
    p.add_argument("--test-n-genes",               type=int,   default=None)
    p.add_argument("--loo",                        action="store_true",
                   help="Leave-one-out Beta fitting: refit the distribution "
                        "excluding each bulk sample before testing it. Removes "
                        "self-contamination bias; costs ~n_samples× more fits.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# BED / mapping I/O
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# BAM helpers
# ---------------------------------------------------------------------------

def _parse_region(region: str) -> Tuple[str, int, int]:
    chrom, se = region.split(":")
    start, end = map(int, se.split("-"))
    return chrom, start, end


def _build_ir_anchor_structures(
    jxn_coords: List[Tuple[int, int]],
) -> Tuple[List[int], List[int], List[List[int]], List[List[int]], List[int], List[int]]:
    """Build sorted lookup structures for full-IR and IPA anchor checks (unchanged)."""
    ss1_map: Dict[int, List[int]] = defaultdict(list)
    ss2_map: Dict[int, List[int]] = defaultdict(list)
    ss1_hi:  Dict[int, int]       = {}
    ss2_hi:  Dict[int, int]       = {}
    for i, (ss1, ss2) in enumerate(jxn_coords):
        lo1 = ss1 - 3;  ss1_map[lo1].append(i);  ss1_hi[lo1] = ss1 - 1
        lo2 = ss2;        ss2_map[lo2].append(i);  ss2_hi[lo2] = ss2 + 3
    ss1_sorted = sorted(ss1_map);  ss2_sorted = sorted(ss2_map)
    return (
        ss1_sorted, ss2_sorted,
        [ss1_map[p] for p in ss1_sorted],
        [ss2_map[p] for p in ss2_sorted],
        [ss1_hi[p]  for p in ss1_sorted],
        [ss2_hi[p]  for p in ss2_sorted],
    )


# CIGAR operators that count as "covered" at a position
# (M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8)
_COV_OPS_EXON_SET   = frozenset((0, 1, 2, 7, 8))   # M, I, D, =, X
_COV_OPS_INTRON_SET = frozenset((0, 1, 2, 3, 7, 8)) # M, I, D, N, =, X

# Sentinel value for positions not covered by the read
_NO_OP = 255

# Pre-built uint8 lookup arrays (index = CIGAR op code)
_IS_EXON_OP   = np.zeros(9, dtype=np.uint8)
_IS_INTRON_OP = np.zeros(9, dtype=np.uint8)
for _op in _COV_OPS_EXON_SET:
    if _op < 9: _IS_EXON_OP[_op]   = 1
for _op in _COV_OPS_INTRON_SET:
    if _op < 9: _IS_INTRON_OP[_op] = 1



def collect_read_data(
    bam_path: str,
    region: str,
    include_monoexonic: bool = False,
    collect_softclips: bool = False,
    strand: str = "+",
    genome_seq: Optional[str] = None,
    gene_region_start: int = 0,
) -> Tuple[Dict[Tuple[int, int], int], int,
           List[Tuple[List[Tuple[int, int]], List[int], Optional[Tuple]]],
           ]:
    """
    Step 1 BAM walk: collect junction counts, read blocks, splice positions,
    and soft-clip info for every qualifying read.

    Returns:
        jxn_raw      – {(ss1, ss2): count}
        gene_cov     – qualifying read count
        reads        – list of (blocks, splice_junctions, softclip_tuple_or_None)
                       blocks: [(start, end), ...] 0-based half-open per exon
                       splice_junctions: [(ss1, ss2), ...] 1-based per intron
                       softclip_tuple: (sc, g, pos3, leading) or None
    """
    chrom_r, region_start, region_end = _parse_region(region)
    jxn_raw:  Dict[Tuple[int, int], int] = defaultdict(int)
    gene_cov  = 0
    reads:    List = []

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
                elif op == 4:   # S
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
                reads.append((blocks, [], None))
                continue

            # Soft-clip info
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

            reads.append((blocks, read_jxns, sc_tuple))

    return jxn_raw, gene_cov, reads


def compute_coverage_metrics(
    reads:              List,
    jxn_coords_for_cov: List[Tuple[int, int]],
    strand:             str,
    genome_seq:         Optional[str],
    gene_region_start:  int,
    alu:                Optional[dict],
    chrom:              str,
) -> Dict[Tuple[int, int], List[int]]:
    """
    Step 2: vectorized coverage computation over pre-collected read data.

    For each read (blocks, splice_junctions, sc_tuple):
      ss1_cov: 3-base exon window [ss1-4,ss1-2] (0-based) fully within a block
               AND read overlaps intron [ss1-1, ss1+1] (0-based)
               i.e. ref_start <= ss1-2 and ref_end > ss1+1
      ss2_cov: 3-base exon window [ss2,ss2+2] (0-based) fully within a block
               AND read overlaps intron [ss2-3, ss2-1] (0-based)
               i.e. ref_start <= ss2-3 and ref_end > ss2+2
      junction_cov: ss1_cov OR ss2_cov
      ss1_ir:  all 6 positions [ss1-4, ss1+1] (0-based) within a single block
      ss2_ir:  all 6 positions [ss2-3, ss2+2] (0-based) within a single block
      full_ir: single block spans [ss1-4, ss2+2] (0-based, i.e. block_start<=ss1-4
               and block_end>=ss2+3)
      ipa:     passes 5ss_cov (ss1_cov on +, ss2_cov on -)
               AND ref_end (0-based excl) satisfies ss1 <= ref_end <= ss2
               (1-based: last covered base is within intron)
               AND no splice junction at or past 5ss
               AND poly-A soft-clip

    Returns: {(ss1,ss2): [ss1_cov, ss2_cov, jxn_cov, ss1_ir, ss2_ir, full_ir, ipa]}
    """
    n_jxns = len(jxn_coords_for_cov)
    if n_jxns == 0:
        return {}

    # Per-junction accumulator arrays (indexed by junction position in jxn_coords_for_cov)
    ss1_cov  = np.zeros(n_jxns, dtype=np.int32)
    ss2_cov  = np.zeros(n_jxns, dtype=np.int32)
    jxn_cov  = np.zeros(n_jxns, dtype=np.int32)
    ss1_ir   = np.zeros(n_jxns, dtype=np.int32)
    ss2_ir   = np.zeros(n_jxns, dtype=np.int32)
    full_ir  = np.zeros(n_jxns, dtype=np.int32)
    ipa_arr  = np.zeros(n_jxns, dtype=np.int32)

    # Pre-compute per-junction coordinate arrays (0-based)
    ss1_arr  = np.array([ss1     for ss1, ss2 in jxn_coords_for_cov], dtype=np.int64)
    ss2_arr  = np.array([ss2     for ss1, ss2 in jxn_coords_for_cov], dtype=np.int64)
    # ss1 coverage: exon window [ss1-4, ss1-2], intron overlap: ref_start<=ss1-2, ref_end>ss1+1
    ss1_exon_lo = ss1_arr - 4   # 0-based start of exon window
    ss1_exon_hi = ss1_arr - 2   # 0-based end of exon window (inclusive)
    ss1_intr_hi = ss1_arr + 1   # 0-based end of intron overlap (inclusive)
    # ss2 coverage: exon window [ss2, ss2+2], intron overlap: ref_start<=ss2-3, ref_end>ss2+2
    ss2_exon_lo = ss2_arr       # 0-based start of exon window
    ss2_exon_hi = ss2_arr + 2   # 0-based end of exon window (inclusive)
    ss2_intr_lo = ss2_arr - 3   # 0-based: ref_start must be <= ss2-3 (covers intron positions ss2-2..ss2 in 1-based)
    # ss1_ir: single block spans [ss1-4, ss1+1] (0-based inclusive)
    ss1_ir_lo = ss1_arr - 4
    ss1_ir_hi = ss1_arr + 1    # 0-based inclusive → block_end > ss1+1
    # ss2_ir: single block spans [ss2-3, ss2+2] (0-based inclusive)
    ss2_ir_lo = ss2_arr - 3
    ss2_ir_hi = ss2_arr + 2    # 0-based inclusive → block_end > ss2+2
    # full_ir: single block spans [ss1-4, ss2+2] (0-based inclusive)
    fir_lo = ss1_arr - 4
    fir_hi = ss2_arr + 2       # 0-based inclusive → block_end > ss2+2

    # 5ss is ss1 for + strand, ss2 for - strand (1-based)
    five_ss_arr = ss1_arr if strand == "+" else ss2_arr

    for blocks, splice_jxns, sc_tuple in reads:
        if not blocks:
            continue
        ref_start = blocks[0][0]   # 0-based
        ref_end   = blocks[-1][1]  # 0-based exclusive

        # Check ss1_cov for all junctions:
        # exon window [ss1-4, ss1-2] must be fully within a block
        # AND read overlaps intron (ref_end > ss1+1)
        # First gate: read must reach exon window and intron start
        cand_ss1 = ((ref_start <= ss1_exon_lo) &   # wait — exon window at ss1-4..ss1-2
                    (ref_end   >  ss1_intr_hi))     # ref_end > ss1+1 (0-based excl)
        # Actually ref_start just needs to be <= ss1_exon_lo is wrong —
        # we need the block to CONTAIN [ss1_exon_lo, ss1_exon_hi]
        # i.e. block_start <= ss1_exon_lo and block_end > ss1_exon_hi
        # Check across all blocks
        h1 = np.zeros(n_jxns, dtype=bool)
        h2 = np.zeros(n_jxns, dtype=bool)
        h1_ir_mask = np.zeros(n_jxns, dtype=bool)
        h2_ir_mask = np.zeros(n_jxns, dtype=bool)
        fir_mask   = np.zeros(n_jxns, dtype=bool)

        for block_start, block_end in blocks:
            bs = np.int64(block_start)
            be = np.int64(block_end)   # 0-based exclusive
            # ss1_cov: block contains exon window AND read overlaps intron
            mask1 = (bs <= ss1_exon_lo) & (be > ss1_exon_hi) & (ref_end > ss1_intr_hi)
            h1   |= mask1
            # ss2_cov: block contains exon window AND read overlaps intron
            mask2 = (bs <= ss2_exon_lo) & (be > ss2_exon_hi) & (ref_start <= ss2_intr_lo)
            h2   |= mask2
            # ss1_ir: single block spans [ss1-4, ss1+1]
            h1_ir_mask |= (bs <= ss1_ir_lo) & (be > ss1_ir_hi)
            # ss2_ir: single block spans [ss2-3, ss2+2]
            h2_ir_mask |= (bs <= ss2_ir_lo) & (be > ss2_ir_hi)
            # full_ir: single block spans [ss1-4, ss2+2]
            fir_mask   |= (bs <= fir_lo)    & (be > fir_hi)

        ss1_cov  += h1.astype(np.int32)
        ss2_cov  += h2.astype(np.int32)
        jxn_cov  += (h1 | h2).astype(np.int32)
        ss1_ir   += h1_ir_mask.astype(np.int32)
        ss2_ir   += h2_ir_mask.astype(np.int32)
        full_ir  += fir_mask.astype(np.int32)

        # IPA: must pass 5ss_cov, terminate within intron, no downstream splice
        if sc_tuple is None:
            continue
        # 5ss_cov: ss1_cov for + strand, ss2_cov for - strand
        five_ss_cov = h1 if strand == "+" else h2
        if not five_ss_cov.any():
            continue
        # ref_end (0-based excl) must satisfy ss1 <= ref_end <= ss2
        # i.e. last covered base (1-based) is within intron
        # 0-based: ss1-1 < ref_end <= ss2  →  ss1 <= ref_end <= ss2 (since ref_end is excl)
        within_intron = (ss1_arr <= ref_end) & (ref_end <= ss2_arr)
        cand_ipa = five_ss_cov & within_intron
        if not cand_ipa.any():
            continue
        # No splice junction at or past 5ss (1-based)
        splice_set = set(j[0] if strand == "+" else j[1] for j in splice_jxns)
        sc, g, pos3, leading = sc_tuple
        is_polya = _is_oligo_dt_priming(sc, g, leading)
        if not is_polya and alu is not None:
            frag = sc[-10:] if leading else sc[:10]
            if len(frag) == 10 and (frag.count("A")/10 >= 0.8 or frag.count("T")/10 >= 0.8):
                if _in_alu(chrom, pos3, alu):
                    is_polya = True
        if not is_polya:
            continue
        for ji in np.where(cand_ipa)[0]:
            five_ss = int(five_ss_arr[ji])
            if not any(s >= five_ss for s in splice_set):
                ipa_arr[ji] += 1

    # Build result dict
    result = {}
    for i, jc in enumerate(jxn_coords_for_cov):
        result[jc] = [
            int(ss1_cov[i]), int(ss2_cov[i]), int(jxn_cov[i]),
            int(ss1_ir[i]),  int(ss2_ir[i]),  int(full_ir[i]),
            int(ipa_arr[i]),
        ]
    return result


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


def compute_junction_coverage_approx(
    junction_counts: Dict[str, int],
    all_jxns: Optional[List[str]] = None,
) -> Dict[str, int]:
    """
    Compute approximate junction coverage (ss1_usage + ss2_usage - junction_usage)
    for every junction.

    all_jxns: cohort-wide junction list.  When provided, coverage is computed
    for every junction in all_jxns even if it has 0 reads in this sample —
    using the ss1/ss2 usage tallied from whatever junctions ARE present.
    Without all_jxns, only junctions in junction_counts are returned (old behaviour).
    """
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

    # Build result for cohort-wide junction list if provided
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


# ---------------------------------------------------------------------------
# Per-sample worker
# ---------------------------------------------------------------------------

def process_sample(
    sample_name: str,
    bulk_bam:    Optional[str],
    hap1_bam:    Optional[str],
    hap2_bam:    Optional[str],
    region: str,
    gene: str,
    all_jxns: List[str],
    approx_only: bool,
    coverage_threshold: int,
    phasing_threshold: float,
    PSI_rescale_factor: float,
    strand: str = "+",
    include_monoexonic: bool = False,
    genome_path: Optional[str] = None,
    alu: Optional[dict] = None,
    no_ss_ir: bool = False,
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
        raw: List[Tuple[str, int, dict]] = []
        for label, bam in (("bulk", bulk_bam), ("hap1", hap1_bam), ("hap2", hap2_bam)):
            if not bam or not isinstance(bam, str): continue
            jxn_counts, gene_cov, _ = collect_read_data(
                bam, region, include_monoexonic)
            if gene_cov == 0: continue
            raw.append((label, gene_cov, jxn_counts))
        if not raw: return empty, time.time() - t0

        if not all_jxns:
            uni: set = set()
            for _, _, jc in raw: uni.update(jc)
            all_jxns = sorted(uni)

        jxn_index = {j: i for i, j in enumerate(all_jxns)}
        n_jxns = len(all_jxns)
        chunks = []
        for label, gene_cov, jxn_counts in raw:
            cov_dict  = compute_junction_coverage_approx(jxn_counts)
            u = np.zeros(n_jxns, dtype=np.int32)
            c = np.zeros(n_jxns, dtype=np.int32)
            for jxn, cnt in jxn_counts.items():
                if jxn in jxn_index: u[jxn_index[jxn]] = cnt
            for jxn, cnt in cov_dict.items():
                if jxn in jxn_index: c[jxn_index[jxn]] = cnt
            chunks.append((label, u, c))

        # Per-metric phasing filter for PSI_approx — denominator = junction_coverage_approx
        label_cov = {label: c for label, _, c in chunks}
        lp_psi_approx = np.zeros(n_jxns, dtype=bool)
        if phasing_threshold > 0 and ("hap1" in label_cov or "hap2" in label_cov):
            bulk_cov  = label_cov.get("bulk", np.zeros(n_jxns, dtype=np.int32))
            hap1_cov  = label_cov.get("hap1", np.zeros(n_jxns, np.int32))
            hap2_cov  = label_cov.get("hap2", np.zeros(n_jxns, np.int32))
            hap_total = hap1_cov + hap2_cov
            lp_psi_approx = (
                (bulk_cov * phasing_threshold > hap_total) |
                (hap1_cov < coverage_threshold) |
                (hap2_cov < coverage_threshold)
            )

        def _lc_approx(arr):
            out = arr.astype(object); out[~np.isfinite(arr)] = "low_coverage"; return out

        rows = []
        for label, u, c in chunks:
            is_hap = label in ("hap1", "hap2")
            lp_col = lp_psi_approx if is_hap else np.zeros(n_jxns, dtype=bool)
            cf = c.astype(np.float32); uf = u.astype(np.float32)
            with np.errstate(divide="ignore", invalid="ignore"):
                psi = np.where(cf >= coverage_threshold, uf / cf, np.nan)
            rsc = np.where(np.isfinite(psi), psi*(1-2*PSI_rescale_factor)+PSI_rescale_factor, np.nan)
            rows.append(pd.DataFrame({
                "junction":          all_jxns,
                "junction_usage":    u,
                "junction_coverage_approx": c,
                "junction_PSI_approx":          _lc_approx(psi),
                "rescaled_junction_PSI_approx": _lc_approx(rsc),
                "phasing":           label,
                "low_phased_junction_PSI_approx": lp_col,
            }))
        df = pd.concat(rows, ignore_index=True)
        df["sample"] = sample_name; df["region"] = region; df["gene"] = gene
        return df, time.time() - t0

    # ---- Full metrics path ----
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
    do_ss_ir = not no_ss_ir

    # Build jxn_coords_for_cov from all_jxns
    jxn_coords_for_cov: List[Tuple[int, int]] = []
    for jxn in all_jxns:
        parts = jxn.split("_")
        jxn_coords_for_cov.append((int(parts[-2]), int(parts[-1])))

    # Step 1: collect jxn_raw and read data for all BAMs
    # For bulk, use pre-collected jxn_raw_bulk if available (from discovery step)
    raw_step1: List = []
    for label, bam in (("bulk", bulk_bam), ("hap1", hap1_bam), ("hap2", hap2_bam)):
        if not bam or not isinstance(bam, str): continue
        if label == "bulk" and jxn_raw_bulk is not None:
            # Re-use jxn_raw from Step 1; still need to walk for read data
            jxn_raw_lbl = jxn_raw_bulk
            _, gene_cov, reads = collect_read_data(
                bam, region, include_monoexonic,
                collect_softclips=do_ipa,
                strand=strand,
                genome_seq=genome_seq,
                gene_region_start=gene_start,
            )
        else:
            jxn_raw_lbl, gene_cov, reads = collect_read_data(
                bam, region, include_monoexonic,
                collect_softclips=do_ipa,
                strand=strand,
                genome_seq=genome_seq,
                gene_region_start=gene_start,
            )
        if gene_cov == 0: continue
        raw_step1.append((label, gene_cov, jxn_raw_lbl, reads))
    if not raw_step1: return empty, time.time() - t0

    # Step 2: compute coverage metrics vectorized over read data
    chunks_wide: List = []
    for label, gene_cov, jxn_raw_lbl, reads in raw_step1:
        # Fingerprints (read diversity)
        jxn_fingerprints: Dict[Tuple[int, int], Dict] = defaultdict(lambda: defaultdict(int))
        for blocks, splice_jxns, sc_tuple in reads:
            if splice_jxns:
                fp = (blocks[0][0], blocks[-1][1], tuple(splice_jxns))
                for jc in splice_jxns:
                    jxn_fingerprints[jc][fp] += 1

        # Shannon N_eff diversity
        import math
        jxn_diversity: Dict[Tuple[int, int], float] = {}
        for jc, fp_counts in jxn_fingerprints.items():
            total = sum(fp_counts.values())
            if total == 0:      jxn_diversity[jc] = 0.0
            elif len(fp_counts) == 1: jxn_diversity[jc] = 1.0
            else:
                h = -sum((c/total)*math.log(c/total) for c in fp_counts.values())
                jxn_diversity[jc] = round(math.exp(h), 2)

        cov_metrics = compute_coverage_metrics(
            reads, jxn_coords_for_cov, strand,
            genome_seq, gene_start, alu, chrom,
        )

        # Build cov_result dict keyed by junction string
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
            m = cov_metrics.get((ss1, ss2), [0]*7)
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
                "junction_read_diversity": (jxn_diversity.get((ss1, ss2), float("nan"))),
            }

        jxn_raw_str = {f"{chrom}_{s}_{e}": c for (s, e), c in jxn_raw_lbl.items()}
        approx_cov  = compute_junction_coverage_approx(jxn_raw_str, all_jxns)
        chunks_wide.append((label, cov_result, jxn_raw_str, approx_cov))

    # ----- Per-metric phasing filter -----
    # For each metric we need its denominator arrays per label.
    # We build arrays for: junction_coverage, junction_coverage_approx,
    # 5ss_coverage, 3ss_coverage.
    # Then for each metric, combine hap1 and hap2 denominators against bulk.

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

    # Collect coverage arrays per denominator column
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

    def _compute_lp(denom_col: str) -> np.ndarray:
        """Return boolean array: True where junction fails phasing for this denominator."""
        arrs = denom_arrays.get(denom_col, {})
        if not arrs or phasing_threshold <= 0:
            return np.zeros(n_jxns, dtype=bool)
        b = arrs.get("bulk", np.zeros(n_jxns, np.int32))
        h1 = arrs.get("hap1", np.zeros(n_jxns, np.int32))
        h2 = arrs.get("hap2", np.zeros(n_jxns, np.int32))
        return (
            (b * phasing_threshold > (h1 + h2)) |
            (h1 < coverage_threshold) |
            (h2 < coverage_threshold)
        )

    # One low_phased array per metric
    lp_per_metric: Dict[str, np.ndarray] = {}
    for mc, denom_col in _METRIC_DENOMINATOR.items():
        lp_per_metric[mc] = _compute_lp(denom_col)

    def _rescale(arr_f):
        return np.where(np.isfinite(arr_f),
                        arr_f * (1 - 2*PSI_rescale_factor) + PSI_rescale_factor, np.nan)

    def _fill_lc(arr):
        out = arr.astype(object)
        out[~np.isfinite(arr)] = "low_coverage"
        return out

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
        if strand == "-":
            s1u, s2u   = s2u,  s1u
            s1c, s2c   = s2c,  s1c
            s1ir, s2ir = s2ir, s1ir
        fir  = np.array([cov_result.get(j, {}).get("junction_full_ir", 0) for j in all_jxns], dtype=np.int32)
        ipa  = np.array([cov_result.get(j, {}).get("junction_ipa",     0) for j in all_jxns], dtype=np.int32)
        five_ss_arr  = [cov_result.get(j, {}).get("five_ss",  f"{chrom}_0") for j in all_jxns]
        three_ss_arr = [cov_result.get(j, {}).get("three_ss", f"{chrom}_0") for j in all_jxns]

        is_hap = label in ("hap1", "hap2")

        jcf  = jc.astype(np.float32)
        s1cf = s1c.astype(np.float32)
        s2cf = s2c.astype(np.float32)
        jc_approx_f = jc_approx.astype(np.float32)

        with np.errstate(divide="ignore", invalid="ignore"):
            psi_approx_v = np.where(jc_approx_f >= coverage_threshold,
                                    ju.astype(np.float32) / jc_approx_f, np.nan)
            psi_v    = np.where(jcf  >= coverage_threshold, ju.astype(np.float32) / jcf,  np.nan)
            irr1     = np.where(s1cf >= coverage_threshold, s1ir.astype(np.float32) / s1cf, np.nan)
            irr2     = np.where(s2cf >= coverage_threshold, s2ir.astype(np.float32) / s2cf, np.nan)
            full_irr = np.where(jcf  >= coverage_threshold, fir.astype(np.float32) / jcf,  np.nan)
            ipar     = np.where(s1cf >= coverage_threshold, ipa.astype(np.float32) / s1cf,  np.nan)

        row_dict = {
            "junction":                        all_jxns,
            "5ss":                             five_ss_arr,
            "3ss":                             three_ss_arr,
            "junction_usage":                  ju,
            "junction_read_diversity":         [
                0.0 if ju[i] == 0
                else cov_result.get(j, {}).get("junction_read_diversity", float("nan"))
                for i, j in enumerate(all_jxns)
            ],
            "junction_coverage_approx":        jc_approx,
            "junction_PSI_approx":             _fill_lc(psi_approx_v),
            "rescaled_junction_PSI_approx":    _fill_lc(_rescale(psi_approx_v)),
            "junction_coverage":               jc,
            "junction_PSI":                    _fill_lc(psi_v),
            "rescaled_junction_PSI":           _fill_lc(_rescale(psi_v)),
            "5ss_usage":         s1u,
            "5ss_coverage":      s1c,
            "5ss_IR_ratio":                    _fill_lc(irr1),
            "rescaled_5ss_IR_ratio":           _fill_lc(_rescale(irr1)),
            "3ss_usage":         s2u,
            "3ss_coverage":      s2c,
            "3ss_IR_ratio":                    _fill_lc(irr2),
            "rescaled_3ss_IR_ratio":           _fill_lc(_rescale(irr2)),
            "junction_full_IR_count":          fir,   # singular
            "junction_full_IR_ratio":          _fill_lc(full_irr),
            "rescaled_junction_full_IR_ratio": _fill_lc(_rescale(full_irr)),
            "junction_IPA_count":              ipa,   # singular
            "junction_IPA_ratio":              _fill_lc(ipar),
            "rescaled_junction_IPA_ratio":     _fill_lc(_rescale(ipar)),
            "phasing":                         label,
        }
        # Attach per-metric low_phased columns
        for mc in _METRIC_DENOMINATOR:
            col_name = f"low_phased_{mc}"
            lp_arr = lp_per_metric[mc] if is_hap else np.zeros(n_jxns, dtype=bool)
            row_dict[col_name] = lp_arr

        rows.append(pd.DataFrame(row_dict))

    df = pd.concat(rows, ignore_index=True)
    df["sample"] = sample_name; df["region"] = region; df["gene"] = gene
    return df, time.time() - t0


# ---------------------------------------------------------------------------
# Junction discovery worker
# ---------------------------------------------------------------------------

def _discover_junctions_worker(
    sample_name: str, bulk_bam: str, region: str,
    include_monoexonic: bool = False,
) -> Tuple[str, Dict[Tuple[int, int], int], int]:
    if not os.path.exists(bulk_bam):
        print(f"[WARNING] Bulk BAM not found for sample '{sample_name}': {bulk_bam}. Skipping.")
        return sample_name, {}, 0
    jxn_raw, gene_cov, _ = collect_read_data(bulk_bam, region, include_monoexonic)
    return sample_name, jxn_raw, gene_cov


# ---------------------------------------------------------------------------
# Beta distribution fitting
# ---------------------------------------------------------------------------

def _fit_one_beta(x: np.ndarray, tol: float, n_threshold: int):
    x = x[np.isfinite(x)]
    n = len(x)
    if n < n_threshold: return n, "low_n", "low_n", "low_n"
    var = float(np.var(x))
    if var < tol:
        m = float(np.clip(np.median(x), tol, 1.0 - tol))
        k = max(m * (1.0 - m) / tol - 1.0, tol)
        a = m * k; b = (1.0 - m) * k
        return n, float(a), float(b), float(m)
    try:
        a, b, *_ = beta.fit(x, floc=0, fscale=1)
        return n, float(a), float(b), float(a / (a + b))
    except Exception:
        return n, "error", "error", "error"


def _fit_beta_rows(args):
    """Worker: fit a sub-block of rows. Returns list of result tuples."""
    mat_block, tol, n_threshold = args
    return [_fit_one_beta(mat_block[i], tol, n_threshold)
            for i in range(len(mat_block))]


def fit_beta_dist_chunk(mat, feat_names, tol, n_threshold, threads: int = 1):
    """
    Fit beta distributions for every row in mat.

    Parallelised across rows using a process pool when threads > 1.
    Results are collected via executor.map, which preserves submission order,
    so output is deterministic regardless of the number of workers.
    """
    n = len(mat)
    if n == 0:
        return pd.DataFrame(columns=["n", "alpha", "beta_param", "expected"])

    n_workers = min(threads, n)
    if n_workers <= 1:
        results = [_fit_one_beta(mat[i], tol, n_threshold) for i in range(n)]
    else:
        chunk_size = max(1, (n + n_workers - 1) // n_workers)
        chunks = [
            (mat[i : i + chunk_size], tol, n_threshold)
            for i in range(0, n, chunk_size)
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
            # map() preserves submission order — results are deterministic
            results = []
            for block in ex.map(_fit_beta_rows, chunks):
                results.extend(block)

    return pd.DataFrame(results, index=feat_names,
                        columns=["n", "alpha", "beta_param", "expected"])


# ---------------------------------------------------------------------------
# Beta-binomial test
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Beta-binomial test
# ---------------------------------------------------------------------------

def _betabinom_test_rows(args):
    """
    Worker: run betabinom test on a sub-DataFrame.
    Receives a dict of arrays (picklable) rather than a DataFrame.
    Returns a 1-D array of p-values (NaN where not testable).
    """
    usage, coverage, alpha_v, beta_v, val, cov_thresh = args
    valid = (np.isfinite(usage) & np.isfinite(coverage) &
             np.isfinite(alpha_v) & np.isfinite(beta_v) & np.isfinite(val) &
             (coverage >= cov_thresh))
    p = np.full(len(usage), np.nan)
    if valid.any():
        u   = np.round(usage[valid]).astype(int)
        cv  = np.round(coverage[valid]).astype(int)
        lte = betabinom.cdf(u, cv, alpha_v[valid], beta_v[valid])
        gte = betabinom.cdf(cv - u, cv, beta_v[valid], alpha_v[valid])
        p[valid] = np.clip(2.0 * np.minimum(lte, gte), 0.0, 1.0)
    return p


def beta_binomial_test_chunk(
    df: pd.DataFrame,
    coverage_threshold: int,
    metric_col: str,
    usage_col: str,
    coverage_col: str,
    p_col: str,
    threads: int = 1,
) -> pd.DataFrame:
    """
    Run betabinom tests for every row in df.

    Parallelised across rows using a process pool when threads > 1.
    Results are collected via executor.map (order-preserving) so p-values
    are identical regardless of worker count.
    """
    df = df.copy()
    usage    = df[usage_col].to_numpy(dtype=float)
    coverage = df[coverage_col].to_numpy(dtype=float)
    alpha_v  = pd.to_numeric(df[f"alpha_{metric_col}"], errors="coerce").to_numpy()
    beta_v   = pd.to_numeric(df[f"beta_{metric_col}"],  errors="coerce").to_numpy()
    val      = df[metric_col].to_numpy(dtype=float)

    n = len(df)
    n_workers = min(threads, n)

    if n_workers <= 1:
        p = _betabinom_test_rows(
            (usage, coverage, alpha_v, beta_v, val, coverage_threshold)
        )
    else:
        chunk_size = max(1, (n + n_workers - 1) // n_workers)
        chunks = [
            (
                usage   [i : i + chunk_size],
                coverage[i : i + chunk_size],
                alpha_v [i : i + chunk_size],
                beta_v  [i : i + chunk_size],
                val     [i : i + chunk_size],
                coverage_threshold,
            )
            for i in range(0, n, chunk_size)
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
            # map() preserves submission order — p-values are deterministic
            p = np.concatenate(list(ex.map(_betabinom_test_rows, chunks)))

    df[p_col] = p
    return df


# ---------------------------------------------------------------------------
# Per-metric beta pipeline
# ---------------------------------------------------------------------------

def _run_one_metric(
    combined_df:        pd.DataFrame,
    metric_col:         str,
    rescaled_col:       str,
    usage_col:          str,
    coverage_col:       str,
    id_col:             str,
    coverage_threshold: int,
    PSI_rescale_factor: float,
    n_threshold:        int,
    threads:            int,
    loo:                bool = False,
) -> pd.DataFrame:
    is_ss_metric = metric_col in ("5ss_IR_ratio", "3ss_IR_ratio")

    if is_ss_metric:
        ss_pos_col = "5ss" if metric_col == "5ss_IR_ratio" else "3ss"
        combined_df = combined_df.copy()
        combined_df["_ss_id"] = combined_df["region"].str.split(":").str[0] + "_" + \
                                  combined_df[ss_pos_col].astype(str)
        fit_id_col = "_ss_id"
    else:
        fit_id_col = id_col

    bulk_sub = combined_df[combined_df["phasing"] == "bulk"][
        ["sample", fit_id_col, coverage_col, rescaled_col]
    ].copy()

    low_cov = bulk_sub[coverage_col].to_numpy(dtype=float) < coverage_threshold
    rvals   = pd.to_numeric(bulk_sub[rescaled_col], errors="coerce").to_numpy(dtype=float)
    rvals[low_cov] = np.nan

    bulk_wide = (
        bulk_sub.assign(**{rescaled_col: rvals})
        .drop_duplicates(subset=["sample", fit_id_col])
        .pivot(index=fit_id_col, columns="sample", values=rescaled_col)
        .astype(np.float32)
    )

    n_col        = f"n_{metric_col}"
    alpha_col    = f"alpha_{metric_col}"
    beta_col     = f"beta_{metric_col}"
    expected_col = f"expected_{metric_col}"
    delta_col    = f"delta_{metric_col}"
    p_col        = f"p_value_{metric_col}"
    p1_col       = f"p1_{metric_col}"
    p99_col      = f"p99_{metric_col}"

    if bulk_wide.empty:
        combined_df[n_col] = 0
        for col in (alpha_col, beta_col, expected_col, p1_col, p99_col, delta_col, p_col):
            combined_df[col] = "low_n"
        return combined_df

    mat        = bulk_wide.to_numpy(dtype=np.float32)
    feat_names = bulk_wide.index.tolist()
    samples    = bulk_wide.columns.tolist()

    # Count non-NaN values per junction (row) — only fit rows with enough data
    n_valid  = np.sum(~np.isnan(mat), axis=1)

    if not loo:
        # ------------------------------------------------------------------
        # Standard path: fit once per junction using all bulk samples
        # ------------------------------------------------------------------
        fit_mask    = n_valid >= n_threshold
        low_n_names = [feat_names[i] for i in range(len(feat_names)) if not fit_mask[i]]
        fit_names   = [feat_names[i] for i in range(len(feat_names)) if fit_mask[i]]

        if fit_names:
            beta_df = fit_beta_dist_chunk(
                mat[fit_mask], fit_names, PSI_rescale_factor, n_threshold, threads
            ).reset_index().rename(
                columns={"index": fit_id_col, "n": n_col, "alpha": alpha_col,
                         "beta_param": beta_col, "expected": expected_col})
        else:
            beta_df = pd.DataFrame(columns=[fit_id_col, n_col, alpha_col, beta_col, expected_col])

        if low_n_names:
            low_n_df = pd.DataFrame({
                fit_id_col:   low_n_names,
                n_col:        [int(n_valid[i]) for i in range(len(feat_names)) if not fit_mask[i]],
                alpha_col:    "low_n",
                beta_col:     "low_n",
                expected_col: "low_n",
            })
            beta_df = pd.concat([beta_df, low_n_df], ignore_index=True)

        combined_df = combined_df.merge(beta_df, on=fit_id_col, how="left")
        combined_df[n_col] = combined_df[n_col].fillna(0)
        for col in (alpha_col, beta_col, expected_col):
            combined_df[col] = combined_df[col].fillna("low_n")

    else:
        # ------------------------------------------------------------------
        # LOO path: for each bulk sample, refit excluding that sample, then
        # attach the LOO parameters only to that sample's rows in combined_df.
        #
        # bulk_wide rows = junctions/splice-sites  (feat_names)
        # bulk_wide cols = samples
        #
        # We build a long DataFrame of (fit_id, sample, alpha, beta, expected, n)
        # and merge it on (fit_id_col, sample).
        # ------------------------------------------------------------------
        loo_records = []   # list of dicts, one per (feature, sample)

        n_samples = len(samples)
        for col_i, samp in enumerate(samples):
            # LOO matrix: drop column col_i
            loo_cols = [j for j in range(n_samples) if j != col_i]
            loo_mat  = mat[:, loo_cols]  # shape: (n_features, n_samples-1)

            loo_n_valid = np.sum(~np.isnan(loo_mat), axis=1)
            fit_mask    = loo_n_valid >= n_threshold

            fit_names   = [feat_names[i] for i in range(len(feat_names)) if fit_mask[i]]
            fit_indices = [i                for i in range(len(feat_names)) if fit_mask[i]]

            if fit_names:
                fitted = fit_beta_dist_chunk(
                    loo_mat[fit_mask], fit_names,
                    PSI_rescale_factor, n_threshold, threads,
                )
                for feat, row in fitted.iterrows():
                    loo_records.append({
                        fit_id_col:   feat,
                        "sample":     samp,
                        n_col:        int(row["n"]),
                        alpha_col:    row["alpha"],
                        beta_col:     row["beta_param"],
                        expected_col: row["expected"],
                    })

            # Features that don't pass n_threshold even after LOO
            for i in range(len(feat_names)):
                if not fit_mask[i]:
                    loo_records.append({
                        fit_id_col:   feat_names[i],
                        "sample":     samp,
                        n_col:        int(loo_n_valid[i]),
                        alpha_col:    "low_n",
                        beta_col:     "low_n",
                        expected_col: "low_n",
                    })

        loo_df = pd.DataFrame(loo_records)

        # Merge on both fit_id_col and sample so each row gets its own LOO params
        combined_df = combined_df.merge(loo_df, on=[fit_id_col, "sample"], how="left")
        combined_df[n_col] = combined_df[n_col].fillna(0)
        for col in (alpha_col, beta_col, expected_col):
            combined_df[col] = combined_df[col].fillna("low_n")

        # Non-bulk rows (hap1/hap2) won't have matched in the LOO table (which
        # only covers bulk samples).  Fill them from the all-sample global fit
        # so haplotype rows still get distribution parameters for their tests.
        hap_mask = combined_df["phasing"].isin(["hap1", "hap2"])
        if hap_mask.any():
            fit_mask_global = n_valid >= n_threshold
            fit_names_global = [feat_names[i] for i in range(len(feat_names)) if fit_mask_global[i]]
            if fit_names_global:
                global_beta_df = fit_beta_dist_chunk(
                    mat[fit_mask_global], fit_names_global,
                    PSI_rescale_factor, n_threshold, threads,
                ).reset_index().rename(
                    columns={"index": fit_id_col, "n": f"_g_{n_col}",
                             "alpha": f"_g_{alpha_col}", "beta_param": f"_g_{beta_col}",
                             "expected": f"_g_{expected_col}"})
                combined_df = combined_df.merge(global_beta_df, on=fit_id_col, how="left")
                for col, gcol in ((alpha_col,    f"_g_{alpha_col}"),
                                  (beta_col,     f"_g_{beta_col}"),
                                  (expected_col, f"_g_{expected_col}"),
                                  (n_col,        f"_g_{n_col}")):
                    needs_fill = hap_mask & (
                        (combined_df[col] == "low_n") | combined_df[col].isna()
                    )
                    combined_df.loc[needs_fill, col] = combined_df.loc[needs_fill, gcol]
                drop_cols = [c for c in combined_df.columns if c.startswith("_g_")]
                combined_df = combined_df.drop(columns=drop_cols)

    alpha_vals   = combined_df[alpha_col]
    is_low_n     = alpha_vals == "low_n"
    is_error     = alpha_vals == "error"
    has_fit      = ~is_low_n & ~is_error
    coverage_arr = combined_df[coverage_col].to_numpy(dtype=float)
    has_cov      = coverage_arr >= coverage_threshold

    # Identify low-phased hap rows upfront (before any testing)
    lp_col_name = f"low_phased_{metric_col}"
    is_hap = combined_df["phasing"].isin(["hap1", "hap2"])
    if lp_col_name in combined_df.columns:
        is_low_phased = combined_df[lp_col_name].astype(bool)
    else:
        is_low_phased = pd.Series(False, index=combined_df.index)
    # low_phased only applies when the row also has coverage (otherwise low_coverage wins)
    is_low_phased_hap = is_hap & is_low_phased & has_cov

    import math
    # Use ceil(n * percentile) as index for conservative bounds:
    # e.g. n=107: ceil(107*0.01)=2 → take 3rd value (0-based index 2),
    # excluding the bottom 2; p99 takes the (n - ceil(n*0.01) - 1)th value.
    #
    # Under LOO, p1/p99 are also computed leave-one-out so the empirical range
    # doesn't include the sample being tested.
    if not loo:
        p1_vals  = np.full(len(feat_names), np.nan)
        p99_vals = np.full(len(feat_names), np.nan)
        for row_i in range(len(mat)):
            row  = mat[row_i]
            vals = np.sort(row[~np.isnan(row)])
            n    = len(vals)
            if n == 0:
                continue
            if n <= 10:
                p1_vals[row_i]  = vals[0]
                p99_vals[row_i] = vals[-1]
            else:
                k = math.ceil(n * 0.01)
                p1_vals[row_i]  = vals[k]
                p99_vals[row_i] = vals[n - 1 - k]
        p1_series  = pd.Series(p1_vals,  index=feat_names, name=p1_col)
        p99_series = pd.Series(p99_vals, index=feat_names, name=p99_col)
        combined_df = combined_df.merge(
            p1_series.reset_index().rename(columns={"index": fit_id_col}),
            on=fit_id_col, how="left")
        combined_df = combined_df.merge(
            p99_series.reset_index().rename(columns={"index": fit_id_col}),
            on=fit_id_col, how="left")
        for pc in (p1_col, p99_col):
            combined_df[pc] = combined_df[pc].where(combined_df[pc].notna(), other="low_n")
    else:
        # LOO p1/p99: for each (feature, sample) pair, compute empirical
        # percentiles from the remaining n-1 bulk values and merge on both keys.
        loo_pct_records = []
        for col_i, samp in enumerate(samples):
            loo_cols = [j for j in range(len(samples)) if j != col_i]
            loo_mat  = mat[:, loo_cols]
            for row_i, feat in enumerate(feat_names):
                row  = loo_mat[row_i]
                vals = np.sort(row[~np.isnan(row)])
                n    = len(vals)
                if n == 0:
                    p1v = np.nan; p99v = np.nan
                elif n <= 10:
                    p1v = float(vals[0]); p99v = float(vals[-1])
                else:
                    k = math.ceil(n * 0.01)
                    p1v = float(vals[k]); p99v = float(vals[n - 1 - k])
                loo_pct_records.append({
                    fit_id_col: feat,
                    "sample":   samp,
                    p1_col:     p1v,
                    p99_col:    p99v,
                })
        loo_pct_df = pd.DataFrame(loo_pct_records)
        combined_df = combined_df.merge(loo_pct_df, on=[fit_id_col, "sample"], how="left")

        # Hap rows: use global (all-sample) percentiles as fallback
        hap_mask = combined_df["phasing"].isin(["hap1", "hap2"])
        if hap_mask.any():
            global_p1  = np.full(len(feat_names), np.nan)
            global_p99 = np.full(len(feat_names), np.nan)
            for row_i in range(len(mat)):
                row  = mat[row_i]
                vals = np.sort(row[~np.isnan(row)])
                n    = len(vals)
                if n == 0: continue
                if n <= 10:
                    global_p1[row_i]  = vals[0]; global_p99[row_i] = vals[-1]
                else:
                    k = math.ceil(n * 0.01)
                    global_p1[row_i]  = vals[k]; global_p99[row_i] = vals[n - 1 - k]
            gp1_s  = pd.Series(global_p1,  index=feat_names, name=f"_g_{p1_col}")
            gp99_s = pd.Series(global_p99, index=feat_names, name=f"_g_{p99_col}")
            combined_df = combined_df.merge(
                gp1_s.reset_index().rename(columns={"index": fit_id_col}),
                on=fit_id_col, how="left")
            combined_df = combined_df.merge(
                gp99_s.reset_index().rename(columns={"index": fit_id_col}),
                on=fit_id_col, how="left")
            for pc, gpc in ((p1_col, f"_g_{p1_col}"), (p99_col, f"_g_{p99_col}")):
                needs_fill = hap_mask & combined_df[pc].isna()
                combined_df.loc[needs_fill, pc] = combined_df.loc[needs_fill, gpc]
            combined_df = combined_df.drop(
                columns=[c for c in combined_df.columns if c.startswith("_g_")])

        for pc in (p1_col, p99_col):
            combined_df[pc] = combined_df[pc].where(combined_df[pc].notna(), other="low_n")

    rescaled_v = pd.to_numeric(combined_df[rescaled_col], errors="coerce")
    p1_v       = pd.to_numeric(combined_df[p1_col],       errors="coerce")
    p99_v      = pd.to_numeric(combined_df[p99_col],      errors="coerce")
    delta_num  = np.where(
        rescaled_v > p99_v, rescaled_v - p99_v,
        np.where(rescaled_v < p1_v, rescaled_v - p1_v, 0.0)
    )
    delta_v = pd.Series(delta_num, index=combined_df.index).astype(object)
    delta_v[is_low_n]                       = "low_n"
    delta_v[is_error]                       = "error"
    delta_v[has_fit & ~has_cov]             = "low_coverage"
    delta_v[has_fit & is_low_phased_hap]    = "low_phased_coverage"
    combined_df[delta_col] = delta_v

    # Assign p_value sentinel strings in priority order:
    #   1. low_coverage  — row coverage < threshold (takes priority over everything)
    #   2. low_phased_coverage — hap row that fails phasing check for this metric
    #   3. low_n / error — distribution could not be fit
    #   4. actual test   — all other rows with a fit and sufficient coverage
    combined_df[p_col] = np.nan  # default; will be filled below
    combined_df.loc[is_low_n,                    p_col] = "low_n"
    combined_df.loc[is_error,                    p_col] = "error"
    combined_df.loc[has_fit & ~has_cov,          p_col] = "low_coverage"
    combined_df.loc[has_fit & is_low_phased_hap, p_col] = "low_phased_coverage"

    # Testable: has a fitted distribution, sufficient coverage, and is NOT low-phased
    testable = combined_df[has_fit & has_cov & ~is_low_phased_hap].copy()
    if len(testable) > 0:
        orig_index = testable.index
        tested = beta_binomial_test_chunk(
            testable.reset_index(drop=True),
            coverage_threshold, metric_col, usage_col, coverage_col, p_col,
            threads,
        )
        combined_df.loc[orig_index, p_col] = tested[p_col].values

    if is_ss_metric:
        combined_df = combined_df.drop(columns=["_ss_id"])

    return combined_df


def run_all_metrics(
    combined_df:        pd.DataFrame,
    coverage_threshold: int,
    PSI_rescale_factor: float,
    n_threshold:        int,
    threads:            int,
    has_ipa:            bool,
    no_ss_ir:           bool = False,
    loo:                bool = False,
) -> pd.DataFrame:
    metrics = [
        ("junction_PSI_approx", "rescaled_junction_PSI_approx", "junction_usage", "junction_coverage_approx", "junction"),
        ("junction_PSI",        "rescaled_junction_PSI",        "junction_usage", "junction_coverage",        "junction"),
        ("junction_full_IR_ratio", "rescaled_junction_full_IR_ratio", "junction_full_IR_count", "junction_coverage", "junction"),
    ]
    if not no_ss_ir:
        metrics.insert(2, ("5ss_IR_ratio", "rescaled_5ss_IR_ratio", "5ss_usage", "5ss_coverage", "junction"))
        metrics.insert(3, ("3ss_IR_ratio", "rescaled_3ss_IR_ratio", "3ss_usage", "3ss_coverage", "junction"))
    if has_ipa:
        metrics.append(
            ("junction_IPA_ratio", "rescaled_junction_IPA_ratio", "junction_IPA_count", "5ss_coverage", "junction")
        )
    for metric_col, rescaled_col, usage_col, coverage_col, id_col in metrics:
        if metric_col not in combined_df.columns:
            continue
        combined_df = _run_one_metric(
            combined_df, metric_col, rescaled_col, usage_col, coverage_col,
            id_col, coverage_threshold, PSI_rescale_factor, n_threshold, threads,
            loo=loo,
        )
    return combined_df


# ---------------------------------------------------------------------------
# Output column definitions
# ---------------------------------------------------------------------------

# NOTE: junction_type inserted after 3ss — silently skipped if absent (no --gtf)
_OUTPUT_COLS = [
    "sample", "gene", "gene_rank", "region", "phasing", "junction", "5ss", "3ss",
    "junction_type",   # only present when --gtf provided
    "junction_usage",
    "junction_read_diversity",
    "junction_coverage_approx",
    "junction_PSI_approx", "rescaled_junction_PSI_approx",
    "n_junction_PSI_approx", "alpha_junction_PSI_approx", "beta_junction_PSI_approx",
    "expected_junction_PSI_approx", "p1_junction_PSI_approx", "p99_junction_PSI_approx",
    "delta_junction_PSI_approx", "p_value_junction_PSI_approx", "padj_junction_PSI_approx",
    "junction_coverage",
    "junction_PSI", "rescaled_junction_PSI",
    "n_junction_PSI", "alpha_junction_PSI", "beta_junction_PSI",
    "expected_junction_PSI", "p1_junction_PSI", "p99_junction_PSI",
    "delta_junction_PSI", "p_value_junction_PSI", "padj_junction_PSI",
    "5ss_usage", "5ss_coverage",
    "5ss_IR_ratio", "rescaled_5ss_IR_ratio",
    "n_5ss_IR_ratio", "alpha_5ss_IR_ratio", "beta_5ss_IR_ratio",
    "expected_5ss_IR_ratio", "p1_5ss_IR_ratio", "p99_5ss_IR_ratio",
    "delta_5ss_IR_ratio", "p_value_5ss_IR_ratio", "padj_5ss_IR_ratio",
    "3ss_usage", "3ss_coverage",
    "3ss_IR_ratio", "rescaled_3ss_IR_ratio",
    "n_3ss_IR_ratio", "alpha_3ss_IR_ratio", "beta_3ss_IR_ratio",
    "expected_3ss_IR_ratio", "p1_3ss_IR_ratio", "p99_3ss_IR_ratio",
    "delta_3ss_IR_ratio", "p_value_3ss_IR_ratio", "padj_3ss_IR_ratio",
    "junction_full_IR_count",   # singular
    "junction_full_IR_ratio", "rescaled_junction_full_IR_ratio",
    "n_junction_full_IR_ratio", "alpha_junction_full_IR_ratio", "beta_junction_full_IR_ratio",
    "expected_junction_full_IR_ratio", "p1_junction_full_IR_ratio", "p99_junction_full_IR_ratio",
    "delta_junction_full_IR_ratio", "p_value_junction_full_IR_ratio", "padj_junction_full_IR_ratio",
    "junction_IPA_count",       # singular
    "junction_IPA_ratio", "rescaled_junction_IPA_ratio",
    "n_junction_IPA_ratio", "alpha_junction_IPA_ratio", "beta_junction_IPA_ratio",
    "expected_junction_IPA_ratio", "p1_junction_IPA_ratio", "p99_junction_IPA_ratio",
    "delta_junction_IPA_ratio", "p_value_junction_IPA_ratio", "padj_junction_IPA_ratio",
]


def select_output_columns(df: pd.DataFrame, *_) -> pd.DataFrame:
    """Select and reorder to canonical output schema. Missing columns silently skipped."""
    present = [c for c in _OUTPUT_COLS if c in df.columns]
    return df[present]


# ---------------------------------------------------------------------------
# GTF parsing — junctions for QC figures and junction_type classification
# ---------------------------------------------------------------------------

def parse_gtf_junctions(
    gtf_path: str,
    gene_names: List[str],
) -> Dict[str, Dict]:
    """
    Returns {gene: {"canonical_junctions": set, "all_junctions": set}}
    Junction strings are "chrom_ss1_ss2".
    """
    import gzip as _gz
    import re as _re

    gene_set = set(gene_names)

    def _attr(attr_str, key):
        m = _re.search(rf'{key}' + r'\s+"([^"]+)"', attr_str)
        return m.group(1) if m else None

    def _strip_ver(s):
        return _re.sub(r'\.\d+$', '', s) if s else s

    tx_canonical: set = set()
    gene_exons: Dict[str, Dict[str, List[Tuple[str, int, int]]]] = defaultdict(lambda: defaultdict(list))

    open_fn = _gz.open if gtf_path.endswith(".gz") else open

    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            feat = parts[2]
            if feat not in ("transcript", "exon"): continue

            chrom = parts[0]
            start = int(parts[3]) - 1   # 0-based
            end   = int(parts[4])        # half-open
            attrs = parts[8]

            gname = _attr(attrs, "gene_name") or _attr(attrs, "gene_symbol")
            gid   = _strip_ver(_attr(attrs, "gene_id"))
            gene_key = None
            if gname and gname in gene_set:
                gene_key = gname
            elif gid and gid in gene_set:
                gene_key = gid
            if gene_key is None:
                continue

            tx_id = _attr(attrs, "transcript_id")
            if tx_id is None: continue

            if feat == "transcript":
                if "Ensembl_canonical" in attrs:
                    tx_canonical.add(tx_id)
            elif feat == "exon":
                gene_exons[gene_key][tx_id].append((chrom, start, end))

    result: Dict[str, Dict] = {}
    for gene_key, tx_dict in gene_exons.items():
        canonical_jxns: set = set()
        all_jxns_set:   set = set()
        for tx_id, exons in tx_dict.items():
            if len(exons) < 2: continue
            exons_sorted = sorted(exons, key=lambda x: x[1])
            for i in range(len(exons_sorted) - 1):
                chrom_e, _, end_e   = exons_sorted[i]
                chrom_n, start_n, _ = exons_sorted[i + 1]
                if chrom_e != chrom_n: continue
                jxn = f"{chrom_e}_{end_e + 1}_{start_n}"
                all_jxns_set.add(jxn)
                if tx_id in tx_canonical:
                    canonical_jxns.add(jxn)
        result[gene_key] = {
            "canonical_junctions": canonical_jxns,
            "all_junctions":       all_jxns_set,
        }
    return result


def assign_junction_types(
    df: pd.DataFrame,
    gtf_junctions: Dict[str, Dict],
) -> pd.DataFrame:
    """
    Add a 'junction_type' column with values 'canonical', 'annotated', or 'novel'.
    - canonical : junction is in the Ensembl canonical transcript for that gene
    - annotated : junction is in any transcript for that gene (but not canonical)
    - novel     : junction is not in any annotated transcript

    Operates per-gene. Returns df with junction_type column added.
    """
    df = df.copy()
    jtype = pd.Series("novel", index=df.index, dtype=object)

    for gene, gdf in df.groupby("gene"):
        if gene not in gtf_junctions:
            continue  # leave as "novel"
        canonical_set = gtf_junctions[gene]["canonical_junctions"]
        annotated_set = gtf_junctions[gene]["all_junctions"]
        jxns = df.loc[gdf.index, "junction"]
        # annotated (not canonical)
        is_annotated = jxns.isin(annotated_set) & ~jxns.isin(canonical_set)
        # canonical
        is_canonical = jxns.isin(canonical_set)
        jtype.loc[gdf.index[is_annotated]] = "annotated"
        jtype.loc[gdf.index[is_canonical]] = "canonical"

    df["junction_type"] = jtype
    return df


# ---------------------------------------------------------------------------
# Event classification
# ---------------------------------------------------------------------------

_IR_IPA_EVENTS = frozenset(("5ss_IR", "3ss_IR", "full_IR", "IPA"))


def _classify_events_for_metric(
    metric_col:         str,
    outlier_set:        set,
    grp_coords:         Dict[Tuple, Dict[str, Tuple[int, int]]],
    strand_map:         Dict[str, str],
    padj_threshold:     float,
    delta_threshold:    float,
    positive_delta_set: set,
    delta_map:          Dict[Tuple, float],
) -> Dict[Tuple, List[str]]:
    """
    Classify splicing events for one metric.

    delta_map: {(sample, gene, phasing, jxn): delta_value} for all rows with
    a numeric delta.  Used to determine direction for alt_5ss/alt_3ss/exon
    classification.

    positive_delta_set: subset of outlier_set where delta >= +delta_threshold.
    For IR/IPA events (5ss_IR, 3ss_IR, full_IR, IPA) the row must be in
    positive_delta_set to receive the event label.

    Event rules:
      alt_3ss: query junction J is an outlier; there exists another junction J'
               in the same (sample, gene, phasing) group sharing J's 5ss but
               with a different 3ss, and delta(J') has the opposite sign to
               delta(J).  J' need not be an outlier itself.
      alt_5ss: symmetric — same 3ss, different 5ss, opposite-sign delta.
      exon_skipping / exon_inclusion:
               Three outlier junctions form a cassette triple:
                 long:  ss1a → ss2b
                 left:  ss1a → ss2a  (ss2a < ss2b)
                 right: ss1b → ss2b  (ss1b > ss1a)
               with ss2a == ss1b (the skipped exon boundaries meet).
               If long goes up (delta > 0) and both short go down → exon_skipping
               If long goes down (delta < 0) and both short go up  → exon_inclusion
    """
    allowed_events = _METRIC_EVENTS.get(metric_col, [])
    if not allowed_events:
        return {}

    outlier_by_grp: Dict[Tuple, set] = defaultdict(set)
    for key in outlier_set:
        sample, gene, phasing, jxn = key
        outlier_by_grp[(sample, gene, phasing)].add(jxn)

    event_map: Dict[Tuple, List[str]] = defaultdict(list)

    def _add(sample, gene, phasing, jxn, ev):
        key = (sample, gene, phasing, jxn)
        if ev not in event_map[key]:
            event_map[key].append(ev)

    def _delta(sample, gene, phasing, jxn):
        return delta_map.get((sample, gene, phasing, jxn), float("nan"))

    def _sign(d):
        """Return 1, -1, or 0."""
        if d > 0: return 1
        if d < 0: return -1
        return 0

    for grp_key, outlier_jxns in outlier_by_grp.items():
        sample, gene, phasing = grp_key
        strand = strand_map.get(gene, "+")
        coords = grp_coords[grp_key]  # ALL junctions for this group

        for jxn in outlier_jxns:
            if jxn not in coords:
                continue
            ss1, ss2 = coords[jxn]
            five_ss  = ss1 if strand == "+" else ss2
            three_ss = ss2 if strand == "+" else ss1

            row_key  = (sample, gene, phasing, jxn)
            is_pos_delta = row_key in positive_delta_set

            # IR/IPA events: only assign when delta is positive
            if "5ss_IR"  in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "5ss_IR")
            if "3ss_IR"  in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "3ss_IR")
            if "full_IR" in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "full_IR")
            if "IPA"     in allowed_events and is_pos_delta: _add(sample, gene, phasing, jxn, "IPA")

            if "alt_5ss" in allowed_events or "alt_3ss" in allowed_events or \
               "alt_5ss_approx" in allowed_events or "alt_3ss_approx" in allowed_events:
                d_jxn = _delta(sample, gene, phasing, jxn)
                s_jxn = _sign(d_jxn)
                alt_5ss_label = next((e for e in allowed_events if e.startswith("alt_5ss")), None)
                alt_3ss_label = next((e for e in allowed_events if e.startswith("alt_3ss")), None)
                # Partner must also be an outlier
                for partner_jxn in outlier_jxns:
                    if partner_jxn == jxn or partner_jxn not in coords:
                        continue
                    p_ss1, p_ss2 = coords[partner_jxn]
                    p_five  = p_ss1 if strand == "+" else p_ss2
                    p_three = p_ss2 if strand == "+" else p_ss1
                    d_partner = _delta(sample, gene, phasing, partner_jxn)
                    s_partner = _sign(d_partner)
                    opposite = (s_jxn != 0 and s_partner != 0 and s_jxn != s_partner)

                    if alt_3ss_label:
                        if p_five == five_ss and p_three != three_ss and opposite:
                            _add(sample, gene, phasing, jxn,         alt_3ss_label)
                            _add(sample, gene, phasing, partner_jxn, alt_3ss_label)

                    if alt_5ss_label:
                        if p_three == three_ss and p_five != five_ss and opposite:
                            _add(sample, gene, phasing, jxn,         alt_5ss_label)
                            _add(sample, gene, phasing, partner_jxn, alt_5ss_label)

            skip_label    = next((e for e in allowed_events if e.startswith("exon_skipping")),  None)
            incl_label    = next((e for e in allowed_events if e.startswith("exon_inclusion")), None)
            if skip_label or incl_label:
                d_long = _delta(sample, gene, phasing, jxn)
                s_long = _sign(d_long)
                if s_long == 0:
                    continue
                # jxn is the long junction (ss1a → ss2b).
                # Find outlier short junctions sharing ss1 (left) and ss2 (right).
                # No requirement that left's ss2 == right's ss1 — the skipped exon
                # can have any size; all that matters is the outer coordinates match.
                left_candidates = [
                    (j, s1, s2) for j, (s1, s2) in coords.items()
                    if j in outlier_jxns and j != jxn
                    and s1 == ss1 and s2 < ss2
                ]
                right_candidates = [
                    (j, s1, s2) for j, (s1, s2) in coords.items()
                    if j in outlier_jxns and j != jxn
                    and s2 == ss2 and s1 > ss1
                ]
                for jl, ls1, ls2 in left_candidates:
                    for jr, rs1, rs2 in right_candidates:
                        d_left  = _delta(sample, gene, phasing, jl)
                        d_right = _delta(sample, gene, phasing, jr)
                        s_left  = _sign(d_left)
                        s_right = _sign(d_right)
                        # Both short junctions must go opposite to the long
                        if s_left == 0 or s_right == 0:
                            continue
                        if s_left == s_long or s_right == s_long:
                            continue
                        if skip_label and s_long > 0:
                            for j in (jxn, jl, jr):
                                _add(sample, gene, phasing, j, skip_label)
                        if incl_label and s_long < 0:
                            for j in (jxn, jl, jr):
                                _add(sample, gene, phasing, j, incl_label)

    return dict(event_map)


def classify_all_events(
    sig_df:              pd.DataFrame,
    final_df:            pd.DataFrame,
    strand_map:          Dict[str, str],
    padj_threshold:      float,
    delta_threshold:     float,
    threads:             int,
    has_ipa:             bool,
    p_cols:              List[str],
    unreliable_hap_outliers: Dict[str, set],
) -> pd.DataFrame:
    if sig_df.empty:
        sig_df = sig_df.copy()
        sig_df["event_type"] = ""
        return sig_df

    computed_metrics = [p.replace("p_value_", "") for p in p_cols]
    active_metrics   = [mc for mc in computed_metrics if _METRIC_EVENTS.get(mc)]
    if not active_metrics:
        sig_df = sig_df.copy(); sig_df["event_type"] = "other"; return sig_df

    # Split grp_coords by gene to keep per-task payloads small
    all_genes  = sorted(final_df["gene"].unique())
    n_genes    = len(all_genes)
    batch_size = max(1, n_genes // (threads * 4)) if threads > 1 else n_genes
    gene_batches = [all_genes[i:i + batch_size]
                    for i in range(0, n_genes, batch_size)]

    # grp_coords_by_gene: {gene: {(sample, gene, phasing): {jxn: (ss1, ss2)}}}
    grp_coords_by_gene: Dict[str, Dict] = {}
    for g in all_genes:
        grp_coords_by_gene[g] = {}
    for _, row in final_df[["sample", "gene", "phasing", "junction"]].iterrows():
        g   = row["gene"]
        key = (row["sample"], g, row["phasing"])
        jxn = row["junction"]
        parts = jxn.split("_")
        grp_coords_by_gene[g][key] = grp_coords_by_gene[g].get(key, {})
        grp_coords_by_gene[g][key][jxn] = (int(parts[-2]), int(parts[-1]))

    # Per-metric sets / maps (full, then sliced per batch at dispatch time)
    metric_outlier_sets:   Dict[str, set] = {}
    metric_pos_delta_sets: Dict[str, set] = {}
    metric_delta_maps:     Dict[str, Dict[Tuple, float]] = {}

    _ir_ipa_metrics = frozenset((
        "5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio", "junction_IPA_ratio"
    ))
    _has_jxn_type = "junction_type" in sig_df.columns

    for mc in computed_metrics:
        ocol    = f"outlier_{mc}"
        delta_c = f"delta_{mc}"
        if ocol not in sig_df.columns:
            metric_outlier_sets[mc]   = set()
            metric_pos_delta_sets[mc] = set()
        else:
            mask = sig_df[ocol].astype(bool)

            # For IR/IPA metrics: additionally require positive delta and,
            # if --gtf was provided, canonical or annotated junction_type.
            if mc in _ir_ipa_metrics:
                if delta_c in sig_df.columns:
                    dv   = pd.to_numeric(sig_df[delta_c], errors="coerce")
                    mask = mask & dv.ge(delta_threshold)
                if _has_jxn_type:
                    mask = mask & sig_df["junction_type"].isin(["canonical", "annotated"])

            metric_outlier_sets[mc] = set(zip(
                sig_df.loc[mask, "sample"], sig_df.loc[mask, "gene"],
                sig_df.loc[mask, "phasing"], sig_df.loc[mask, "junction"],
            ))
            if delta_c in sig_df.columns:
                dv       = pd.to_numeric(sig_df[delta_c], errors="coerce")
                pos_mask = mask & dv.ge(delta_threshold)
            else:
                pos_mask = mask
            metric_pos_delta_sets[mc] = set(zip(
                sig_df.loc[pos_mask, "sample"], sig_df.loc[pos_mask, "gene"],
                sig_df.loc[pos_mask, "phasing"], sig_df.loc[pos_mask, "junction"],
            ))
        if delta_c not in final_df.columns:
            metric_delta_maps[mc] = {}
        else:
            dv    = pd.to_numeric(final_df[delta_c], errors="coerce")
            valid = dv.notna()
            metric_delta_maps[mc] = dict(zip(
                zip(final_df.loc[valid, "sample"], final_df.loc[valid, "gene"],
                    final_df.loc[valid, "phasing"], final_df.loc[valid, "junction"]),
                dv[valid],
            ))

    def _slice_batch(gene_batch, mc):
        gene_set = set(gene_batch)
        gc = {}
        for g in gene_batch:
            gc.update(grp_coords_by_gene.get(g, {}))
        os_ = {k for k in metric_outlier_sets.get(mc, set())   if k[1] in gene_set}
        ps_ = {k for k in metric_pos_delta_sets.get(mc, set()) if k[1] in gene_set}
        dm_ = {k: v for k, v in metric_delta_maps.get(mc, {}).items() if k[1] in gene_set}
        return gc, os_, ps_, dm_

    all_event_maps: List[Dict[Tuple, List[str]]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futs = {}
        for mc in active_metrics:
            for batch in gene_batches:
                gc, os_, ps_, dm_ = _slice_batch(batch, mc)
                if not os_:
                    continue
                fut = ex.submit(
                    _classify_events_for_metric, mc,
                    os_, gc, strand_map, padj_threshold, delta_threshold, ps_, dm_,
                )
                futs[fut] = (mc, batch[0])
        for fut in concurrent.futures.as_completed(futs):
            mc, g0 = futs[fut]
            try:
                all_event_maps.append(fut.result())
            except Exception as e:
                print(f"[WARNING] Classification error for {mc} (batch @{g0}): {e}")
                traceback.print_exc()

    merged: Dict[Tuple, set] = defaultdict(set)
    for event_map in all_event_maps:
        for key, events in event_map.items():
            merged[key].update(events)

    sig_df = sig_df.copy()
    sig_records = sig_df[["sample", "gene", "phasing", "junction"]].to_dict("records")
    event_strs = []
    for rec in sig_records:
        key = (rec["sample"], rec["gene"], rec["phasing"], rec["junction"])
        evs = merged.get(key, set())
        event_strs.append(",".join(sorted(evs)) if evs else "none")
    sig_df["event_type"] = event_strs
    return sig_df


# ---------------------------------------------------------------------------
# process_one_region
# ---------------------------------------------------------------------------

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


def process_one_region(
    gene: str,
    region: str,
    region_df: pd.DataFrame,
    approx_only: bool,
    coverage_threshold: int,
    PSI_rescale_factor: float,
    n_threshold: int,
    phasing_threshold: float,
    threads: int,
    strand: str = "+",
    include_monoexonic: bool = False,
    min_jxn_reads: int = 20,
    genome_path: Optional[str] = None,
    alu: Optional[dict] = None,
    no_ss_ir: bool = False,
    loo: bool = False,
) -> Optional[pd.DataFrame]:
    print(f"\n{'='*70}")
    print(f"  Gene: {gene}   Region: {region}")
    print(f"{'='*70}")

    bulk_rows = region_df[region_df["bulk"].notna()]

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as pool:

        if approx_only:
            print(f"  [1/2] PSI_approx + junction discovery ({len(region_df)} samples) ...")
            t0 = time.time()
            sample_dfs = []
            all_jxns_set: set = set()
            futures = {}
            sorted_rows = sorted(region_df.itertuples(index=False),
                                 key=lambda r: _row_bam_size(pd.Series(r._asdict())),
                                 reverse=True)
            for row in sorted_rows:
                row = pd.Series(row._asdict())
                sname = row["sample"]
                bulk  = row["bulk"] if pd.notna(row.get("bulk","")) else None
                hap1  = row["hap1"] if pd.notna(row.get("hap1","")) else None
                hap2  = row["hap2"] if pd.notna(row.get("hap2","")) else None
                futures[pool.submit(
                    process_sample, sname, bulk, hap1, hap2,
                    region, gene, [], True, coverage_threshold,
                    phasing_threshold, PSI_rescale_factor, strand,
                    include_monoexonic,
                )] = sname
            _lines: List[str] = []
            for fut in concurrent.futures.as_completed(futures):
                sname = futures[fut]
                try:
                    sdf, elapsed = fut.result()
                    _lines.append(f"       {sname}: {elapsed:.2f}s")
                    if len(sdf):
                        sample_dfs.append(sdf)
                        all_jxns_set.update(sdf["junction"].unique())
                except Exception as e:
                    print(f"  [ERROR] {sname}: {e}"); traceback.print_exc()
                    print(f"  [WARNING] Skipping {gene}."); return None, 0.0
            if not sample_dfs:
                print("  No data. Skipping."); return None, 0.0
            for _line in sorted(_lines): print(_line)
            combined = pd.concat(sorted(sample_dfs, key=lambda d: d["sample"].iloc[0]), ignore_index=True)
            bulk_max: Dict[str, int] = defaultdict(int)
            for sdf_b in sample_dfs:
                for _, r in sdf_b[sdf_b["phasing"]=="bulk"].iterrows():
                    jxn = r["junction"]; u = int(r["junction_usage"])
                    if u > bulk_max[jxn]: bulk_max[jxn] = u
            kept_jxns = {j for j in all_jxns_set if bulk_max.get(j, 0) >= min_jxn_reads}
            n_dropped = len(all_jxns_set) - len(kept_jxns)
            if n_dropped > 0:
                combined = combined[combined["junction"].isin(kept_jxns)]
            print(f"       → {len(all_jxns_set)} junctions found, {len(kept_jxns)} kept "
                  f"(≥{min_jxn_reads}), {n_dropped} dropped, "
                  f"{len(combined)} rows ({time.time()-t0:.2f}s)")

        else:
            print(f"  [1/3] Junction discovery (bulk BAMs, {len(bulk_rows)} samples) ...")
            t0 = time.time()
            jxn_count_union: Dict[Tuple[int,int], int] = defaultdict(int)
            # Also save per-sample jxn_raw to reuse in Step 2 (avoids re-walking bulk BAMs)
            sample_jxn_raw: Dict[str, Dict[Tuple[int,int], int]] = {}
            disc_futures = {
                pool.submit(_discover_junctions_worker,
                            row["sample"], row["bulk"], region, include_monoexonic): row["sample"]
                for _, row in bulk_rows.iterrows()
            }
            for fut in concurrent.futures.as_completed(disc_futures):
                try:
                    sname, jxn_raw, _ = fut.result()
                    sample_jxn_raw[sname] = jxn_raw
                    for coord, cnt in jxn_raw.items():
                        jxn_count_union[coord] = max(jxn_count_union[coord], cnt)
                except Exception as e:
                    print(f"  [ERROR] Discovery: {e}"); traceback.print_exc()

            chrom = region.split(":")[0]
            all_jxns = sorted(f"{chrom}_{s}_{e}"
                               for (s,e), cnt in jxn_count_union.items()
                               if cnt >= min_jxn_reads)
            n_tot = len(jxn_count_union); n_kept = len(all_jxns)
            print(f"       → {n_tot} found, {n_kept} kept "
                  f"(with ≥{min_jxn_reads} reads in ≥1 sample) "
                  f"({time.time()-t0:.2f}s)")
            if not all_jxns:
                print("  No junctions passed filter. Skipping."); return None, 0.0

            print(f"  [2/3] Calculating junction usage ({len(region_df)} samples) ...")
            t0 = time.time()
            sample_dfs = []
            psi_futures = {}
            sorted_rows_psi = sorted(region_df.itertuples(index=False),
                                     key=lambda r: _row_bam_size(pd.Series(r._asdict())),
                                     reverse=True)
            for row in sorted_rows_psi:
                row = pd.Series(row._asdict())
                sname = row["sample"]
                bulk  = row["bulk"] if pd.notna(row.get("bulk","")) else None
                hap1  = row["hap1"] if pd.notna(row.get("hap1","")) else None
                hap2  = row["hap2"] if pd.notna(row.get("hap2","")) else None
                # Pass bulk jxn_raw from Step 1 to avoid re-walking the bulk BAM
                bulk_jxn_raw = sample_jxn_raw.get(sname)
                psi_futures[pool.submit(
                    process_sample, sname, bulk, hap1, hap2,
                    region, gene, all_jxns, approx_only, coverage_threshold,
                    phasing_threshold, PSI_rescale_factor, strand,
                    include_monoexonic, genome_path, alu, no_ss_ir,
                    bulk_jxn_raw,
                )] = sname
            _psi_lines: List[str] = []
            _psi_t0 = time.time()
            for fut in concurrent.futures.as_completed(psi_futures):
                sname = psi_futures[fut]
                try:
                    sdf, elapsed = fut.result()
                    _psi_lines.append(f"       {sname}: {elapsed:.2f}s")
                    if len(sdf): sample_dfs.append(sdf)
                except Exception as e:
                    print(f"  [ERROR] {sname}: {e}"); traceback.print_exc()
                    print(f"  [WARNING] Skipping {gene}."); return None, 0.0
            if not sample_dfs:
                print("  No data. Skipping."); return None, 0.0
            combined = pd.concat(sorted(sample_dfs, key=lambda d: d["sample"].iloc[0]), ignore_index=True)
            _psi_wall = time.time() - _psi_t0
            for _line in sorted(_psi_lines): print(_line)
            print(f"       → {len(combined):,} rows (junctions × samples × phasings) "
                  f"({_psi_wall:.2f}s)")

    step = "[3/3]" if not approx_only else "[2/2]"
    print(f"  {step} Beta fitting + tests ...")
    t0 = time.time()
    has_ipa = (not approx_only) and genome_path is not None
    combined = combined.sort_values(
        ["junction", "sample", "phasing"], ignore_index=True
    )
    if approx_only:
        combined = _run_one_metric(
            combined, "junction_PSI_approx", "rescaled_junction_PSI_approx",
            "junction_usage", "junction_coverage_approx", "junction",
            coverage_threshold, PSI_rescale_factor, n_threshold, threads,
            loo=loo,
        )
    else:
        combined = run_all_metrics(combined, coverage_threshold, PSI_rescale_factor,
                                   n_threshold, threads, has_ipa, no_ss_ir,
                                   loo=loo)

    p_cols_present = [c for c in combined.columns if c.startswith("p_value_")]
    n_fit = 0; n_tests = 0
    bulk_combined = combined[combined["phasing"] == "bulk"]
    for p_col in p_cols_present:
        a_col = p_col.replace("p_value_", "alpha_")
        if a_col not in combined.columns: continue
        is_ss = any(x in p_col for x in ("5ss_IR_ratio", "3ss_IR_ratio"))
        if is_ss:
            ss_col = "5ss" if "5ss" in p_col else "3ss"
            if ss_col in bulk_combined.columns:
                n_fit += int(
                    bulk_combined[bulk_combined[a_col].apply(
                        lambda x: x not in ("low_n", "error")
                    )][ss_col].nunique()
                )
        else:
            n_fit += int(
                bulk_combined[bulk_combined[a_col].apply(
                    lambda x: x not in ("low_n", "error")
                )]["junction"].nunique()
            )
        n_tests += int(pd.to_numeric(combined[p_col], errors="coerce").notna().sum())
    print(f"       → fit {n_fit:,} distributions and performed {n_tests:,} tests "
          f"({time.time()-t0:.2f}s)")
    return combined, 0.0


# ---------------------------------------------------------------------------
# QC figures
# ---------------------------------------------------------------------------

# QC figure definitions: each entry is
#   (coverage_col, file_suffix, [bar_metrics])
# The file will be: {prefix}_qc_{file_suffix}_{subset}.pdf
_QC_FIGURES = [
    ("junction_coverage_approx", "junction_coverage_approx", ["junction_PSI_approx"]),
    ("junction_coverage",        "junction_coverage",        ["junction_PSI", "junction_full_IR_ratio"]),
    ("5ss_coverage",             "5ss_coverage",             ["5ss_IR_ratio", "junction_IPA_ratio"]),
    ("3ss_coverage",             "3ss_coverage",             ["3ss_IR_ratio"]),
]
# Coverage columns that belong to splice-site (rather than junction) metrics
_SS_COVERAGE_COLS = frozenset(("5ss_coverage", "3ss_coverage"))


def make_qc_figure(
    final_df: pd.DataFrame,
    gtf_junctions: Dict[str, Dict],
    cov_col: str,
    companions: List[str],
    file_suffix: str,
    jxn_subset: str,
    coverage_threshold: int,
    out_pdf: str,
    has_ipa: bool,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        print("[WARNING] matplotlib not available; skipping QC figures.")
        return

    color = "#d95d5b" if jxn_subset == "canonical" else "#4c8fca"
    cmap  = mcolors.LinearSegmentedColormap.from_list("c", ["white", color])

    companions = [m for m in companions
                  if has_ipa or m != "junction_IPA_ratio"]

    df = final_df  # already bulk-only, pre-sliced by caller
    genes = [g for g in gtf_junctions if g in df["gene"].unique()]
    if not genes:
        print(f"[WARNING] No overlapping genes for {file_suffix} ({jxn_subset}). Skipping.")
        return
    samples = sorted(df["sample"].unique())
    jxn_key = "canonical_junctions" if jxn_subset == "canonical" else "all_junctions"
    genes = [g for g in genes if len(gtf_junctions[g][jxn_key]) > 0]
    if not genes:
        print(f"[WARNING] No genes with {jxn_subset} junctions for {file_suffix}. Skipping.")
        return

    n_g = len(genes); n_s = len(samples)
    cov_mat = np.full((n_g, n_s), np.nan)
    cnt_mat = np.full((n_g, n_s), np.nan)
    fit_mat = np.full((n_g, len(companions)), np.nan)
    jxn_n   = []
    sample_idx = {s: i for i, s in enumerate(samples)}
    gene_idx   = {g: i for i, g in enumerate(genes)}

    for gene in genes:
        gi      = gene_idx[gene]
        jxn_set = gtf_junctions[gene][jxn_key]
        jxn_n.append(len(jxn_set))
        gdf     = df[df["gene"] == gene]
        for ci, cm in enumerate(companions):
            acol = f"alpha_{cm}"
            if acol not in gdf.columns:
                fit_mat[gi, ci] = 0; continue
            sub = gdf[gdf["junction"].isin(jxn_set)].drop_duplicates("junction")
            if sub.empty:
                fit_mat[gi, ci] = 0
            else:
                n_fitted = sub[acol].apply(lambda x: x not in ("low_n", "error", None)).sum()
                fit_mat[gi, ci] = float(n_fitted) / len(jxn_set) if jxn_set else np.nan
        for sample in samples:
            si = sample_idx[sample]
            sub = gdf[(gdf["sample"] == sample) & gdf["junction"].isin(jxn_set)]
            if not jxn_set or sub.empty:
                cov_mat[gi, si] = 0.0; cnt_mat[gi, si] = 0.0; continue
            cov = pd.to_numeric(sub[cov_col], errors="coerce")
            n   = int((cov >= coverage_threshold).sum())
            cov_mat[gi, si] = float(n) / len(jxn_set)
            cnt_mat[gi, si] = float(n)

    gene_order   = np.lexsort((-np.nanmean(cnt_mat, axis=1),  -np.nanmean(cov_mat, axis=1)))
    sample_order = np.lexsort((-np.nanmean(cnt_mat, axis=0),  -np.nanmean(cov_mat, axis=0)))
    cov_mat_s = cov_mat[np.ix_(gene_order, sample_order)]
    fit_mat_s = fit_mat[gene_order]
    genes_s   = [genes[i]  for i in gene_order]
    jxn_n_s   = [jxn_n[i]  for i in gene_order]
    ylabels   = [f"{g}  (n={n})" for g, n in zip(genes_s, jxn_n_s)]

    heat_w = min(max(2.0, n_s * 0.05), 10.0)
    bar_w  = 1.5
    n_bars = len(companions)
    fig_w  = heat_w + bar_w * n_bars + 0.8
    fig_h  = max(2.0, n_g * 0.18 + 1.2)
    width_ratios = [heat_w] + [bar_w] * n_bars
    fig, axes = plt.subplots(1, 1 + n_bars, figsize=(fig_w, fig_h),
                              gridspec_kw={"width_ratios": width_ratios, "wspace": 0.05})
    if 1 + n_bars == 1: axes = [axes]

    ax0 = axes[0]
    feat_word = "splice sites" if cov_col in _SS_COVERAGE_COLS else "junctions"
    im = ax0.imshow(cov_mat_s, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ax0.set_title(f"Proportion of {feat_word} with {cov_col} ≥ {coverage_threshold}", fontsize=8, pad=3)
    ax0.set_xticks([]); ax0.set_yticks(range(n_g))
    ax0.set_yticklabels(ylabels, fontsize=6)
    fig.colorbar(im, ax=ax0, fraction=0.03, pad=0.01)

    y_pos = np.arange(n_g)
    for ci, cm in enumerate(companions):
        ax = axes[ci + 1]
        ax.barh(y_pos, fit_mat_s[:, ci], 0.6, color=color, alpha=0.85)
        ax.set_xlim(0, 1)
        ax.set_title(f"Proportion of {feat_word} with\n{cm} modeled", fontsize=7, pad=3)
        ax.set_yticks(y_pos); ax.set_yticklabels([])
        ax.tick_params(left=False); ax.spines["left"].set_visible(False)
        ax.invert_yaxis()
    for ax in axes:
        ax.set_ylim(n_g - 0.5, -0.5)

    fig.suptitle(f"{cov_col} — {jxn_subset.capitalize()} junctions", fontsize=9, y=1.01)
    t_qc_fig = time.time()
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"  QC figure → {out_pdf} ({time.time()-t_qc_fig:.2f}s)")


# ---------------------------------------------------------------------------
# Outlier heatmap + box plots
# ---------------------------------------------------------------------------

def make_outlier_heatmap(
    sig_df: pd.DataFrame,
    metric_col: str,
    delta_col: str,
    padj_col: str,
    padj_threshold: float,
    delta_threshold: float,
    out_pdf: str,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        return

    padj_v  = pd.to_numeric(sig_df[padj_col],  errors="coerce")
    delta_v = pd.to_numeric(sig_df[delta_col], errors="coerce")
    mask    = padj_v.le(padj_threshold) & delta_v.abs().ge(delta_threshold)
    df      = sig_df[mask].copy()
    df["_abs_delta"] = delta_v[mask].abs()
    if df.empty: return

    agg = (df.groupby(["gene", "sample"])["_abs_delta"]
             .max().reset_index().rename(columns={"_abs_delta": "max_delta"}))
    genes   = sorted(agg["gene"].unique())
    samples = sorted(agg["sample"].unique())
    mat = np.full((len(genes), len(samples)), np.nan)
    g_idx = {g: i for i, g in enumerate(genes)}
    s_idx = {s: i for i, s in enumerate(samples)}
    for _, row in agg.iterrows():
        mat[g_idx[row["gene"]], s_idx[row["sample"]]] = row["max_delta"]

    gene_order   = np.lexsort((-np.nanmean(mat, axis=1), -(~np.isnan(mat)).sum(axis=1)))
    sample_order = np.lexsort((-np.nanmean(mat, axis=0), -(~np.isnan(mat)).sum(axis=0)))
    mat_s   = mat[np.ix_(gene_order, sample_order)]
    genes_s = [genes[i] for i in gene_order]
    n_g = len(genes_s); n_s = len(samples)

    cmap = mcolors.LinearSegmentedColormap.from_list("wd", ["#ffffff", "#912321"])
    cmap.set_bad(color="#f2f3f4")
    heat_w = min(max(2.0, n_s * 0.05), 12.0)
    fig_h  = max(2.0, n_g * 0.18 + 1.2)
    fig, ax = plt.subplots(figsize=(heat_w, fig_h))
    im = ax.imshow(mat_s, aspect="auto", cmap=cmap,
                   vmin=0, vmax=np.nanmax(mat_s), interpolation="nearest")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01, label=f"|delta {metric_col}|")
    ax.set_yticks(range(n_g)); ax.set_yticklabels(genes_s, fontsize=6)
    ax.set_xticks([])
    ax.set_title(f"Outlier heatmap: {metric_col}  (padj<={padj_threshold}, |delta|>={delta_threshold})",
                 fontsize=9)
    fig.tight_layout()
    t_heat = time.time()
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"  Outlier heatmap → {out_pdf} ({time.time()-t_heat:.2f}s)")


def make_hit_boxplots(
    fmt_updated: pd.DataFrame,
    outlier_map: dict,
    metric_col: str,
    rescaled_col: str,
    outdir: str,
    prefix_name: str,
    tmp_dir: str,
    padj_threshold: float = 0.01,
    delta_threshold: float = 0.1,
    is_ir_ipa: bool = False,
    has_jxn_type_filter: bool = False,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        return

    gene_jxn_map = outlier_map.get(metric_col, {})
    if not gene_jxn_map: return

    ind_dir = os.path.join(tmp_dir, metric_col)
    os.makedirs(ind_dir, exist_ok=True)
    out_pdf = os.path.join(outdir, f"{prefix_name}_boxplots_{metric_col}.pdf")
    n_hits = sum(len(jxns) for jxns in gene_jxn_map.values())
    n_written = 0

    t_box = time.time()
    with PdfPages(out_pdf) as pdf:
        # --- Title page ---
        fig_t, ax_t = plt.subplots(figsize=(4.5, 4.0))
        ax_t.axis("off")
        if is_ir_ipa:
            delta_str = f">= {delta_threshold}"
            jt_line   = "\njunction_type: canonical, annotated" if has_jxn_type_filter else ""
        else:
            delta_str = f"|{delta_threshold}|"
            jt_line   = ""
        title_text = (
            f"Metric: {metric_col}\n"
            f"padj threshold: <= {padj_threshold}\n"
            f"delta threshold: {delta_str}"
            f"{jt_line}"
        )
        ax_t.text(0.5, 0.5, title_text, transform=ax_t.transAxes,
                  fontsize=12, va="center", ha="center",
                  bbox=dict(boxstyle="round,pad=0.6", facecolor="#ffffff", edgecolor="#0068a9"))
        fig_t.tight_layout()
        pdf.savefig(fig_t, bbox_inches="tight")
        plt.close(fig_t)
        for gene in sorted(gene_jxn_map):
            for junction in sorted(gene_jxn_map[gene]):
                outlier_samples = gene_jxn_map[gene][junction]
                jxn_df = fmt_updated[
                    (fmt_updated["gene"] == gene) &
                    (fmt_updated["junction"] == junction)
                ]
                all_bulk   = jxn_df[jxn_df["phasing"] == "bulk"]
                bulk_vals  = pd.to_numeric(all_bulk[rescaled_col], errors="coerce").dropna()
                if len(bulk_vals) == 0: continue

                safe_jxn  = junction.replace(":", "_").replace("/", "_")
                out_path  = os.path.join(ind_dir, f"{gene}_{safe_jxn}.pdf")
                try:
                    fig, ax = plt.subplots(figsize=(4.5, 4.0))
                    ax.boxplot(bulk_vals.values, positions=[0], widths=0.4,
                               patch_artist=True, showfliers=False,
                               boxprops=dict(facecolor="#4c8fca", color="#91c4e9"),
                               medianprops=dict(color="#d95d5b", linewidth=2),
                               whiskerprops=dict(color="#91c4e9"),
                               capprops=dict(color="#91c4e9"))
                    non_out_vals = pd.to_numeric(
                        all_bulk[~all_bulk["sample"].isin(outlier_samples)][rescaled_col],
                        errors="coerce").dropna()
                    rng = np.random.default_rng(seed=42)
                    ax.scatter(rng.uniform(-0.18, 0.18, size=len(non_out_vals)),
                               non_out_vals.values, color="#91c4e9", alpha=0.35, s=10, zorder=3, linewidths=0)
                    colors = {"bulk": "#d95d5b", "hap1": "#ea9a9c", "hap2": "#ea9a9c"}
                    for si, sample in enumerate(sorted(outlier_samples)):
                        x = 0.08 + si * 0.06
                        sample_rows = jxn_df[jxn_df["sample"] == sample]
                        # Only plot haplotype values if BOTH hap1 and hap2 are available
                        hap1_row = sample_rows[sample_rows["phasing"] == "hap1"]
                        hap2_row = sample_rows[sample_rows["phasing"] == "hap2"]
                        hap1_val = pd.to_numeric(hap1_row[rescaled_col], errors="coerce").iloc[0] if not hap1_row.empty else float("nan")
                        hap2_val = pd.to_numeric(hap2_row[rescaled_col], errors="coerce").iloc[0] if not hap2_row.empty else float("nan")
                        plot_haps = pd.notna(hap1_val) and pd.notna(hap2_val)
                        for phasing in (("bulk", "hap1", "hap2") if plot_haps else ("bulk",)):
                            ph_row = sample_rows[sample_rows["phasing"] == phasing]
                            if ph_row.empty: continue
                            val = pd.to_numeric(ph_row[rescaled_col], errors="coerce").iloc[0]
                            if pd.isna(val): continue
                            ax.scatter(x, val, color=colors[phasing],
                                       s=40 if phasing == "bulk" else 22, zorder=5)
                            label = sample if phasing == "bulk" else f"{sample} ({phasing})"
                            ax.annotate(label, xy=(x, val), xytext=(x + 0.03, val),
                                        fontsize=5, va="center", ha="left", color=colors[phasing],
                                        arrowprops=dict(arrowstyle="-", color=colors[phasing], lw=0.5))
                    ax.set_xlim(-0.6, max(0.8, 0.08 + len(outlier_samples) * 0.06 + 0.3))
                    ax.set_xticks([]); ax.set_ylabel(f"rescaled {metric_col}", fontsize=8)
                    ax.set_ylim(-0.05, 1.05)
                    ax.set_title(f"{gene}: {junction}\n{metric_col}", fontsize=8)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    with PdfPages(out_path) as ind_pdf:
                        ind_pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig); n_written += 1
                except Exception as e:
                    print(f"  [WARNING] Box plot failed for {gene} {junction}: {e}")
                    plt.close("all")

    print(f"  Box plots → {out_pdf}  ({n_written}/{n_hits} pages, {time.time()-t_box:.2f}s)")


# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

def main() -> None:
    print("\n" + "*"*80)
    print("  Splice Junction Outlier Detection")
    print("*"*80 + "\n")

    args        = parse_args()
    approx_only = args.approx

    prefix      = args.outprefix.rstrip("/")
    outdir      = os.path.dirname(os.path.abspath(prefix))
    prefix_name = os.path.basename(prefix)
    base_results = os.path.join(outdir, f"{prefix_name}_results")
    qc_dir       = os.path.join(outdir, f"{prefix_name}_qc")
    tmp_dir      = os.path.join(outdir, f"{prefix_name}_tmp")

    # Parse threshold combinations — empty list means no outlier identification
    threshold_pairs: List[Tuple[float, float]] = []
    for tok in args.thresholds:
        try:
            p_str, d_str = tok.split(":")
            threshold_pairs.append((float(p_str), float(d_str)))
        except Exception:
            raise ValueError(f"Invalid threshold format '{tok}'. Expected padj:delta, e.g. 0.01:0.1")
    if not threshold_pairs:
        print("[INFO] No --thresholds provided; will write per-gene TSVs and QC figures only.")

    def _results_dir():
        os.makedirs(base_results, exist_ok=True); return base_results

    def _threshold_dir(padj: float, delta: float) -> str:
        name = f"{prefix_name}_padj{padj}_delta{delta}"
        d    = os.path.join(outdir, name)
        os.makedirs(d, exist_ok=True)
        return d

    gene_info  = load_bed(args.bed)
    mapping_df = load_and_validate_mapping(args.mapping_file)

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

    # IPA setup
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

    # Parse GTF upfront if provided
    gtf_junctions: Optional[Dict[str, Dict]] = None
    if args.gtf:
        print(f"\nParsing GTF for annotated junctions (upfront) ...")
        t_gtf = time.time()
        gene_names = [g for g, _, _, _ in gene_groups]
        gtf_junctions = parse_gtf_junctions(args.gtf, gene_names)
        n_matched = sum(1 for g in gene_names if g in gtf_junctions)
        print(f"  → matched {n_matched}/{len(gene_names)} genes ({time.time()-t_gtf:.2f}s)")
    else:
        print("[INFO] --gtf not provided; junction_type column will be omitted.")

    n_genes = len(gene_groups)
    _metrics = ["PSI_approx"]
    if not approx_only:
        _metrics += ["PSI", "full_IR_ratio"]
        if not args.no_ss_IR:
            _metrics += ["5ss_IR_ratio", "3ss_IR_ratio"]
        if genome_path:
            _metrics.append("IPA_ratio")
    print(f"\nWill process {n_genes} gene(s)")
    print(f"Metrics: {', '.join(_metrics)}")
    print(f"Threads per gene: {args.threads}\n")

    def _write_tsv(df, path): df.to_csv(path, sep="\t", index=False)

    strand_map = {g: info[2] for g, info in gene_info.items()}

    def _fdr_and_write(
        all_results: List[pd.DataFrame],
        results_subdir: str,
        p_cols: List[str],
        padj_threshold: float,
        delta_threshold: float,
    ) -> Optional[pd.DataFrame]:
        if not all_results: return None

        total_rows = sum(len(r) for r in all_results)
        _hdr = "  Correcting p-values and identifying outliers..."
        print(f"\n{'=' * 70}")
        print(_hdr)
        print(f"{'=' * 70}")
        print(f"\n  Applying FDR-BH ({total_rows:,} rows across {len(p_cols)} metric(s)) ...")
        t_fdr = time.time()

        final_df = pd.concat(all_results, ignore_index=True)
        final_df = final_df.sort_values(
            ["gene", "junction", "sample", "phasing"], ignore_index=True
        )

        # Assign junction_type if GTF was provided
        if gtf_junctions is not None:
            final_df = assign_junction_types(final_df, gtf_junctions)

        fmt_df = final_df.copy()

        for p_col in p_cols:
            padj_col = p_col.replace("p_value_", "padj_")
            if p_col not in fmt_df.columns: continue
            p_vals = pd.to_numeric(fmt_df[p_col], errors="coerce")
            is_ss  = any(x in p_col for x in ("5ss_IR_ratio", "3ss_IR_ratio"))
            if is_ss:
                ss_pos_col = "5ss" if "5ss" in p_col else "3ss"
                dedup_key = list(zip(fmt_df["sample"], fmt_df["gene"],
                                     fmt_df["phasing"], fmt_df[ss_pos_col]))
                seen: Dict[tuple, int] = {}
                dedup_idx = []
                for i, k in enumerate(dedup_key):
                    if k not in seen:
                        seen[k] = i; dedup_idx.append(i)
                dedup_p = p_vals.iloc[dedup_idx]
                valid_mask = dedup_p.notna()
                padj_dedup = np.full(len(dedup_p), np.nan)
                if valid_mask.sum() > 0:
                    _, pv, _, _ = multipletests(dedup_p[valid_mask].to_numpy(), method="fdr_bh")
                    padj_dedup[valid_mask.to_numpy()] = pv
                key_to_padj = {k: padj_dedup[j] for j, k in enumerate(
                    [dedup_key[i] for i in dedup_idx])}
                p_str = fmt_df[p_col]
                fmt_df[padj_col] = [
                    key_to_padj.get(k, p_str.iloc[i])
                    for i, k in enumerate(dedup_key)
                ]
            else:
                p_str  = fmt_df[p_col]
                valid  = p_vals.notna()
                padj_vals = p_str.copy().astype(object)
                padj_arr  = np.full(valid.sum(), np.nan)
                if valid.sum() > 0:
                    _, pv, _, _ = multipletests(p_vals[valid].to_numpy(), method="fdr_bh")
                    padj_arr = pv
                padj_vals[valid] = padj_arr
                # Rows with sentinel strings (low_n, low_coverage, etc.) should mirror p_value
                is_sentinel = ~valid & p_str.notna()
                padj_vals[is_sentinel] = p_str[is_sentinel]
                fmt_df[padj_col] = padj_vals

        fmt_df = select_output_columns(fmt_df)
        print(f"       FDR done ({time.time()-t_fdr:.2f}s)")
        return fmt_df

    results_subdir = _results_dir()
    has_ipa_metric = (not approx_only) and (genome_path is not None)
    all_results: List[pd.DataFrame] = []

    for gene, region, strand, region_df in gene_groups:
        t_gene = time.time()
        try:
            res, _ = process_one_region(
                gene=gene, region=region, region_df=region_df,
                approx_only=approx_only,
                coverage_threshold=args.coverage_threshold,
                PSI_rescale_factor=args.PSI_rescale_factor,
                n_threshold=args.n_threshold,
                phasing_threshold=args.phasing_threshold,
                threads=args.threads, strand=strand,
                include_monoexonic=args.include_monoexonic,
                min_jxn_reads=args.min_jxn_reads,
                genome_path=genome_path, alu=alu,
                no_ss_ir=args.no_ss_IR,
                loo=args.loo,
            )
        except Exception as e:
            print(f"[ERROR] Gene {gene}: {e}"); traceback.print_exc()
            res = None
        if res is not None and len(res):
            all_results.append(res)
        print(f"  {gene} complete ({time.time() - t_gene:.0f}s)")

    if not all_results:
        print("  No results."); return

    sample_df = all_results[0]
    p_cols = [c for c in sample_df.columns if c.startswith("p_value_")]

    # FDR correction (done once, shared across all thresholds)
    # Use dummy thresholds for initial n_sample_outlier if none provided
    _pt = threshold_pairs[0][0] if threshold_pairs else 0.01
    _dt = threshold_pairs[0][1] if threshold_pairs else 0.1
    final_df = _fdr_and_write(all_results, results_subdir, p_cols, _pt, _dt)
    if final_df is None: return

    # ---- Write per-gene TSVs once (threshold-independent) ----
    fmt_base = select_output_columns(final_df)
    n_genes_write = fmt_base["gene"].nunique()
    print(f"  Writing {n_genes_write} per-gene results files ...")
    t_w = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = [
            ex.submit(_write_tsv, gdf,
                      os.path.join(results_subdir, f"{g}.tsv"))
            for g, gdf in fmt_base.groupby("gene")
        ]
        for fut in concurrent.futures.as_completed(futs):
            fut.result()
    print(f"       File writes done ({time.time()-t_w:.2f}s)")

    # ---- Per-threshold outlier identification, classification, and output ----
    for padj_threshold, delta_threshold in threshold_pairs:
        t_thr = time.time()
        thr_dir = _threshold_dir(padj_threshold, delta_threshold)
        print(f"\n{'='*70}")
        print(f"  Threshold: padj <= {padj_threshold}, |delta| >= {delta_threshold}")
        print(f"  Output: {thr_dir}")
        print(f"{'='*70}")

        # ---- Outlier mask (any metric) ----
        outlier_mask = pd.Series(False, index=final_df.index)
        for p_col in p_cols:
            padj_c  = p_col.replace("p_value_", "padj_")
            delta_c = p_col.replace("p_value_", "delta_")
            if padj_c not in final_df.columns or delta_c not in final_df.columns: continue
            padj_v  = pd.to_numeric(final_df[padj_c],  errors="coerce")
            delta_v = pd.to_numeric(final_df[delta_c], errors="coerce")
            outlier_mask |= (padj_v.le(padj_threshold) & delta_v.abs().ge(delta_threshold))

        sig_df = final_df[outlier_mask].copy()

        _IR_IPA_METRICS = frozenset((
            "5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio", "junction_IPA_ratio"
        ))
        print(f"  Identifying outliers with padj <= {padj_threshold} and |delta| >= {delta_threshold} ...")
        for p_col in p_cols:
            mc     = p_col.replace("p_value_", "")
            padj_c = f"padj_{mc}"; delta_c = f"delta_{mc}"
            if padj_c not in final_df.columns or delta_c not in final_df.columns: continue
            pv = pd.to_numeric(final_df[padj_c], errors="coerce")
            dv = pd.to_numeric(final_df[delta_c], errors="coerce")
            mask_mc = pv.le(padj_threshold) & dv.abs().ge(delta_threshold)
            n_rows  = int(mask_mc.sum())
            n_jxns  = int(final_df[mask_mc]["junction"].nunique())
            if mc == "5ss_IR_ratio" and "5ss" in final_df.columns:
                n_ss = int(final_df[mask_mc]["5ss"].nunique())
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions, {n_ss} unique 5ss)")
            elif mc == "3ss_IR_ratio" and "3ss" in final_df.columns:
                n_ss = int(final_df[mask_mc]["3ss"].nunique())
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions, {n_ss} unique 3ss)")
            else:
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions)")

        # ---- Unreliable haplotype outlier filter per metric ----
        # A haplotype (sample, junction) pair is unreliable when EITHER:
        #   1. bulk is NOT sandwiched: NOT (hap1 <= bulk <= hap2 OR hap1 >= bulk >= hap2)
        #   2. one haplotype dominates asymmetrically:
        #      max(|hap1-bulk|, |hap2-bulk|) / min(|hap1-bulk|, |hap2-bulk|) > 10
        # If either condition holds, the pair is unreliable and outlier_{mc} is set to False.
        unreliable_hap_outliers: Dict[str, set] = {}

        if not sig_df.empty:
            key_col = list(zip(sig_df["sample"], sig_df["junction"]))
            sig_df["_key"] = key_col

            for p_col in p_cols:
                mc         = p_col.replace("p_value_", "")
                rescaled_c = f"rescaled_{mc}"
                if rescaled_c not in final_df.columns:
                    unreliable_hap_outliers[mc] = set()
                    continue
                phasing_df = final_df[final_df["phasing"].isin(["bulk", "hap1", "hap2"])][
                    ["sample", "junction", "phasing", rescaled_c]
                ].copy()
                phasing_df[rescaled_c] = pd.to_numeric(phasing_df[rescaled_c], errors="coerce")
                pivot = phasing_df.pivot_table(
                    index=["sample", "junction"], columns="phasing",
                    values=rescaled_c, aggfunc="first"
                )
                missing = [c for c in ("bulk", "hap1", "hap2") if c not in pivot.columns]
                if missing:
                    unreliable_hap_outliers[mc] = set()
                    continue
                pivot = pivot.dropna(subset=["bulk", "hap1", "hap2"])
                if pivot.empty:
                    unreliable_hap_outliers[mc] = set()
                    continue

                b  = pivot["bulk"]
                h1 = pivot["hap1"]
                h2 = pivot["hap2"]

                # Condition 1: bulk is sandwiched between hap1 and hap2
                sandwiched = ((h1 <= b) & (b <= h2)) | ((h2 <= b) & (b <= h1))

                # Condition 2: ratio of larger to smaller deviation <= 10
                d1   = (b - h1).abs()
                d2   = (b - h2).abs()
                dmax = np.maximum(d1, d2)
                dmin = np.minimum(d1, d2)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(dmin > 0, dmax / dmin, np.inf)
                symmetric = ratio <= 10

                unreliable_hap_outliers[mc] = set(pivot.index[~sandwiched | ~symmetric])

            sig_df = sig_df.drop(columns=["_key"])
            print(f"  Removing unreliable haplotype-only outliers ...")
            for p_col in p_cols:
                mc         = p_col.replace("p_value_", "")
                unreliable = unreliable_hap_outliers.get(mc, set())
                padj_c = f"padj_{mc}"; delta_c = f"delta_{mc}"
                if not unreliable: continue
                if padj_c not in sig_df.columns or delta_c not in sig_df.columns: continue
                hap_rows_mc = sig_df[sig_df["phasing"].isin(["hap1", "hap2"])]
                pv = pd.to_numeric(hap_rows_mc[padj_c], errors="coerce")
                dv = pd.to_numeric(hap_rows_mc[delta_c], errors="coerce")
                passes_mc = pv.le(padj_threshold) & dv.abs().ge(delta_threshold)
                keys = list(zip(
                    hap_rows_mc.loc[passes_mc, "sample"],
                    hap_rows_mc.loc[passes_mc, "junction"],
                ))
                n_tot = sum(1 for k in keys if k in unreliable)
                print(f"       {mc}: removed {n_tot:,} unreliable haplotype outlier rows")

            n_total_jxns = sig_df["junction"].nunique() if not sig_df.empty else 0
            print(f"  {len(sig_df):,} total outlier rows ({n_total_jxns} unique junctions)")
        else:
            unreliable_hap_outliers = {p.replace("p_value_", ""): set() for p in p_cols}


        # ---- outlier_{metric} boolean columns ----
        sig_df = sig_df.copy()
        for p_col in p_cols:
            mc     = p_col.replace("p_value_", "")
            padj_c = f"padj_{mc}"; delta_c = f"delta_{mc}"
            ocol   = f"outlier_{mc}"
            if padj_c not in sig_df.columns or delta_c not in sig_df.columns:
                sig_df[ocol] = False; continue
            pv = pd.to_numeric(sig_df[padj_c], errors="coerce")
            dv = pd.to_numeric(sig_df[delta_c], errors="coerce")
            passes = pv.le(padj_threshold) & dv.abs().ge(delta_threshold)
            unreliable = unreliable_hap_outliers.get(mc, set())
            if unreliable:
                is_hap       = sig_df["phasing"].isin(["hap1", "hap2"])
                keys         = list(zip(sig_df["sample"], sig_df["junction"]))
                is_unreliable = pd.Series([k in unreliable for k in keys], index=sig_df.index)
                sig_df[ocol] = passes & ~(is_hap & is_unreliable)
            else:
                sig_df[ocol] = passes

        # Per-metric breakdown of final outlier counts
        _has_jxn_type = "junction_type" in sig_df.columns
        for p_col in p_cols:
            mc      = p_col.replace("p_value_", "")
            ocol    = f"outlier_{mc}"; delta_c = f"delta_{mc}"
            if ocol not in sig_df.columns: continue
            mask_mc = sig_df[ocol].astype(bool)
            n_rows  = int(mask_mc.sum())
            n_jxns  = int(sig_df[mask_mc]["junction"].nunique())
            is_5ss  = mc == "5ss_IR_ratio" and "5ss" in sig_df.columns
            is_3ss  = mc == "3ss_IR_ratio" and "3ss" in sig_df.columns
            ss_col  = "5ss" if is_5ss else ("3ss" if is_3ss else None)
            ss_lbl  = "5ss" if is_5ss else ("3ss" if is_3ss else None)
            if mc in _IR_IPA_METRICS and delta_c in sig_df.columns:
                dv       = pd.to_numeric(sig_df[delta_c], errors="coerce")
                pos_mask = mask_mc & dv.ge(delta_threshold)
                n_pos    = int(sig_df[pos_mask]["junction"].nunique())
                ss_total = f", {int(sig_df[mask_mc][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                ss_pos   = f", {int(sig_df[pos_mask][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                if _has_jxn_type:
                    can_ann = sig_df["junction_type"].isin(["canonical", "annotated"])
                    n_can   = int(sig_df[pos_mask & can_ann]["junction"].nunique())
                    ss_can  = f", {int(sig_df[pos_mask & can_ann][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                    print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions{ss_total} → {n_pos} unique junctions{ss_pos} with delta > 0 → {n_can} unique junctions{ss_can} canonical or annotated)")
                else:
                    print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions{ss_total} → {n_pos} unique junctions{ss_pos} with delta > 0)")
            elif ss_col:
                n_ss = int(sig_df[mask_mc][ss_col].nunique())
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions, {n_ss} unique {ss_lbl})")
            else:
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions)")

                # ---- Compute n_sample_outlier for outlier TSV (threshold-specific) ----
        fmt_updated = final_df.copy()
        for p_col in p_cols:
            metric_name = p_col.replace("p_value_", "")
            n_col  = f"n_sample_outlier_{metric_name}"
            ocol   = f"outlier_{metric_name}"
            if sig_df.empty or ocol not in sig_df.columns:
                fmt_updated[n_col] = 0
                continue
            counts = (
                sig_df[sig_df[ocol].astype(bool)]
                .groupby(["gene", "junction"])["sample"]
                .nunique().rename(n_col).reset_index()
            )
            fmt_updated = fmt_updated.merge(counts, on=["gene", "junction"], how="left")
            fmt_updated[n_col] = fmt_updated[n_col].fillna(0).astype(int)

        # ---- Event classification ----
        n_sig = len(sig_df)
        print(f"  Outliers identified ({time.time() - t_thr:.2f}s)")
        print(f"  Classifying events ...")
        t_cls = time.time()
        if n_sig > 0:
            sig_df = classify_all_events(
                sig_df, final_df, strand_map,
                padj_threshold, delta_threshold, args.threads, has_ipa_metric,
                p_cols, unreliable_hap_outliers,
            )
        else:
            sig_df["event_type"] = "none"
        print(f"       done ({time.time()-t_cls:.2f}s)  "
              f"{sig_df['event_type'].value_counts().to_dict() if n_sig else {}}")

        # ---- Sort outliers ----
        if n_sig > 0:
            delta_cols = [c for c in sig_df.columns if c.startswith("delta_")]
            if delta_cols:
                delta_num = sig_df[delta_cols].apply(lambda col: pd.to_numeric(col, errors="coerce"))
                sig_df["_max_abs_delta"] = delta_num.abs().max(axis=1)

                # Per-row: does it have a named (non-'other') event, and what's the
                # max abs delta of the metric(s) that produced it?
                all_metric_event_sets = {
                    mc: frozenset(evs) for mc, evs in _METRIC_EVENTS.items()
                }
                def _named_event_delta(row):
                    et  = row.get("event_type", "other") or "other"
                    evs = set(e.strip() for e in et.split(",")) - {"other", ""}
                    if not evs:
                        return 0, 0.0
                    best = 0.0
                    for mc, mc_evs in all_metric_event_sets.items():
                        if evs & mc_evs:
                            dc = f"delta_{mc}"
                            if dc in row:
                                v = pd.to_numeric(row[dc], errors="coerce")
                                if pd.notna(v):
                                    best = max(best, abs(v))
                    return 1, best

                if "event_type" in sig_df.columns:
                    named_info = sig_df.apply(_named_event_delta, axis=1, result_type="expand")
                    sig_df["_has_named_event"]   = named_info[0]
                    sig_df["_named_event_delta"]  = named_info[1]
                else:
                    sig_df["_has_named_event"]   = 0
                    sig_df["_named_event_delta"]  = 0.0

                gene_agg = (
                    sig_df.groupby(["sample", "gene"])
                    .agg(
                        _gene_has_named  =("_has_named_event",   "max"),
                        _gene_named_delta=("_named_event_delta", "max"),
                        _gene_max_delta  =("_max_abs_delta",     "max"),
                    )
                    .reset_index()
                )
                gene_agg["_sort_key"] = list(zip(
                    -gene_agg["_gene_has_named"],
                    -gene_agg["_gene_named_delta"],
                    -gene_agg["_gene_max_delta"],
                ))
                gene_agg["_gene_rank"] = (
                    gene_agg.groupby("sample")["_sort_key"]
                    .rank(ascending=True, method="min")
                )
                sig_df = sig_df.merge(
                    gene_agg[["sample", "gene", "_gene_rank"]],
                    on=["sample", "gene"], how="left"
                )
                sig_df = sig_df.sort_values(
                    ["sample", "_gene_rank", "junction"],
                    ascending=[True, True, True],
                ).drop(columns=["_max_abs_delta", "_has_named_event", "_named_event_delta"])
                sig_df = sig_df.rename(columns={"_gene_rank": "gene_rank"})

        # ---- Build outlier TSV column order ----
        # _OUTPUT_COLS + n_sample_outlier_{metric} + outlier_{metric} after each padj_{metric}
        # + event_type at end
        # Merge n_sample_outlier columns from fmt_updated into sig_df
        n_sample_cols = [f"n_sample_outlier_{p.replace('p_value_', '')}" for p in p_cols]
        n_sample_cols_present = [c for c in n_sample_cols if c in fmt_updated.columns]
        if n_sample_cols_present:
            sig_df = sig_df.merge(
                fmt_updated[["sample", "gene", "junction", "phasing"] + n_sample_cols_present]
                .drop_duplicates(subset=["sample", "gene", "junction", "phasing"]),
                on=["sample", "gene", "junction", "phasing"], how="left"
            )

        computed_metrics_ordered = [p.replace("p_value_", "") for p in p_cols]
        outlier_tsv_cols = []
        for col in _OUTPUT_COLS:
            outlier_tsv_cols.append(col)
            for mc in computed_metrics_ordered:
                if col == f"padj_{mc}":
                    n_sc = f"n_sample_outlier_{mc}"
                    if n_sc in sig_df.columns:
                        outlier_tsv_cols.append(n_sc)
                    outlier_tsv_cols.append(f"outlier_{mc}")
        outlier_tsv_cols.append("event_type")
        outlier_tsv_cols = [c for c in outlier_tsv_cols if c in sig_df.columns]

        # ---- Parallel output: TSVs + heatmaps + box plots ----
        t_out = time.time()
        out_jobs = []

        # TSV jobs — need the full prepared slice
        _outliers_data = sig_df[outlier_tsv_cols].copy()
        if "event_type" in sig_df.columns:
            filt = sig_df[sig_df["event_type"] != "none"]
            _outliers_filt_data = filt[outlier_tsv_cols].copy() if not filt.empty else pd.DataFrame(columns=outlier_tsv_cols)
        else:
            _outliers_filt_data = pd.DataFrame(columns=outlier_tsv_cols)

        def _write_outliers(data=_outliers_data):
            t0 = time.time()
            out = os.path.join(thr_dir, f"{prefix_name}_outliers.tsv")
            data.to_csv(out, sep="\t", index=False)
            print(f"  Outliers → {out} ({time.time()-t0:.2f}s)")

        def _write_outliers_filtered(data=_outliers_filt_data):
            t0 = time.time()
            out_filt = os.path.join(thr_dir, f"{prefix_name}_outliers_filtered.tsv")
            data.to_csv(out_filt, sep="\t", index=False)
            print(f"  Outliers (filtered) → {out_filt} ({time.time()-t0:.2f}s)")


        # ---- Build outlier_map for box plots ----
        # Built from filt (outliers_filtered rows: event_type != "none", outlier_{mc} True),
        # so boxplots exactly match outliers_filtered.tsv.
        outlier_map: Dict[str, Dict[str, Dict[str, set]]] = {}
        for p_col in p_cols:
            mc   = p_col.replace("p_value_", "")
            ocol = f"outlier_{mc}"
            gene_jxn_map: Dict[str, Dict[str, set]] = defaultdict(lambda: defaultdict(set))
            if not filt.empty and ocol in filt.columns:
                passing = filt[filt[ocol].astype(bool)]
                for _, row in passing.iterrows():
                    gene_jxn_map[row["gene"]][row["junction"]].add(row["sample"])
            outlier_map[mc] = gene_jxn_map

        out_jobs.append(_write_outliers)
        out_jobs.append(_write_outliers_filtered)

        # Heatmap + box plot jobs per metric — pre-slice to only needed columns/rows
        if n_sig > 0:
            for mc in [p.replace("p_value_", "") for p in p_cols]:
                dc = f"delta_{mc}"; pc = f"padj_{mc}"
                if dc not in sig_df.columns or pc not in sig_df.columns:
                    continue

                # Heatmap: only gene, sample, padj, delta — rows passing threshold
                pv = pd.to_numeric(sig_df[pc], errors="coerce")
                dv = pd.to_numeric(sig_df[dc], errors="coerce")
                heat_mask = pv.le(padj_threshold) & dv.abs().ge(delta_threshold)
                heat_df = sig_df.loc[heat_mask, ["gene", "sample", dc, pc]].copy()

                def _make_heatmap(df=heat_df, mc=mc, dc=dc, pc=pc):
                    make_outlier_heatmap(
                        df, mc, dc, pc,
                        padj_threshold, delta_threshold,
                        os.path.join(thr_dir, f"{prefix_name}_outlier_heatmap_{mc}.pdf"),
                    )
                out_jobs.append(_make_heatmap)

                rc = f"rescaled_{mc}"
                if rc in fmt_updated.columns:
                    # Box plots: only gene, junction, phasing, sample, rescaled col
                    # filtered to genes/junctions in outlier_map[mc]
                    gj_map = outlier_map.get(mc, {})
                    if gj_map:
                        bp_genes = set(gj_map.keys())
                        bp_jxns  = set(j for jxns in gj_map.values() for j in jxns)
                        bp_df = fmt_updated.loc[
                            fmt_updated["gene"].isin(bp_genes) &
                            fmt_updated["junction"].isin(bp_jxns),
                            ["gene", "junction", "phasing", "sample", rc]
                        ].copy()
                        def _make_boxplot(df=bp_df, mc=mc, rc=rc):
                            make_hit_boxplots(
                                df, outlier_map, mc, rc, thr_dir, prefix_name, tmp_dir,
                                padj_threshold=padj_threshold,
                                delta_threshold=delta_threshold,
                                is_ir_ipa=(mc in _IR_IPA_METRICS),
                                has_jxn_type_filter=("junction_type" in fmt_updated.columns
                                                     and mc in _IR_IPA_METRICS),
                            )
                        out_jobs.append(_make_boxplot)

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = [ex.submit(fn) for fn in out_jobs]
            for fut in concurrent.futures.as_completed(futs):
                try:
                    fut.result()
                except Exception as e:
                    print(f"[WARNING] Output job failed: {e}"); traceback.print_exc()

        print(f"  All outputs written ({time.time()-t_out:.2f}s)")

        if os.path.exists(tmp_dir):
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)

        # ---- Outlier summary for this threshold ----
        try:
            _METRIC_DISPLAY = {
                "junction_PSI":            "PSI",
                "junction_PSI_approx":     "PSI_approx",
                "5ss_IR_ratio":            "5ss_IR",
                "3ss_IR_ratio":            "3ss_IR",
                "junction_full_IR_ratio":  "full_IR",
                "junction_IPA_ratio":      "IPA",
            }
            _METRIC_ORDER   = list(_METRIC_DISPLAY.keys())
            computed_mets   = {c.replace("p_value_", "") for c in p_cols}
            metrics_present = [m for m in _METRIC_ORDER if m in computed_mets]

            if "sample" in sig_df.columns and "gene" in sig_df.columns:
                rows_any = sig_df[["sample", "gene"]].drop_duplicates()
                any_ser  = rows_any.groupby("sample")["gene"].nunique().rename("_any")
                metric_sers = []
                for mc in metrics_present:
                    ocol = f"outlier_{mc}"
                    if ocol not in sig_df.columns: continue
                    ms = (sig_df[sig_df[ocol].astype(bool)][["sample", "gene"]]
                          .drop_duplicates()
                          .groupby("sample")["gene"].nunique().rename(mc))
                    metric_sers.append(ms)
                summary = pd.concat([any_ser] + metric_sers, axis=1)
                summary = summary.fillna(0)
                for col in summary.columns:
                    summary[col] = pd.to_numeric(summary[col], errors="coerce").fillna(0).astype(int)
                summary = summary.sort_values("_any", ascending=False).reset_index()
                header_metrics = ["Genes"] + [_METRIC_DISPLAY[m] for m in metrics_present
                                               if m in summary.columns]
                col_w  = max(10, max(len(h) for h in header_metrics) + 2)
                header = f"  {'Sample':<40}" + "".join(f"{h:>{col_w}}" for h in header_metrics)
                sep    = f"  {'-'*40}" + ("-"*col_w) * len(header_metrics)
                print("\n" + "="*len(sep.rstrip()))
                print(f"  Outlier summary — padj<={padj_threshold}, |delta|>={delta_threshold}")
                print(f"  (samples ranked by total genes with outlier)")
                print("="*len(sep.rstrip()))
                print(header); print(sep)
                for _, row in summary.iterrows():
                    vals = f"  {row['sample']:<40}{row['_any']:>{col_w}}"
                    for m in metrics_present:
                        if m in summary.columns:
                            vals += f"{row.get(m, 0):>{col_w}}"
                    print(vals)
                print("="*len(sep.rstrip()))
        except Exception as e:
            print(f"[WARNING] Outlier summary failed: {e}"); traceback.print_exc()

        print(f"\n  Finished threshold: padj <= {padj_threshold}, |delta| >= {delta_threshold} ({time.time()-t_thr:.2f}s)")

    # end threshold loop

    # ---- QC figures (uses already-parsed gtf_junctions) ----
    if args.gtf and gtf_junctions is not None:
        gene_names = [g for g, _, _, _ in gene_groups]
        print(f"\n{'='*56}")
        print(f"  QC")
        print(f"{'='*56}")
        t_qc = time.time()
        print(f"\n  {'Gene':<20} {'Canonical jxns':>16} {'Annotated jxns':>16}")
        print(f"  {'-'*20} {'-'*16} {'-'*16}")
        for g in gene_names:
            if g in gtf_junctions:
                n_can = len(gtf_junctions[g]["canonical_junctions"])
                n_ann = len(gtf_junctions[g]["all_junctions"])
                print(f"  {g:<20} {n_can:>16} {n_ann:>16}")
            else:
                print(f"  {g:<20} {'NOT FOUND':>16} {'NOT FOUND':>16}")

        results_subdir = _results_dir()

        # Use final_df directly — bulk rows only, pre-sliced per figure
        has_ipa_plot = (not approx_only) and (genome_path is not None)
        os.makedirs(qc_dir, exist_ok=True)
        bulk_df = final_df[final_df["phasing"] == "bulk"]

        qc_jobs = []
        for cov_col, file_suffix, bar_metrics in _QC_FIGURES:
            if approx_only and cov_col != "junction_coverage_approx":
                continue
            if args.no_ss_IR and cov_col in ("5ss_coverage", "3ss_coverage"):
                continue
            companions = [m for m in bar_metrics
                          if has_ipa_plot or m != "junction_IPA_ratio"]
            if not companions:
                continue
            needed_cols = ["gene", "sample", "junction", cov_col]
            for cm in companions:
                ac = f"alpha_{cm}"
                if ac in bulk_df.columns:
                    needed_cols.append(ac)
            needed_cols = list(dict.fromkeys(c for c in needed_cols if c in bulk_df.columns))
            qc_slice = bulk_df[needed_cols].copy()

            for subset in ("canonical", "annotated"):
                out_pdf = os.path.join(qc_dir, f"{prefix_name}_qc_{file_suffix}_{subset}.pdf")
                def _qc_job(df=qc_slice, cc=cov_col, cm=companions, fs=file_suffix,
                            ss=subset, op=out_pdf):
                    make_qc_figure(
                        df, gtf_junctions,
                        cc, cm, fs, ss, args.coverage_threshold, op,
                        has_ipa_plot,
                    )
                qc_jobs.append(_qc_job)

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
                futs = [ex.submit(fn) for fn in qc_jobs]
                for fut in concurrent.futures.as_completed(futs):
                    try:
                        fut.result()
                    except Exception as e:
                        print(f"[WARNING] QC figure failed: {e}"); traceback.print_exc()

        print(f"\n  Finished QC ({time.time()-t_qc:.2f}s)")

    print("\nDone.")


if __name__ == "__main__":
    main()