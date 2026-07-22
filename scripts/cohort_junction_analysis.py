#!/usr/bin/env python3
"""
scripts/cohort_junction_analysis.py

Core per-gene splice-junction metric computation for cohort-level outlier
analysis. Companion script scripts/identify_cohort_junction_outliers.py
performs the downstream statistical testing/filtering stage.

For each gene in the input BED file, reads every cohort sample's bulk/hap1/
hap2 BAMs (from --mapping-file) and computes per-junction, per-sample,
per-phasing coverage and ratio metrics (junction usage/PSI, 5'/3' intron-
retention ratios, intronic-polyadenylation ratio) -- no statistical testing,
no cross-sample comparison.

Writes one raw metrics TSV per gene into --outdir, plus a manifest TSV
(--manifest) with one row per gene in the input BED file: the gene name and
the path to that gene's result file, or the literal string "None" if the
gene produced no output (not present in --mapping-file, no BAMs found, an
error during processing, etc.). The manifest's presence marks this rule as
complete for Snakemake's purposes. scripts/identify_cohort_junction_outliers.py
reads the manifest to find each gene's raw results for statistical testing.
"""

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
from pandas.errors import PerformanceWarning

warnings.filterwarnings("ignore", category=PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


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


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

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
    p.add_argument("--coverage-threshold",         type=int,   default=20)
    p.add_argument("--PSI-rescale-factor",         type=float, default=1e-3)
    p.add_argument("--phasing-threshold",          type=float, default=0.8)
    p.add_argument("--min-jxn-reads",              type=int,   default=20)
    p.add_argument("--include-monoexonic",         action="store_true")
    p.add_argument("--genome",                     default=None)
    p.add_argument("--alu-bed",                    default=None)
    p.add_argument("--threads",                    type=int,   default=1)
    p.add_argument("--test-n-genes",               type=int,   default=None)
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
# Per-gene junction discovery + coverage/usage metric computation
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


def compute_gene_junction_metrics(
    gene: str,
    region: str,
    region_df: pd.DataFrame,
    approx_only: bool,
    coverage_threshold: int,
    PSI_rescale_factor: float,
    phasing_threshold: float,
    threads: int,
    strand: str = "+",
    include_monoexonic: bool = False,
    min_jxn_reads: int = 20,
    genome_path: Optional[str] = None,
    alu: Optional[dict] = None,
) -> Optional[pd.DataFrame]:
    """
    Discovers this gene's splice junctions from every cohort sample's bulk
    BAM, then computes per-junction, per-sample, per-phasing coverage/usage
    metrics (junction usage/PSI, 5'/3' intron-retention ratios, intronic-
    polyadenylation ratio). No statistical testing is performed here --
    see scripts/identify_cohort_junction_outliers.py for that.

    Returns None if no data could be collected for this gene.
    """
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
                    print(f"  [WARNING] Skipping {gene}."); return None
            if not sample_dfs:
                print("  No data. Skipping."); return None
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
            print(f"  [1/2] Junction discovery (bulk BAMs, {len(bulk_rows)} samples) ...")
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
                print("  No junctions passed filter. Skipping."); return None

            print(f"  [2/2] Calculating junction usage ({len(region_df)} samples) ...")
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
                    include_monoexonic, genome_path, alu,
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
                    print(f"  [WARNING] Skipping {gene}."); return None
            if not sample_dfs:
                print("  No data. Skipping."); return None
            combined = pd.concat(sorted(sample_dfs, key=lambda d: d["sample"].iloc[0]), ignore_index=True)
            _psi_wall = time.time() - _psi_t0
            for _line in sorted(_psi_lines): print(_line)
            print(f"       → {len(combined):,} rows (junctions × samples × phasings) "
                  f"({_psi_wall:.2f}s)")

    combined = combined.sort_values(
        ["junction", "sample", "phasing"], ignore_index=True
    )
    return combined


# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

def _write_manifest(gene_info: Dict[str, Tuple[str, str, str]],
                     result_paths: Dict[str, Optional[str]],
                     manifest_path: str) -> None:
    """Writes one row per gene in gene_info (BED order): gene, result_path
    (or the literal string 'None')."""
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

    # ---- Minimum cohort size gate ----
    # Fitting a per-junction distribution across the cohort's bulk samples
    # isn't meaningful with only a handful of samples, so small groups are
    # skipped entirely rather than producing unreliable results.
    n_samples = mapping_df["sample"].nunique()
    if n_samples < args.min_samples:
        msg = (f"SKIPPED: only {n_samples} sample(s) in this group "
               f"(minimum {args.min_samples} required to run cohort junction analysis). "
               f"No genes were processed.")
        print(f"\n[WARNING] {msg}")
        with open(args.note, "w") as fh:
            fh.write(msg + "\n")
        _write_manifest(gene_info, {}, args.manifest)
        print("\nDone (skipped).")
        return

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

    n_genes = len(gene_groups)
    print(f"\nWill process {n_genes} gene(s)")
    print(f"Threads per gene: {args.threads}\n")

    # ---- Process each gene, writing its raw metrics TSV as we go ----
    result_paths: Dict[str, Optional[str]] = {}

    for gene, region, strand, region_df in gene_groups:
        t_gene = time.time()
        try:
            combined = compute_gene_junction_metrics(
                gene=gene, region=region, region_df=region_df,
                approx_only=approx_only,
                coverage_threshold=args.coverage_threshold,
                PSI_rescale_factor=args.PSI_rescale_factor,
                phasing_threshold=args.phasing_threshold,
                threads=args.threads, strand=strand,
                include_monoexonic=args.include_monoexonic,
                min_jxn_reads=args.min_jxn_reads,
                genome_path=genome_path, alu=alu,
            )
        except Exception as e:
            print(f"[ERROR] Gene {gene}: {e}"); traceback.print_exc()
            combined = None

        if combined is not None and len(combined):
            out_path = os.path.join(args.outdir, f"{gene}.tsv")
            combined.to_csv(out_path, sep="\t", index=False)
            result_paths[gene] = out_path
            print(f"  {gene} → {out_path}")
        else:
            result_paths[gene] = None
            print(f"  {gene}: no output")
        print(f"  {gene} complete ({time.time() - t_gene:.0f}s)")

    # ---- Write manifest: one row per gene in the BED file, in BED order ----
    _write_manifest(gene_info, result_paths, args.manifest)
    n_with_data = sum(1 for p in result_paths.values() if p)
    with open(args.note, "w") as fh:
        fh.write(f"Ran normally: {n_samples} sample(s), {len(gene_groups)} gene(s) "
                 f"attempted, {n_with_data} gene(s) produced results.\n")
    print("\nDone.")


if __name__ == "__main__":
    main()
