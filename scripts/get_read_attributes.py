#!/usr/bin/env python3
"""
scripts/get_read_attributes.py
Computes one sample's read-length five-number summary (min, Q1, median, Q3,
max) per target type (on_target/off_target/mapped/unmapped), directly from
its own BAM (+ optional BED for on/off-target classification). Invoked
per-sample by rules/6_sample_qc.smk (_6B).

Read length is taken from each *primary* alignment's SEQ field
(secondary/supplementary records are skipped, since these can be
hard-clipped and would understate the true read length).

Only the five-number summary is kept, not every individual read length --
that's all the cohort-level boxplot (scripts/plot_read_attributes.py) needs,
and it turns a many-GB per-sample read-length table into a handful of
numbers that can cross the pipeline's rule/log boundaries cheaply.

Adapted from the BAM-scanning logic in the old (pre-split) version of
plot_read_attributes.py, pulled from Github 2026.04.08.
"""

import argparse
import bisect
from collections import defaultdict

import pandas as pd
import numpy as np
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description="Get one sample's read-length five-number summary per target type")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bed", default=None,
        help="Optional; on/off-target classification (vs. just mapped/unmapped) is skipped without it")
    parser.add_argument("--outfile", required=True)
    return parser.parse_args()


def load_on_target_intervals(bed):
    """{chrom: [(start, end), ...]} of merged, sorted on-target intervals
    from a BED file, plus a matching {chrom: [start, ...]} for bisecting."""
    raw = defaultdict(list)
    with open(bed) as bedfile:
        for line in bedfile:
            if line.strip() and not line.startswith('#'):
                chrom, start, end = line.split()[:3]
                raw[chrom].append((int(start), int(end)))

    merged = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        m = []
        for s, e in intervals:
            if m and s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        merged[chrom] = m

    starts = {chrom: [s for s, _ in ivs] for chrom, ivs in merged.items()}
    return merged, starts


def is_on_target(chrom, start, end, merged, starts):
    """True if [start, end) overlaps any merged on-target interval on chrom.
    merged/starts intervals are sorted and non-overlapping (see
    load_on_target_intervals), so it's enough to check the interval whose
    start is <= this read's start (the closest candidate from the left) and
    the very next one (the closest candidate from the right): if any
    interval overlapped but wasn't one of those two, it would have to sit
    strictly between them, which is impossible once they're merged."""
    ivs = merged.get(chrom)
    if not ivs:
        return False
    idx = bisect.bisect_right(starts[chrom], start) - 1
    for j in (idx, idx + 1):
        if 0 <= j < len(ivs):
            s, e = ivs[j]
            if s < end and e > start:
                return True
    return False


def get_read_lengths(sample, bam, bed):
    """Single pass over the BAM; returns {target_type: numpy int32 array of
    read lengths}. on-target/off-target status is decided inline per read
    via a binary-search interval lookup against the BED, instead of
    pre-fetching on-target read IDs region-by-region and then re-reading the
    whole file a second time."""
    merged, starts = load_on_target_intervals(bed) if bed else (None, None)

    lengths_by_type = defaultdict(list)
    n_seen = 0
    heartbeat_every = 2_000_000
    with pysam.AlignmentFile(bam, "rb", threads=2) as bamfile:
        for read in bamfile.fetch(until_eof=True):
            n_seen += 1
            if n_seen % heartbeat_every == 0:
                print(f"  ...{sample}: {n_seen:,} alignments read so far", flush=True)

            # Skip secondary/supplementary records: they can be hard-clipped,
            # which would understate the read's true length, and would
            # double-count the same underlying read alongside its primary.
            if read.is_secondary or read.is_supplementary:
                continue

            seq_len = read.query_length or 0

            if read.is_unmapped:
                target_type = "unmapped"
            elif bed:
                on_target = is_on_target(read.reference_name, read.reference_start, read.reference_end,
                                          merged, starts)
                target_type = "on_target" if on_target else "off_target"
            else:
                target_type = "mapped"

            lengths_by_type[target_type].append(seq_len)

    print(f"  {sample}: {n_seen:,} total alignments read.", flush=True)
    return {tt: np.asarray(lengths, dtype=np.int32) for tt, lengths in lengths_by_type.items()}


def main():
    args = parse_args()

    lengths_by_type = get_read_lengths(args.sample, args.bam, args.bed)

    rows = []
    for target_type, arr in lengths_by_type.items():
        if len(arr) == 0:
            continue
        q1, median, q3 = np.percentile(arr, [25, 50, 75])
        rows.append(dict(
            sample=args.sample, targetType=target_type, n=len(arr),
            min=int(arr.min()), q1=q1, median=median, q3=q3, max=int(arr.max()),
        ))

    df = pd.DataFrame(rows, columns=["sample", "targetType", "n", "min", "q1", "median", "q3", "max"])
    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"Saved: {args.outfile}")


if __name__ == "__main__":
    main()
