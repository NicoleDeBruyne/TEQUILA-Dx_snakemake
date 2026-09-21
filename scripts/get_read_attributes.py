
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
    merged, starts = load_on_target_intervals(bed) if bed else (None, None)

    lengths_by_type = defaultdict(list)
    n_seen = 0
    heartbeat_every = 2_000_000
    with pysam.AlignmentFile(bam, "rb", threads=2) as bamfile:
        for read in bamfile.fetch(until_eof=True):
            n_seen += 1
            if n_seen % heartbeat_every == 0:
                print(f"  ...{sample}: {n_seen:,} alignments read so far", flush=True)

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
