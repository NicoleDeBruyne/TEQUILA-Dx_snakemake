#!/usr/bin/env python3
"""
scripts/get_on_target_rate.py
Computes one sample's total/mapped/on-target read counts from its own
BAM+BED pair. Invoked per-sample by rules/6_sample_qc.smk (_6A). Mapping
rate and on-target rate are NOT computed here -- they're cheap ratios that
the cohort-level merge step (scripts/plot_on_target_rates.py, rules/
9_merge_results.smk's _9G) derives from these raw counts, so a re-run
triggered only by cohort membership changing (not this sample's own BAM)
never needs to touch the BAM again.

Adapted from the BAM-scanning logic in the old (pre-split) version of
plot_on_target_rates.py, pulled from Github 2026.04.08.
"""

import argparse

import pandas as pd
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description="Get one sample's total/mapped/on-target read counts")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--outfile", required=True)
    return parser.parse_args()


def count_reads(bam, bed):
    """Distinct query names of primary, mapped alignments overlapping any BED interval."""
    ids = set()
    with pysam.AlignmentFile(bam, "rb") as f:
        with open(bed) as b:
            for line in b:
                if line.strip() and not line.startswith('#'):
                    c, s, e = line.split()[:3]
                    for r in f.fetch(c, int(s), int(e)):
                        if not r.is_unmapped and not r.is_secondary:
                            ids.add(r.query_name)
    return len(ids)


def count_mapped(bam):
    """Distinct query names, split into mapped vs. unmapped, across the whole BAM."""
    m, u = set(), set()
    with pysam.AlignmentFile(bam, "rb") as f:
        for r in f.fetch(until_eof=True):
            (u if r.is_unmapped else m).add(r.query_name)
    return len(m), len(u)


def main():
    args = parse_args()

    mapped, unmapped = count_mapped(args.bam)
    total = mapped + unmapped
    target = count_reads(args.bam, args.bed)

    df = pd.DataFrame([{"sample": args.sample, "total": total, "mapped": mapped, "target": target}])
    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"{args.sample}: total={total}, mapped={mapped}, target={target}")
    print(f"Saved: {args.outfile}")


if __name__ == "__main__":
    main()
