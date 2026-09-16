#!/usr/bin/env python3
"""
scripts/simplify_all_hits.py
Reads merged_all_hits.tsv and writes a companion merged_all_hits_simplified.tsv,
keeping only the core per-sample/per-gene metadata and ranking columns plus
variant_ID, and collapsing every *_jxns column (bulk_jxns, hap1_jxns,
hap2_jxns, cohort_bulk_jxns, cohort_hap1_jxns, cohort_hap2_jxns -- matched
by suffix rather than hardcoded, so this doesn't go stale if a new phasing
tier or cohort-comparison metric is added later) into a single deduplicated
"jxns" column. For a quick read, the goal is usually just "which junctions
were flagged for this gene/sample at all", not which specific tier/phasing/
GTEx-tissue-comparison flagged each one separately -- that detail is still
available in merged_all_hits.tsv itself.

Invoked by rules/6_merge_hits.smk, as part of _6F_final_merge (right after it writes merged_all_hits.tsv).
"""

import argparse
import re

import pandas as pd

from sample_alias import add_alias_map_arg, parse_alias_map, resolve

_KEPT_COLUMNS = [
    'sample', 'gene', 'phenotypes', 'inheritance_patterns', 'haploinsufficient',
    'ranking', 'tier', 'variant', 'pathogenic_variant', 'ASE', 'outlier_junction',
    'cohort_outlier_junction', 'relative_gene_expression', 'cohort_relative_gene_expression',
    'n_cohort', 'variant_ID',
]

# '.' is this pipeline's standard missing-value sentinel (see e.g.
# merge_hits.py's cohort_*/phenotypes fallback-filling); the others are
# just defensive against however pandas/upstream tools might stringify a
# genuinely empty cell.
_MISSING_SENTINELS = {'', '.', 'nan', 'none'}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write a simplified companion to merged_all_hits.tsv, with a single deduplicated jxns column.")
    parser.add_argument('--infile', required=True, help="Path to merged_all_hits.tsv")
    parser.add_argument('--outfile', required=True, help="Path to write merged_all_hits_simplified.tsv")
    parser.add_argument('--alias-outfile',
        help="If given, also write an alias-labeled copy here, with the 'sample' "
             "column resolved through --alias-map (falling back to the real sample "
             "ID for any sample with no alias configured).")
    add_alias_map_arg(parser)
    return parser.parse_args()


def _split_jxns(cell):
    """One cell's junction list -> list of individual junction strings.
    bulk-tier *_jxns columns join distinct junctions with ';', hap1/hap2-
    tier columns join with ',' (see merge_hits.py's build_phased_junction_df)
    -- splitting on both covers either convention regardless of which
    column this particular cell came from."""
    if pd.isna(cell):
        return []
    s = str(cell).strip()
    if s.lower() in _MISSING_SENTINELS:
        return []
    return [j.strip() for j in re.split('[,;]', s)
            if j.strip() and j.strip().lower() not in _MISSING_SENTINELS]


def main():
    args = parse_args()
    df = pd.read_csv(args.infile, sep='\t', dtype=str)

    jxn_cols = [c for c in df.columns if c.endswith('_jxns')]
    if not jxn_cols:
        raise ValueError(f"No *_jxns columns found in {args.infile} -- can't build a merged jxns column.")
    print(f"Merging and deduplicating {len(jxn_cols)} junction column(s): {', '.join(jxn_cols)}")

    def merge_row_jxns(row):
        all_jxns = []
        for col in jxn_cols:
            all_jxns.extend(_split_jxns(row[col]))
        # dict.fromkeys preserves first-seen order while deduplicating --
        # same convention used elsewhere in this pipeline (e.g.
        # merge_hits.py's gnomAD_AF/CLNSIG aggregation) rather than
        # sorting, so junctions still read in roughly bulk-then-haplotype-
        # then-cohort order instead of an arbitrary alphabetical one.
        deduped = list(dict.fromkeys(all_jxns))
        return ';'.join(deduped) if deduped else '.'

    df['jxns'] = df.apply(merge_row_jxns, axis=1)

    missing = [c for c in _KEPT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{args.infile} is missing expected column(s): {missing}")

    out_df = df[_KEPT_COLUMNS + ['jxns']]
    out_df.to_csv(args.outfile, sep='\t', index=False)
    print(f"Saved simplified hits table ({len(out_df)} rows) to {args.outfile}")

    if args.alias_outfile:
        alias_map = parse_alias_map(args.alias_map)
        alias_df = out_df.copy()
        alias_df['sample'] = alias_df['sample'].apply(lambda s: resolve(s, alias_map))
        alias_df.to_csv(args.alias_outfile, sep='\t', index=False)
        print(f"Saved alias-labeled simplified hits table ({len(alias_df)} rows) to {args.alias_outfile}")


if __name__ == "__main__":
    main()
