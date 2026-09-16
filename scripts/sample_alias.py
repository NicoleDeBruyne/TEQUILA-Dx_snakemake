"""Shared helper for the optional per-sample 'alias' feature.

Any script that produces an output containing real sample IDs (a TSV column/
index, or a plot label) can accept an --alias-map argument and use this
module to resolve each sample to its alias, falling back to the sample's own
ID when it has no alias configured. This lets the same computed data be
written out twice -- once under real sample IDs, once under aliases -- in a
single script invocation, without recomputing anything upstream (e.g. BAM
parsing).

CLI convention: --alias-map is a list of "sample=alias" tokens, e.g.:
    --alias-map SAMPLE_A=PT01 SAMPLE_B=PT02
Samples not listed simply map to themselves. An empty/omitted --alias-map
means "no aliases in this cohort" -- callers should treat that as "produce
an alias file that is identical to the real-ID file" rather than skip it,
so the file always exists for any output path Snakemake expects it at.
"""

import argparse


def add_alias_map_arg(parser: argparse.ArgumentParser, required: bool = False):
    """Add the standard --alias-map argument to an argparse parser."""
    parser.add_argument(
        "--alias-map",
        nargs="*",
        default=[],
        required=required,
        metavar="SAMPLE=ALIAS",
        help="Zero or more 'sample=alias' pairs. Samples not listed here "
             "fall back to their own sample ID. Used to additionally emit "
             "an alias-labeled version of this script's output(s).",
    )


def parse_alias_map(pairs):
    """Turn a list of 'sample=alias' strings (as produced by add_alias_map_arg)
    into a {sample: alias} dict. Raises ValueError on a malformed token."""
    alias_map = {}
    for pair in pairs or []:
        if "=" not in pair:
            raise ValueError(
                f"--alias-map entries must look like 'sample=alias', got: {pair!r}"
            )
        sample, alias = pair.split("=", 1)
        sample = sample.strip()
        alias = alias.strip()
        if not sample or not alias:
            raise ValueError(
                f"--alias-map entries must have a non-empty sample and alias, got: {pair!r}"
            )
        alias_map[sample] = alias
    return alias_map


def resolve(sample, alias_map):
    """Alias for one sample, falling back to the sample's own ID."""
    return alias_map.get(sample, sample)


def resolve_all(samples, alias_map):
    """Alias for each sample in an iterable, preserving order."""
    return [resolve(s, alias_map) for s in samples]
