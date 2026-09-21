
import argparse


def add_alias_map_arg(parser: argparse.ArgumentParser, required: bool = False):
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
    return alias_map.get(sample, sample)


def resolve_all(samples, alias_map):
    return [resolve(s, alias_map) for s in samples]
