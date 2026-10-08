
from __future__ import annotations

import os
import sys
import argparse
import warnings
import traceback
import time
import math
import glob as _glob
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from sample_alias import add_alias_map_arg, parse_alias_map, resolve

import numpy as np
import pandas as pd
import concurrent.futures
from pandas.errors import PerformanceWarning

warnings.filterwarnings("ignore", category=PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*Glyph.*missing from font.*", category=UserWarning)
warnings.filterwarnings("ignore", message=".*Adding colorbar to a different Figure.*", category=UserWarning)


_METRIC_EVENTS: Dict[str, List[str]] = {
    "junction_PSI_approx":    ["alt_5ss_approx", "alt_3ss_approx", "exon_skipping_approx", "exon_inclusion_approx", "complex_approx"],
    "junction_PSI":           ["alt_5ss", "alt_3ss", "exon_skipping", "exon_inclusion", "complex"],
    "5ss_IR_ratio":           ["5ss_IR"],
    "3ss_IR_ratio":           ["3ss_IR"],
    "junction_full_IR_ratio": ["full_IR"],
    "junction_IPA_ratio":     ["IPA"],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Reads a combined scored TSV from fit_and_score_cohort_junctions.py, "
                     "applies thresholds, and writes outlier outputs (step 8C)."
    )
    p.add_argument("--outprefix",                  required=True)
    p.add_argument("--bed",                        required=True)
    add_alias_map_arg(p)
    p.add_argument("--approx",                     action="store_true")
    p.add_argument("--has-ipa",                    action="store_true",
                   help="Set if --genome was provided to cohort_junction_analysis.py "
                        "(enables testing junction_IPA_ratio as a metric).")
    p.add_argument("--coverage-threshold",         type=int,   default=20)
    p.add_argument("--phasing-threshold",          type=float, default=0.5,
                   help="Minimum fraction of bulk coverage that must be phased "
                        "(hap1+hap2 denom / bulk denom) for haplotype metrics to be scored.")
    method_group = p.add_mutually_exclusive_group(required=True)
    method_group.add_argument(
        "--bb-thresholds",
        nargs="*", default=None,
        metavar="PADJ:DELTA",
        help="Use beta-binomial testing thresholds. One or more padj:delta "
             "threshold pairs, e.g. 0.05:0.1 0.01:0.1 0.01:0.2 -- each "
             "combination produces its own output subdirectory.")
    method_group.add_argument(
        "--z-thresholds",
        nargs="*", default=None,
        metavar="MODZ:DELTA",
        help="Use modified z-score testing thresholds. One or more modZ:delta "
             "threshold pairs, e.g. 3.5:0.1 5:0.1. "
             "Defaults to 3.5:0.1 if the flag is given with no values.")
    p.add_argument("--gtf",                        default=None,
                   help="GTF/GTF.gz. If provided, adds junction_type column.")
    p.add_argument("--threads",                    type=int,   default=1)
    return p.parse_args()


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


def parse_gtf_junctions(
    gtf_path: str,
    gene_names: List[str],
) -> Dict[str, Dict]:
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
            start = int(parts[3]) - 1
            end   = int(parts[4])
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


def _classify_events_for_metric(
    metric_col:         str,
    outlier_set:        set,
    grp_coords:         Dict[Tuple, Dict[str, Tuple[int, int]]],
    strand_map:         Dict[str, str],
    padj_threshold:     float,
    delta_threshold:    float,
    positive_delta_set: set,
    delta_map:          Dict[Tuple, float],
    canonical_set:      Optional[set] = None,
) -> Dict[Tuple, List[str]]:
    """canonical_set: (gene, junction) pairs that are canonical. A skipping/inclusion or
    alt-splice-site group is only labelled if at least one of its junctions is in it. None
    disables the requirement (no GTF given)."""
    allowed_events = _METRIC_EVENTS.get(metric_col, [])
    if not allowed_events:
        return {}, []
    groups: List[Tuple] = []      # junction keys called together (one splicing group each)

    def _has_canonical(gene, *jxns):
        return canonical_set is None or any((gene, j) in canonical_set for j in jxns)

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
        if d > 0: return 1
        if d < 0: return -1
        return 0

    for grp_key, outlier_jxns in outlier_by_grp.items():
        sample, gene, phasing = grp_key
        strand = strand_map.get(gene, "+")
        coords = grp_coords[grp_key]

        for jxn in outlier_jxns:
            if jxn not in coords:
                continue
            ss1, ss2 = coords[jxn]
            five_ss  = ss1 if strand == "+" else ss2
            three_ss = ss2 if strand == "+" else ss1

            row_key  = (sample, gene, phasing, jxn)
            is_pos_delta = row_key in positive_delta_set

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
                for partner_jxn in outlier_jxns:
                    if partner_jxn == jxn or partner_jxn not in coords:
                        continue
                    p_ss1, p_ss2 = coords[partner_jxn]
                    p_five  = p_ss1 if strand == "+" else p_ss2
                    p_three = p_ss2 if strand == "+" else p_ss1
                    d_partner = _delta(sample, gene, phasing, partner_jxn)
                    s_partner = _sign(d_partner)
                    opposite = (s_jxn != 0 and s_partner != 0 and s_jxn != s_partner)

                    if not _has_canonical(gene, jxn, partner_jxn):
                        continue

                    if alt_3ss_label:
                        if p_five == five_ss and p_three != three_ss and opposite:
                            _add(sample, gene, phasing, jxn,         alt_3ss_label)
                            _add(sample, gene, phasing, partner_jxn, alt_3ss_label)
                            groups.append(((sample, gene, phasing, jxn), (sample, gene, phasing, partner_jxn)))

                    if alt_5ss_label:
                        if p_three == three_ss and p_five != five_ss and opposite:
                            _add(sample, gene, phasing, jxn,         alt_5ss_label)
                            _add(sample, gene, phasing, partner_jxn, alt_5ss_label)
                            groups.append(((sample, gene, phasing, jxn), (sample, gene, phasing, partner_jxn)))

            skip_label    = next((e for e in allowed_events if e.startswith("exon_skipping")),  None)
            incl_label    = next((e for e in allowed_events if e.startswith("exon_inclusion")), None)
            if skip_label or incl_label:
                d_long = _delta(sample, gene, phasing, jxn)
                s_long = _sign(d_long)
                if s_long == 0:
                    continue
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
                        if s_left == 0 or s_right == 0:
                            continue
                        if s_left == s_long or s_right == s_long:
                            continue
                        if not _has_canonical(gene, jxn, jl, jr):
                            continue
                        if skip_label and s_long > 0:
                            for j in (jxn, jl, jr):
                                _add(sample, gene, phasing, j, skip_label)
                            groups.append(tuple((sample, gene, phasing, j) for j in (jxn, jl, jr)))
                        if incl_label and s_long < 0:
                            for j in (jxn, jl, jr):
                                _add(sample, gene, phasing, j, incl_label)
                            groups.append(tuple((sample, gene, phasing, j) for j in (jxn, jl, jr)))

    return dict(event_map), groups


def _resolve_psi_events(events: set) -> set:
    """At most one PSI label and one PSI_approx label per junction, by priority:
    exon_skipping/exon_inclusion > alt_5ss/alt_3ss. A tie at the winning level (skipping AND
    inclusion, or alt_5ss AND alt_3ss) becomes "complex" ("complex_approx" for PSI_approx).
    Labels from other metrics (5ss_IR, 3ss_IR, full_IR, IPA) are left as they are."""
    out = set(events)
    for sfx in ("", "_approx"):
        skip_incl = out & {f"exon_skipping{sfx}", f"exon_inclusion{sfx}"}
        alt       = out & {f"alt_5ss{sfx}", f"alt_3ss{sfx}"}
        top = skip_incl or alt
        if not top:
            continue
        out -= skip_incl | alt
        out.add(f"complex{sfx}" if len(top) > 1 else next(iter(top)))
    return out


def classify_all_events(
    sig_df:              pd.DataFrame,
    final_df:            pd.DataFrame,
    strand_map:          Dict[str, str],
    effect_threshold:    float,
    threads:             int,
    has_ipa:             bool,
    computed_metrics:    List[str],
    effect_col_fn,
    unreliable_hap_outliers: Dict[str, set],
) -> pd.DataFrame:
    if sig_df.empty:
        sig_df = sig_df.copy()
        sig_df["event_type"] = ""
        return sig_df

    active_metrics = [mc for mc in computed_metrics if _METRIC_EVENTS.get(mc)]
    if not active_metrics:
        sig_df = sig_df.copy(); sig_df["event_type"] = "other"; return sig_df

    all_genes  = sorted(final_df["gene"].unique())
    n_genes    = len(all_genes)
    batch_size = max(1, n_genes // (threads * 4)) if threads > 1 else n_genes
    gene_batches = [all_genes[i:i + batch_size]
                    for i in range(0, n_genes, batch_size)]

    metric_outlier_sets:   Dict[str, set] = {}
    metric_pos_delta_sets: Dict[str, set] = {}
    metric_delta_maps:     Dict[str, Dict[Tuple, float]] = {}

    _ir_ipa_metrics = frozenset((
        "5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio", "junction_IPA_ratio"
    ))
    _has_jxn_type = "junction_type" in sig_df.columns

    for mc in computed_metrics:
        ocol    = f"outlier_{mc}"
        delta_c = effect_col_fn(mc)
        if ocol not in sig_df.columns:
            metric_outlier_sets[mc]   = set()
            metric_pos_delta_sets[mc] = set()
        else:
            mask = sig_df[ocol].astype(bool)

            if mc in _ir_ipa_metrics:
                if delta_c in sig_df.columns:
                    dv   = pd.to_numeric(sig_df[delta_c], errors="coerce")
                    mask = mask & dv.ge(effect_threshold)
                if _has_jxn_type:
                    mask = mask & sig_df["junction_type"].eq("canonical")

            metric_outlier_sets[mc] = set(zip(
                sig_df.loc[mask, "sample"], sig_df.loc[mask, "gene"],
                sig_df.loc[mask, "phasing"], sig_df.loc[mask, "junction"],
            ))
            if delta_c in sig_df.columns:
                dv       = pd.to_numeric(sig_df[delta_c], errors="coerce")
                pos_mask = mask & dv.ge(effect_threshold)
            else:
                pos_mask = mask
            metric_pos_delta_sets[mc] = set(zip(
                sig_df.loc[pos_mask, "sample"], sig_df.loc[pos_mask, "gene"],
                sig_df.loc[pos_mask, "phasing"], sig_df.loc[pos_mask, "junction"],
            ))
        # Deltas are only ever looked up for this metric's outlier keys (the junction itself and
        # its partner junctions, which are drawn from the same outlier set), so only those rows
        # of the full scored table are needed -- same result as mapping every row, far less work.
        out_keys = metric_outlier_sets.get(mc, set())
        if delta_c not in final_df.columns or not out_keys:
            metric_delta_maps[mc] = {}
        else:
            sub   = final_df[final_df["junction"].isin({k[3] for k in out_keys})]
            dv    = pd.to_numeric(sub[delta_c], errors="coerce")
            valid = dv.notna()
            metric_delta_maps[mc] = {
                k: v for k, v in zip(
                    zip(sub.loc[valid, "sample"], sub.loc[valid, "gene"],
                        sub.loc[valid, "phasing"], sub.loc[valid, "junction"]),
                    dv[valid],
                ) if k in out_keys
            }

    # Splice-site coordinates, likewise only for outlier keys: classification only looks up the
    # coordinates of outlier junctions (coords[jxn] / "j in outlier_jxns"). Parsed from the
    # junction ID ("chr_ss1_ss2"), exactly as the per-row loop over final_df used to.
    grp_coords_by_gene: Dict[str, Dict] = {g: {} for g in all_genes}
    for sample, g, phasing, jxn in set().union(*metric_outlier_sets.values()):
        parts = jxn.split("_")
        grp_coords_by_gene.setdefault(g, {}).setdefault((sample, g, phasing), {})[jxn] = (
            int(parts[-2]), int(parts[-1]))

    # (gene, junction) pairs that are canonical -- every junction in a skipping/inclusion or
    # alt-splice-site group must come from the outlier sets, so the outlier rows are enough.
    if _has_jxn_type:
        is_can = sig_df["junction_type"].eq("canonical")
        canonical_pairs: Optional[set] = set(zip(sig_df.loc[is_can, "gene"], sig_df.loc[is_can, "junction"]))
    else:
        canonical_pairs = None
        print("[WARNING] No junction_type column (no --gtf): splicing events are not required "
              "to include a canonical junction.")

    def _slice_batch(gene_batch, mc):
        gene_set = set(gene_batch)
        gc = {}
        for g in gene_batch:
            gc.update(grp_coords_by_gene.get(g, {}))
        os_ = {k for k in metric_outlier_sets.get(mc, set())   if k[1] in gene_set}
        ps_ = {k for k in metric_pos_delta_sets.get(mc, set()) if k[1] in gene_set}
        dm_ = {k: v for k, v in metric_delta_maps.get(mc, {}).items() if k[1] in gene_set}
        cs_ = None if canonical_pairs is None else {k for k in canonical_pairs if k[0] in gene_set}
        return gc, os_, ps_, dm_, cs_

    all_event_maps: List[Dict[Tuple, List[str]]] = []
    all_groups: List[Tuple] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futs = {}
        for mc in active_metrics:
            for batch in gene_batches:
                gc, os_, ps_, dm_, cs_ = _slice_batch(batch, mc)
                if not os_:
                    continue
                fut = ex.submit(
                    _classify_events_for_metric, mc,
                    os_, gc, strand_map, effect_threshold, effect_threshold, ps_, dm_, cs_,
                )
                futs[fut] = (mc, batch[0])
        for fut in concurrent.futures.as_completed(futs):
            mc, g0 = futs[fut]
            try:
                emap, grps = fut.result()
                all_event_maps.append(emap)
                all_groups.extend(grps)
            except Exception as e:
                print(f"[WARNING] Classification error for {mc} (batch @{g0}): {e}")
                traceback.print_exc()

    merged: Dict[Tuple, set] = defaultdict(set)
    for event_map in all_event_maps:
        for key, events in event_map.items():
            merged[key].update(events)

    # Event IDs: junctions linked by any splicing group (any metric, before labels are resolved) in
    # the same sample, gene and phasing form one event; a junction with only an IR/IPA label is an
    # event on its own. Numbered E1, E2, ... per sample and gene (bulk, hap1, hap2, then position).
    parent: Dict[Tuple, Tuple] = {}

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for grp in all_groups:
        for k in grp:
            parent.setdefault(k, k)
        for k in grp[1:]:
            parent[_find(k)] = _find(grp[0])
    for key, evs in merged.items():
        if evs:
            parent.setdefault(key, key)
    members: Dict[Tuple, List[Tuple]] = defaultdict(list)
    for k in parent:
        members[_find(k)].append(k)
    _phase_rank = {"bulk": 0, "hap1": 1, "hap2": 2}
    comps = sorted(members.values(), key=lambda m: (
        str(m[0][0]), str(m[0][1]), _phase_rank.get(m[0][2], 3),
        min(int(k[3].split("_")[-2]) for k in m)))
    event_id: Dict[Tuple, str] = {}
    counter: Dict[Tuple, int] = defaultdict(int)
    for m in comps:
        counter[(m[0][0], m[0][1])] += 1
        for k in m:
            event_id[k] = f"E{counter[(m[0][0], m[0][1])]}"

    sig_df = sig_df.copy()
    sig_records = sig_df[["sample", "gene", "phasing", "junction"]].to_dict("records")
    event_strs, event_ids = [], []
    for rec in sig_records:
        key = (rec["sample"], rec["gene"], rec["phasing"], rec["junction"])
        evs = _resolve_psi_events(merged.get(key, set()))
        event_strs.append(",".join(sorted(evs)) if evs else "none")
        event_ids.append(event_id.get(key, ".") if evs else ".")
    sig_df["event_type"] = event_strs
    sig_df["event_id"] = event_ids
    return sig_df


def _per_metric_output_cols(metric_col: str) -> List[str]:
    return [
        f"n_{metric_col}",
        f"alpha_{metric_col}", f"beta_{metric_col}", f"expected_{metric_col}",
        f"p1_{metric_col}", f"p99_{metric_col}",
        f"delta_{metric_col}", f"p_value_{metric_col}", f"padj_{metric_col}",
        f"median_{metric_col}", f"mad_{metric_col}", f"modz_{metric_col}",
    ]


_OUTPUT_COLS = [
    "sample", "gene", "gene_rank", "region", "phasing", "junction", "5ss", "3ss",
    "junction_type",
    "junction_usage",
    "junction_read_diversity",
    "junction_coverage_approx",
    "junction_PSI_approx", "rescaled_junction_PSI_approx",
    *_per_metric_output_cols("junction_PSI_approx"),
    "junction_coverage",
    "junction_PSI", "rescaled_junction_PSI",
    *_per_metric_output_cols("junction_PSI"),
    "5ss_usage", "5ss_coverage",
    "5ss_IR_ratio", "rescaled_5ss_IR_ratio",
    *_per_metric_output_cols("5ss_IR_ratio"),
    "3ss_usage", "3ss_coverage",
    "3ss_IR_ratio", "rescaled_3ss_IR_ratio",
    *_per_metric_output_cols("3ss_IR_ratio"),
    "junction_full_IR_count",
    "junction_full_IR_ratio", "rescaled_junction_full_IR_ratio",
    *_per_metric_output_cols("junction_full_IR_ratio"),
    "junction_IPA_count",
    "junction_IPA_ratio", "rescaled_junction_IPA_ratio",
    *_per_metric_output_cols("junction_IPA_ratio"),
]


def select_output_columns(df: pd.DataFrame, *_) -> pd.DataFrame:
    present = [c for c in _OUTPUT_COLS if c in df.columns]
    return df[present]


def main() -> None:
    print("\n" + "*"*80)
    print("  Cohort Junction Outlier Identification")
    print("*"*80 + "\n")

    args = parse_args()

    if args.bb_thresholds is not None:
        expected_method = "beta_binomial"
        threshold_specs: List = []
        for tok in args.bb_thresholds:
            try:
                p_str, d_str = tok.split(":")
                threshold_specs.append((float(p_str), float(d_str)))
            except Exception:
                raise ValueError(f"Invalid threshold format '{tok}'. Expected padj:delta, e.g. 0.01:0.1")
        if not threshold_specs:
            print("[INFO] --bb-thresholds given with no values; will apply no outlier thresholds.")
    else:
        expected_method = "modified_zscore"
        z_vals = args.z_thresholds if args.z_thresholds else ["3.5:0.1"]
        threshold_specs = []
        for tok in z_vals:
            try:
                z_str, d_str = tok.split(":")
                threshold_specs.append((float(z_str), float(d_str)))
            except Exception:
                raise ValueError(f"Invalid threshold format '{tok}'. Expected modZ:delta, e.g. 3.5:0.1")

    alias_map     = parse_alias_map(args.alias_map)
    prefix        = args.outprefix.rstrip("/")
    outdir        = os.path.dirname(os.path.abspath(prefix))
    prefix_name   = os.path.basename(prefix)
    approx_only   = args.approx

    # Find results directory written by 8B
    results_dirs = _glob.glob(os.path.join(outdir, f"{prefix_name}_results_*"))
    if not results_dirs:
        raise FileNotFoundError(
            f"No results directory found matching {prefix_name}_results_* under {outdir}. "
            "Did fit_and_score_cohort_junctions.py (8B) run successfully?"
        )

    # Read method from file
    method = None
    base_results = None
    for rd in results_dirs:
        method_file = os.path.join(rd, "method.txt")
        if os.path.exists(method_file):
            with open(method_file) as fh:
                method = fh.read().strip()
            base_results = rd
            break

    if method is None:
        raise FileNotFoundError(
            f"Could not find method.txt in any results directory under {outdir}. "
            "Did fit_and_score_cohort_junctions.py (8B) run successfully?"
        )

    if method != expected_method:
        raise ValueError(
            f"Method mismatch: results directory has method='{method}' but "
            f"CLI flags indicate method='{expected_method}'. "
            "Use matching --bb-thresholds / --z-thresholds flags."
        )

    print(f"Method (from scoring run): {method}")

    # Read computed_metrics
    metrics_file = os.path.join(base_results, "computed_metrics.txt")
    if os.path.exists(metrics_file):
        with open(metrics_file) as fh:
            computed_metrics = [l.strip() for l in fh if l.strip()]
    else:
        computed_metrics = []

    # Read all_scored_results.tsv
    all_scored_path = os.path.join(base_results, "all_scored_results.tsv")
    if not os.path.exists(all_scored_path):
        raise FileNotFoundError(
            f"Combined scored results not found: {all_scored_path}. "
            "Did fit_and_score_cohort_junctions.py (8B) run successfully?"
        )

    print(f"Loading scored results from {all_scored_path} ...")
    t_load = time.time()
    final_df = pd.read_csv(all_scored_path, sep="\t", dtype=str)
    for col in final_df.columns:
        try:
            final_df[col] = pd.to_numeric(final_df[col])
        except (ValueError, TypeError):
            pass
    print(f"  Loaded {len(final_df):,} rows ({time.time()-t_load:.2f}s)")

    # If computed_metrics is empty, infer from columns
    if not computed_metrics:
        test_prefix = "p_value_" if method == "beta_binomial" else "modz_"
        computed_metrics = [c.replace(test_prefix, "") for c in final_df.columns
                            if c.startswith(test_prefix)]
        print(f"  Inferred computed_metrics from columns: {computed_metrics}")

    # Load BED for strand_map
    gene_info = load_bed(args.bed)
    strand_map = {g: info[2] for g, info in gene_info.items()}

    # Load GTF if provided
    gtf_junctions: Optional[Dict[str, Dict]] = None
    if args.gtf:
        gene_names = list(final_df["gene"].unique()) if "gene" in final_df.columns else []
        print(f"\nParsing GTF for annotated junctions ...")
        t_gtf = time.time()
        gtf_junctions = parse_gtf_junctions(args.gtf, gene_names)
        n_matched = sum(1 for g in gene_names if g in gtf_junctions)
        print(f"  -> matched {n_matched}/{len(gene_names)} genes ({time.time()-t_gtf:.2f}s)")
    else:
        print("[INFO] --gtf not provided; junction_type column will not be added.")

    tmp_dir = os.path.join(outdir, f"{prefix_name}_tmp")

    def _threshold_dir(subdir_name: str) -> str:
        d = os.path.join(outdir, f"{prefix_name}_{subdir_name}")
        os.makedirs(d, exist_ok=True)
        return d

    def _thr_subdir_name(thr) -> str:
        if method == "beta_binomial":
            padj_threshold, effect_threshold = thr
            return f"padj{padj_threshold}_delta{effect_threshold}"
        z_threshold, effect_threshold = thr
        return f"z{z_threshold}_delta{effect_threshold}"

    def _empty_outlier_outputs(reason: str) -> None:
        print(f"\n[WARNING] {reason}")
        empty_cols = _OUTPUT_COLS + ["event_type"]
        for thr in threshold_specs:
            thr_dir  = _threshold_dir(_thr_subdir_name(thr))
            empty_df = pd.DataFrame(columns=empty_cols)
            empty_df.to_csv(os.path.join(thr_dir, f"{prefix_name}_outliers.tsv"), sep="\t", index=False)
            empty_df.to_csv(os.path.join(thr_dir, f"{prefix_name}_outliers_filtered.tsv"), sep="\t", index=False)
            empty_df.to_csv(os.path.join(thr_dir, f"{prefix_name}_outliers_alias.tsv"), sep="\t", index=False)
        print("  Wrote empty outlier file(s) so downstream outputs still exist.")

    if final_df.empty or not computed_metrics:
        _empty_outlier_outputs("Scored results are empty or no metrics were computed.")
        print("\nDone (nothing to do).")
        return

    _IR_IPA_METRICS = frozenset((
        "5ss_IR_ratio", "3ss_IR_ratio", "junction_full_IR_ratio", "junction_IPA_ratio"
    ))

    for thr in threshold_specs:
        t_thr = time.time()

        if method == "beta_binomial":
            padj_threshold, effect_threshold = thr
            stat_label = "delta"
            thr_desc   = f"padj <= {padj_threshold}, |delta| >= {effect_threshold}"

            def _effect_col(mc): return f"delta_{mc}"

            def _metric_cols_ok(df, mc):
                return f"padj_{mc}" in df.columns and _effect_col(mc) in df.columns

            def _outlier_mask(df, mc):
                if not _metric_cols_ok(df, mc):
                    return pd.Series(False, index=df.index)
                padj_v  = pd.to_numeric(df[f"padj_{mc}"], errors="coerce")
                delta_v = pd.to_numeric(df[_effect_col(mc)], errors="coerce")
                return padj_v.le(padj_threshold) & delta_v.abs().ge(effect_threshold)
        else:
            z_threshold, effect_threshold = thr
            stat_label = "delta"
            thr_desc   = f"|modZ| >= {z_threshold}, |delta| >= {effect_threshold}"

            def _effect_col(mc): return f"delta_{mc}"

            def _metric_cols_ok(df, mc):
                return f"modz_{mc}" in df.columns and _effect_col(mc) in df.columns

            def _outlier_mask(df, mc):
                if not _metric_cols_ok(df, mc):
                    return pd.Series(False, index=df.index)
                zv = pd.to_numeric(df[f"modz_{mc}"], errors="coerce")
                delta_v = pd.to_numeric(df[_effect_col(mc)], errors="coerce")
                return zv.abs().ge(z_threshold) & delta_v.abs().ge(effect_threshold)

        def _pos_mask(df, mc):
            ec = _effect_col(mc)
            if ec not in df.columns:
                return pd.Series(False, index=df.index)
            ev = pd.to_numeric(df[ec], errors="coerce")
            return ev.ge(effect_threshold)

        thr_dir = _threshold_dir(_thr_subdir_name(thr))
        print(f"\n{'='*70}")
        print(f"  Threshold: {thr_desc}")
        print(f"  Output: {thr_dir}")
        print(f"{'='*70}")

        outlier_mask = pd.Series(False, index=final_df.index)
        for mc in computed_metrics:
            outlier_mask |= _outlier_mask(final_df, mc)

        sig_df = final_df[outlier_mask].copy()

        print(f"  Identifying outliers with {thr_desc} ...")
        for mc in computed_metrics:
            if not _metric_cols_ok(final_df, mc): continue
            mask_mc = _outlier_mask(final_df, mc)
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

        unreliable_hap_outliers: Dict[str, set] = {}

        if not sig_df.empty:
            key_col = list(zip(sig_df["sample"], sig_df["junction"]))
            sig_df["_key"] = key_col
            # Only (sample, junction) keys that have an outlier row are ever looked up below, so
            # the bulk/hap1/hap2 pivot is built for those keys only (same values, far fewer rows).
            _sig_keys = pd.MultiIndex.from_tuples(sorted(set(key_col)), names=["sample", "junction"])
            _phase_base = final_df[
                final_df["phasing"].isin(["bulk", "hap1", "hap2"]).to_numpy()
                & pd.MultiIndex.from_arrays([final_df["sample"], final_df["junction"]]).isin(_sig_keys)
            ]

            for mc in computed_metrics:
                rescaled_c = f"rescaled_{mc}"
                if rescaled_c not in final_df.columns:
                    unreliable_hap_outliers[mc] = set()
                    continue
                phasing_df = _phase_base[
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

                sandwiched = ((h1 <= b) & (b <= h2)) | ((h2 <= b) & (b <= h1))

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
            for mc in computed_metrics:
                unreliable = unreliable_hap_outliers.get(mc, set())
                if not unreliable or not _metric_cols_ok(sig_df, mc): continue
                hap_rows_mc = sig_df[sig_df["phasing"].isin(["hap1", "hap2"])]
                passes_mc = _outlier_mask(hap_rows_mc, mc)
                keys = list(zip(
                    hap_rows_mc.loc[passes_mc, "sample"],
                    hap_rows_mc.loc[passes_mc, "junction"],
                ))
                n_tot = sum(1 for k in keys if k in unreliable)
                print(f"       {mc}: removed {n_tot:,} unreliable haplotype outlier rows")

            n_total_jxns = sig_df["junction"].nunique() if not sig_df.empty else 0
            print(f"  {len(sig_df):,} total outlier rows ({n_total_jxns} unique junctions)")
        else:
            unreliable_hap_outliers = {mc: set() for mc in computed_metrics}

        sig_df = sig_df.copy()
        for mc in computed_metrics:
            ocol = f"outlier_{mc}"
            if not _metric_cols_ok(sig_df, mc):
                sig_df[ocol] = False; continue
            passes     = _outlier_mask(sig_df, mc)
            unreliable = unreliable_hap_outliers.get(mc, set())
            if unreliable:
                is_hap        = sig_df["phasing"].isin(["hap1", "hap2"])
                keys          = list(zip(sig_df["sample"], sig_df["junction"]))
                is_unreliable = pd.Series([k in unreliable for k in keys], index=sig_df.index)
                sig_df[ocol]  = passes & ~(is_hap & is_unreliable)
            else:
                sig_df[ocol] = passes

        _has_jxn_type = "junction_type" in sig_df.columns
        for mc in computed_metrics:
            ocol = f"outlier_{mc}"
            if ocol not in sig_df.columns: continue
            mask_mc = sig_df[ocol].astype(bool)
            n_rows  = int(mask_mc.sum())
            n_jxns  = int(sig_df[mask_mc]["junction"].nunique())
            is_5ss  = mc == "5ss_IR_ratio" and "5ss" in sig_df.columns
            is_3ss  = mc == "3ss_IR_ratio" and "3ss" in sig_df.columns
            ss_col  = "5ss" if is_5ss else ("3ss" if is_3ss else None)
            ss_lbl  = "5ss" if is_5ss else ("3ss" if is_3ss else None)
            if mc in _IR_IPA_METRICS and _effect_col(mc) in sig_df.columns:
                pos_mask = mask_mc & _pos_mask(sig_df, mc)
                n_pos    = int(sig_df[pos_mask]["junction"].nunique())
                ss_total = f", {int(sig_df[mask_mc][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                ss_pos   = f", {int(sig_df[pos_mask][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                if _has_jxn_type:
                    can_ann = sig_df["junction_type"].isin(["canonical"])
                    n_can   = int(sig_df[pos_mask & can_ann]["junction"].nunique())
                    ss_can  = f", {int(sig_df[pos_mask & can_ann][ss_col].nunique())} unique {ss_lbl}" if ss_col else ""
                    print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions{ss_total} -> {n_pos} unique junctions{ss_pos} with {stat_label} > 0 -> {n_can} unique junctions{ss_can} canonical)")
                else:
                    print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions{ss_total} -> {n_pos} unique junctions{ss_pos} with {stat_label} > 0)")
            elif ss_col:
                n_ss = int(sig_df[mask_mc][ss_col].nunique())
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions, {n_ss} unique {ss_lbl})")
            else:
                print(f"       {mc}: {n_rows:,} rows ({n_jxns} unique junctions)")

        # Number of distinct outlier samples per (gene, junction) and metric. Only the outlier
        # rows ever receive these columns, so they are counted here and attached to sig_df below
        # (rather than merged into a full copy of the scored table first).
        n_sample_counts: Dict[str, Optional[pd.DataFrame]] = {}
        for mc in computed_metrics:
            n_col = f"n_sample_outlier_{mc}"
            ocol  = f"outlier_{mc}"
            if sig_df.empty or ocol not in sig_df.columns:
                n_sample_counts[n_col] = None
                continue
            n_sample_counts[n_col] = (
                sig_df[sig_df[ocol].astype(bool)]
                .groupby(["gene", "junction"])["sample"]
                .nunique().rename(n_col).reset_index()
            )

        n_sig = len(sig_df)
        print(f"  Outliers identified ({time.time() - t_thr:.2f}s)")
        print(f"  Classifying events ...")
        t_cls = time.time()
        if n_sig > 0:
            sig_df = classify_all_events(
                sig_df, final_df, strand_map,
                effect_threshold, args.threads, args.has_ipa,
                computed_metrics, _effect_col, unreliable_hap_outliers,
            )
        else:
            sig_df["event_type"] = "none"
        print(f"       done ({time.time()-t_cls:.2f}s)  "
              f"{sig_df['event_type'].value_counts().to_dict() if n_sig else {}}")

        if n_sig > 0:
            effect_cols = [c for c in sig_df.columns
                          if c.startswith("delta_") or c.startswith("modz_")]
            if effect_cols:
                effect_num = sig_df[effect_cols].apply(lambda col: pd.to_numeric(col, errors="coerce"))
                sig_df["_max_abs_delta"] = effect_num.abs().max(axis=1)

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
                            dc = _effect_col(mc)
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

        for n_col, counts in n_sample_counts.items():
            if counts is None:
                sig_df[n_col] = 0
            else:
                sig_df = sig_df.merge(counts, on=["gene", "junction"], how="left")
                sig_df[n_col] = sig_df[n_col].fillna(0).astype(int)

        outlier_anchor_col = "padj_" if method == "beta_binomial" else "modz_"
        outlier_tsv_cols = []
        for col in _OUTPUT_COLS:
            outlier_tsv_cols.append(col)
            for mc in computed_metrics:
                if col == f"{outlier_anchor_col}{mc}":
                    n_sc = f"n_sample_outlier_{mc}"
                    if n_sc in sig_df.columns:
                        outlier_tsv_cols.append(n_sc)
                    outlier_tsv_cols.append(f"outlier_{mc}")
        outlier_tsv_cols.append("event_type")
        outlier_tsv_cols.append("event_id")
        outlier_tsv_cols = [c for c in outlier_tsv_cols if c in sig_df.columns]

        t_out = time.time()
        out_jobs = []

        _outliers_data = sig_df[outlier_tsv_cols].copy()
        if "event_type" in sig_df.columns:
            filt = sig_df[sig_df["event_type"] != "none"]
            _outliers_filt_data = filt[outlier_tsv_cols].copy() if not filt.empty else pd.DataFrame(columns=outlier_tsv_cols)
        else:
            filt = sig_df.iloc[0:0]
            _outliers_filt_data = pd.DataFrame(columns=outlier_tsv_cols)

        def _write_outliers(data=_outliers_data):
            t0 = time.time()
            out = os.path.join(thr_dir, f"{prefix_name}_outliers.tsv")
            data.to_csv(out, sep="\t", index=False)
            print(f"  Outliers -> {out} ({time.time()-t0:.2f}s)")

        def _write_outliers_filtered(data=_outliers_filt_data):
            t0 = time.time()
            out_filt = os.path.join(thr_dir, f"{prefix_name}_outliers_filtered.tsv")
            data.to_csv(out_filt, sep="\t", index=False)
            print(f"  Outliers (filtered) -> {out_filt} ({time.time()-t0:.2f}s)")

        def _write_outliers_alias(data=_outliers_data):
            t0 = time.time()
            out_alias = os.path.join(thr_dir, f"{prefix_name}_outliers_alias.tsv")
            alias_data = data.copy()
            if "sample" in alias_data.columns:
                alias_data["sample"] = alias_data["sample"].apply(lambda s: resolve(s, alias_map))
            alias_data.to_csv(out_alias, sep="\t", index=False)
            print(f"  Outliers (alias) -> {out_alias} ({time.time()-t0:.2f}s)")


        out_jobs.append(_write_outliers)
        out_jobs.append(_write_outliers_filtered)
        out_jobs.append(_write_outliers_alias)

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
            computed_mets   = set(computed_metrics)
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
                print(f"  Outlier summary -- {thr_desc}")
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

        print(f"\n  Finished threshold: {thr_desc} ({time.time()-t_thr:.2f}s)")

    print("\nDone.")


if __name__ == "__main__":
    main()
