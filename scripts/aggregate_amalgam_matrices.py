
import argparse
import gzip
import re

import numpy as np
import pandas as pd

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all


# ---------------------------------------------------------------------------
# GTF helpers (same logic as normalize_amalgam_matrix.py)
# ---------------------------------------------------------------------------

_ATTR_RE_CACHE = {}


def _attr(attr_str, key):
    rx = _ATTR_RE_CACHE.get(key)
    if rx is None:
        rx = re.compile(key + r'\s+"([^"]+)"')
        _ATTR_RE_CACHE[key] = rx
    m = rx.search(attr_str)
    return m.group(1) if m else None


def _strip_ver(s):
    return re.sub(r'\.\d+$', '', s) if s else s


def load_gene_id_to_symbol(gtf_path):
    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    mapping = {}
    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            gid = _strip_ver(_attr(attrs, "gene_id"))
            gname = _attr(attrs, "gene_name") or _attr(attrs, "gene_symbol")
            if gid and gname:
                mapping[gid] = gname
    return mapping


def load_targeted_genes(bed_path):
    genes = []
    seen = set()
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            gene = cols[3]
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return sorted(genes)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-sample AMALGAM transcript quantification files into "
                    "cohort-wide raw count matrices (genome-wide and BED-panel-filtered), "
                    "written directly into by_amalgam/.")
    parser.add_argument('--infiles', required=True, nargs='+',
        help="Every sample's <sample>_transcript_quantification.tsv in this group.")
    parser.add_argument('--samples', required=True, nargs='+',
        help="Sample name for each --infiles entry, same order/length.")
    parser.add_argument('--outdir', required=True,
        help="Directory to write output files into (by_amalgam/).")
    parser.add_argument('--bed', required=True,
        help="This group's BED panel (column 4 = target gene symbols).")
    parser.add_argument('--gtf', required=True,
        help="Reference annotation GTF, for gene_id -> gene_name translation.")
    add_alias_map_arg(parser)
    args = parser.parse_args()
    if len(args.infiles) != len(args.samples):
        parser.error("--infiles and --samples must have the same number of entries")
    return args


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Note: read_csv's usecols does NOT reorder columns (they come back in file order,
    # which for AMALGAM is gene_id, transcript_id, count), so the key order is pinned
    # explicitly via KEY everywhere -- otherwise the index levels mismatch and every
    # reindex() silently returns NaN (-> an all-zero transcript matrix).
    KEY = ['transcript_id', 'gene_id']

    # ------------------------------------------------------------------
    # Pass 1: build the union transcript index across all samples
    # ------------------------------------------------------------------
    tx_index = None
    for f in args.infiles:
        ids = pd.read_csv(f, sep='\t', usecols=KEY)[KEY]
        idx = pd.MultiIndex.from_frame(ids)
        tx_index = idx.unique() if tx_index is None else tx_index.union(idx, sort=False)

    # ------------------------------------------------------------------
    # Pass 2: fill transcript matrix and accumulate per-gene counts
    # ------------------------------------------------------------------
    transcript_matrix = pd.DataFrame(index=tx_index)
    gene_tables = []
    for sample, f in zip(args.samples, args.infiles):
        df = pd.read_csv(f, sep='\t', usecols=KEY + ['count'])
        counts = df.set_index(KEY)['count']
        transcript_matrix[sample] = counts.reindex(tx_index)
        # Sanity check: every count in the input must land in the matrix
        if not np.isclose(transcript_matrix[sample].sum(skipna=True), counts.sum()):
            raise ValueError(f"{sample}: transcript counts did not align to the matrix index "
                             f"({transcript_matrix[sample].sum(skipna=True)} vs {counts.sum()})")
        gene = df.groupby('gene_id')['count'].sum().to_frame(name=sample)
        gene_tables.append(gene)
        print(f"Processed {sample}", flush=True)

    transcript_matrix = transcript_matrix.fillna(0).reset_index()
    gene_matrix_gid = pd.concat(gene_tables, axis=1).fillna(0)   # gene_id-keyed

    # ------------------------------------------------------------------
    # Translate gene_id -> gene symbol, collapse, filter to BED panel
    # ------------------------------------------------------------------
    gid_to_symbol = load_gene_id_to_symbol(args.gtf)
    targeted_genes = load_targeted_genes(args.bed)

    # Gene matrix: genome-wide, symbol-keyed
    symbol_index = [gid_to_symbol.get(_strip_ver(gid), gid) for gid in gene_matrix_gid.index]
    gene_matrix_all = gene_matrix_gid.copy()
    gene_matrix_all.index = symbol_index
    gene_matrix_all.index.name = "gene"
    gene_matrix_all = gene_matrix_all.groupby(level=0).sum()

    # Gene matrix: BED-panel genes only
    gene_matrix_filtered = gene_matrix_all.reindex(targeted_genes).fillna(0).astype(int)
    gene_matrix_filtered.index.name = "gene"

    # Transcript matrix: genome-wide (transcript_id/gene_id as columns)
    # Translate gene_id column to symbol in-place
    transcript_matrix_all = transcript_matrix.copy()
    transcript_matrix_all['gene_id'] = [
        gid_to_symbol.get(_strip_ver(gid), gid)
        for gid in transcript_matrix_all['gene_id']
    ]

    # Transcript matrix: BED-panel genes only
    transcript_matrix_filtered = transcript_matrix_all[
        transcript_matrix_all['gene_id'].isin(set(targeted_genes))
    ].copy()

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    import os
    os.makedirs(args.outdir, exist_ok=True)

    out_gene_all      = os.path.join(args.outdir, "gene_amalgam_matrix_raw_all_genes.tsv")
    out_gene_filtered = os.path.join(args.outdir, "gene_amalgam_matrix_raw.tsv")
    out_tx_all        = os.path.join(args.outdir, "transcript_amalgam_matrix_raw_all_genes.tsv")
    out_tx_filtered   = os.path.join(args.outdir, "transcript_amalgam_matrix_raw.tsv")

    gene_matrix_all.to_csv(out_gene_all, sep='\t')
    print(f"Done: {out_gene_all} written.", flush=True)

    gene_matrix_filtered.to_csv(out_gene_filtered, sep='\t')
    print(f"Done: {out_gene_filtered} written.", flush=True)

    transcript_matrix_all.to_csv(out_tx_all, sep='\t', index=False)
    print(f"Done: {out_tx_all} written.", flush=True)

    transcript_matrix_filtered.to_csv(out_tx_filtered, sep='\t', index=False)
    print(f"Done: {out_tx_filtered} written.", flush=True)


if __name__ == '__main__':
    main()
