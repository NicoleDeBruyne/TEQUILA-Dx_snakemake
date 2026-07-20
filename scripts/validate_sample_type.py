#!/usr/bin/env python3

import argparse
import gc
import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="PSI-based tissue identity classifier for targeted RNA-seq panels."
    )
    parser.add_argument(
        "--matrix-refs", required=True, nargs="+",
        help="Junction count matrices for reference tissues (one per tissue, TSV)"
    )
    parser.add_argument(
        "--ref-names", required=True, nargs="+",
        help="Display names for each reference matrix (same order as --matrix-refs)"
    )
    parser.add_argument(
        "--matrix-query", required=True, nargs="+",
        help="One or more junction count matrices for query samples (TSV, col1 = chr_ss1_ss2). "
             "PSI is always calculated from the count matrix itself."
    )
    parser.add_argument(
        "--query-names", required=True, nargs="+",
        help="Display names for each query matrix (same order as --matrix-query)."
    )
    parser.add_argument(
        "--ref-colors", default=None, nargs="+",
        help="Hex colors for each reference tissue (same order as --ref-names). "
             "E.g. '#2166ac' '#d6604d'. Falls back to built-in palette if omitted."
    )
    parser.add_argument(
        "--query-colors", default=None, nargs="+",
        help="Hex colors for each query matrix (same order as --query-names). "
             "E.g. '#111111' '#e31a1c'. Falls back to built-in palette if omitted."
    )
    parser.add_argument(
        "--bed", required=True,
        help="BED file defining targeted regions (chr, start, end)"
    )
    parser.add_argument(
        "--outprefix", required=True,
        help="Output file prefix"
    )
    parser.add_argument(
        "--min-coverage", type=int, default=20,
        help="Minimum denominator (ss1+ss2-jxn) to call a PSI value (default: 20)"
    )
    parser.add_argument(
        "--min-ref-samples", type=int, default=50,
        help="Minimum number of reference samples with valid PSI to retain a "
             "junction in a reference matrix (default: 50)"
    )
    parser.add_argument(
        "--n-variable-junctions", type=int, default=0,
        help="Override automatic elbow detection and use exactly this many "
             "top-variance junctions (0 = auto-detect via elbow, default: 0)"
    )
    args = parser.parse_args()

    if len(args.matrix_refs) != len(args.ref_names):
        parser.error("--matrix-refs and --ref-names must have the same number of entries")

    if len(args.query_names) != len(args.matrix_query):
        parser.error("--query-names must have the same number of entries as --matrix-query")

    if args.ref_colors is not None:
        if len(args.ref_colors) != len(args.ref_names):
            parser.error("--ref-colors must have the same number of entries as --ref-names")

    if args.query_colors is not None:
        if len(args.query_colors) != len(args.query_names):
            parser.error("--query-colors must have the same number of entries as --query-names")

    return args


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _normalize_junction_index(index: pd.Index) -> pd.Index:
    """Normalize junction IDs to chr_start_end.

    GTEx-derived matrices (--matrix-refs) use "chr:start-end:strand" (e.g.
    "chr1:11212-12009:+"), while this pipeline's own sample-derived matrices
    (--matrix-query) already use "chr_start_end" (e.g. "chr1_11212_12009") --
    see get_splice_junction_counts_by_region.py / make_junction_count_matrix.py.
    Only entries in the colon/dash format are rewritten; entries already in
    chr_start_end format (no ':') pass through unchanged, so this is safe to
    apply uniformly to both matrix types regardless of source.
    """
    idx = index.to_series(index=range(len(index))).astype(str)
    is_colon_fmt = idx.str.contains(':', regex=False)
    if not is_colon_fmt.any():
        return index

    colon_idx = idx[is_colon_fmt]
    chrom_coord = colon_idx.str.split(':', n=2, expand=True)
    coords = chrom_coord[1].str.split('-', n=1, expand=True)
    normalized = chrom_coord[0] + '_' + coords[0] + '_' + coords[1]

    idx = idx.copy()
    idx.loc[is_colon_fmt] = normalized
    return pd.Index(idx.values, name=index.name)


def read_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = _normalize_junction_index(df.index)
    df = df.apply(pd.to_numeric, errors="coerce", downcast="float")
    # Ensure float32 (downcast may still leave float64 for columns with NaN-only etc.)
    df = df.astype(np.float32)
    return df


def read_bed(path: str) -> pd.DataFrame:
    bed = pd.read_csv(
        path, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2], names=["chrom", "start", "end"]
    )
    bed["start"] = bed["start"].astype(int)
    bed["end"]   = bed["end"].astype(int)
    return bed


def build_bed_dict(bed: pd.DataFrame) -> dict:
    """
    Pre-build the chrom → (starts_array, ends_array) lookup used by
    filter_to_bed. Intervals are sorted by start for vectorized searching.
    Construct once and reuse across all matrices.
    """
    bed_dict = {}
    for chrom, sub in bed.groupby("chrom", sort=False):
        starts = np.sort(sub["start"].to_numpy())
        ends   = sub["end"].to_numpy()[np.argsort(sub["start"].to_numpy())]
        bed_dict[chrom] = (starts, ends)
    return bed_dict


def deduplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Append .1, .2, … to any duplicated column names."""
    seen = {}
    new_cols = []
    duplicates = set()
    for col in df.columns:
        if col in seen:
            duplicates.add(col)
            seen[col] += 1
            new_cols.append(f"{col}.{seen[col]}")
        else:
            seen[col] = 0
            new_cols.append(col)
    for dup in sorted(duplicates):
        print(f"  [WARNING] Duplicate sample name '{dup}' — suffixes appended (.1, .2, …)")
    df.columns = new_cols
    return df


# ---------------------------------------------------------------------------
# Junction parsing & BED overlap
# ---------------------------------------------------------------------------

def parse_junctions(index: pd.Index) -> pd.DataFrame:
    """Split chr_ss1_ss2 index into component columns (vectorized)."""
    split = index.to_series(index=range(len(index))).str.rsplit("_", n=2, expand=True)

    if split.shape[1] < 3:
        # Every entry failed to split into 3 parts
        print(f"  WARNING: {len(index)} junction IDs could not be parsed and will be skipped.")
        return pd.DataFrame(
            columns=["chr", "ss1", "ss2", "ss1_key", "ss2_key"],
            index=pd.Index([], name="junction"),
        )

    bad_mask = split[0].isna() | split[1].isna() | split[2].isna()
    n_bad = int(bad_mask.sum())
    if n_bad:
        print(f"  WARNING: {n_bad} junction IDs could not be parsed and will be skipped.")

    split = split.loc[~bad_mask.values]
    chrom = split[0].astype(str)
    ss1   = split[1].astype(str)
    ss2   = split[2].astype(str)

    out = pd.DataFrame({
        "chr":     chrom.values,
        "ss1":     ss1.values,
        "ss2":     ss2.values,
        "ss1_key": (chrom + "_" + ss1).values,
        "ss2_key": (chrom + "_" + ss2).values,
    }, index=index[~bad_mask.values])
    out.index.name = "junction"
    return out


def filter_to_bed(psi_df: pd.DataFrame, bed_dict: dict) -> pd.DataFrame:
    """
    Drop rows whose junction does not have BOTH splice sites within a BED
    interval on the same chromosome.  Operates on an already-computed PSI
    DataFrame so no count data needs to be kept around.

    Vectorized: for each chromosome, splice-site positions are checked
    against sorted (start, end) interval arrays using searchsorted, avoiding
    a per-row Python loop over hundreds of thousands of junctions.
    """
    jxn_info = parse_junctions(psi_df.index)

    def in_bed_vectorized(chrom_arr: np.ndarray, pos_arr: np.ndarray) -> np.ndarray:
        result = np.zeros(len(pos_arr), dtype=bool)
        for chrom, (starts, ends) in bed_dict.items():
            mask = (chrom_arr == chrom)
            if not mask.any():
                continue
            p = pos_arr[mask]
            # For each position, find the rightmost interval with start <= p
            idx = np.searchsorted(starts, p, side="right") - 1
            valid = idx >= 0
            hit = np.zeros(len(p), dtype=bool)
            if valid.any():
                hit[valid] = p[valid] < ends[idx[valid]]
            result[mask] = hit
        return result

    chrom_arr = jxn_info["chr"].to_numpy()
    ss1_arr   = jxn_info["ss1"].to_numpy(dtype=np.int64)
    ss2_arr   = jxn_info["ss2"].to_numpy(dtype=np.int64)

    keep_ss1 = in_bed_vectorized(chrom_arr, ss1_arr)
    keep_ss2 = in_bed_vectorized(chrom_arr, ss2_arr)
    keep = keep_ss1 & keep_ss2

    return psi_df.loc[jxn_info.index[keep]]


# ---------------------------------------------------------------------------
# PSI calculation
# ---------------------------------------------------------------------------

def compute_psi(counts: pd.DataFrame,
                min_coverage: int = 20,
                min_samples: int = 50) -> pd.DataFrame:
    """
    Compute PSI = jxn_count / (ss1_count + ss2_count - jxn_count).

    Cells where denominator < min_coverage are set to NaN.
    Junctions with fewer than min_samples valid (non-NaN) PSI values are
    dropped.  Pass min_samples=1 for query matrices (keep a junction if at
    least one sample has sufficient coverage).
    """
    jxn_info   = parse_junctions(counts.index)
    counts      = counts.loc[jxn_info.index]          # drop unparseable rows
    counts_arr  = counts.values.astype(np.float32)
    ss1_keys    = jxn_info["ss1_key"].values
    ss2_keys    = jxn_info["ss2_key"].values

    # Accumulate per-splice-site read totals across all junctions sharing that site
    ss_to_rows = defaultdict(list)
    for row_i, (k1, k2) in enumerate(zip(ss1_keys, ss2_keys)):
        ss_to_rows[k1].append(row_i)
        ss_to_rows[k2].append(row_i)

    ss_sum_cache = {
        ss: np.nansum(counts_arr[rows, :], axis=0)
        for ss, rows in ss_to_rows.items()
    }
    del ss_to_rows

    # Build denominator directly without materializing separate
    # ss1_counts / ss2_counts full-size arrays (saves 2x matrix memory).
    denominator = np.empty_like(counts_arr)
    for row_i, (k1, k2) in enumerate(zip(ss1_keys, ss2_keys)):
        denominator[row_i, :] = ss_sum_cache[k1] + ss_sum_cache[k2]
    del ss_sum_cache
    denominator -= counts_arr

    with np.errstate(invalid="ignore", divide="ignore"):
        psi = np.where(denominator >= min_coverage,
                       counts_arr / denominator,
                       np.nan).astype(np.float32)
    del counts_arr, denominator

    psi_df = pd.DataFrame(psi, index=counts.index, columns=counts.columns)
    del psi

    valid = psi_df.notna().sum(axis=1)
    psi_df = psi_df.loc[valid >= min_samples]

    del jxn_info
    gc.collect()

    return psi_df


# ---------------------------------------------------------------------------
# Per-matrix load-and-filter pipeline
# ---------------------------------------------------------------------------

def load_ref_matrix(path: str, name: str,
                    bed_dict: dict,
                    min_coverage: int,
                    min_ref_samples: int) -> pd.DataFrame:
    """
    Read one reference count matrix, compute PSI, apply coverage filter,
    then restrict to BED junctions.  The raw count matrix is released from
    memory before returning.
    """
    print(f"Reading in reference matrix: {path} ({name})...")
    counts = read_matrix(path)
    print(f"  {counts.shape[0]:,} junctions x {counts.shape[1]:,} samples")

    psi = compute_psi(counts, min_coverage=min_coverage,
                      min_samples=min_ref_samples)
    del counts   # release raw counts immediately
    gc.collect()

    print(f"  {len(psi):,} junctions with coverage ≥{min_coverage} in ≥{min_ref_samples} samples")

    psi = filter_to_bed(psi, bed_dict)
    print(f"  {len(psi):,} junctions located within BED regions")

    return psi


def load_query_matrix(path: str, name: str,
                      bed_dict: dict,
                      min_coverage: int,
                      seen_columns: set,
                      common_jxns: pd.Index,
                      variable_jxns: pd.Index) -> tuple[pd.DataFrame, dict]:
    """
    Read one query count matrix, compute PSI, apply coverage filter
    (keep junction if ≥1 sample has valid PSI), restrict to BED junctions,
    then report overlap with the reference common junction set and the
    tissue-discriminating variable junction set.
    Deduplicates column names against previously seen names.

    Returns (psi_df, sample_to_query mapping for this matrix).
    """
    print(f"Reading in query matrix: {path}...")
    counts = read_matrix(path)

    # Deduplicate column names against all previously loaded query samples
    new_cols = []
    col_seen_local = {}
    for col in counts.columns:
        if col in seen_columns or col in col_seen_local:
            suffix = col_seen_local.get(col, 0) + 1
            col_seen_local[col] = suffix
            new_col = f"{col}.{suffix}"
            print(f"  [WARNING] Duplicate sample name '{col}' — renamed to '{new_col}'")
            new_cols.append(new_col)
        else:
            col_seen_local[col] = 0
            new_cols.append(col)
    counts.columns = new_cols
    seen_columns.update(new_cols)

    print(f"  {counts.shape[0]:,} junctions x {counts.shape[1]:,} samples")

    psi = compute_psi(counts, min_coverage=min_coverage, min_samples=1)
    del counts
    gc.collect()

    print(f"  {len(psi):,} junctions with coverage ≥{min_coverage} in ≥1 sample")

    psi = filter_to_bed(psi, bed_dict)
    print(f"  {len(psi):,} junctions located within BED regions")

    n_in_ref = psi.index.isin(common_jxns).sum()
    print(f"  {n_in_ref:,} junctions with adequate coverage in all reference matrices")

    n_discriminating = psi.index.isin(variable_jxns).sum()
    print(f"  {n_discriminating:,} discriminating junctions")

    # Downstream code (PCA, distance scoring) only ever uses variable_jxns,
    # so retain just those rows. This is the dominant memory saving: full
    # BED-filtered query matrices can have 10,000-150,000+ junctions, while
    # variable_jxns is typically ~100.
    psi = psi.reindex(variable_jxns)
    gc.collect()

    sample_to_query = {col: name for col in psi.columns}
    return psi, sample_to_query


# ---------------------------------------------------------------------------
# Tissue-discriminating junction selection
# ---------------------------------------------------------------------------

def _detect_elbow(values: np.ndarray) -> int:
    """
    Kneedle algorithm: index of maximum perpendicular distance from the line
    joining the first and last points of the (normalised) curve.
    """
    n = len(values)
    if n < 3:
        return n - 1
    if values.max() - values.min() < 1e-10:
        return n - 1

    x = np.linspace(0, 1, n)
    y = (values - values.min()) / (values.max() - values.min())
    dx, dy = x[-1] - x[0], y[-1] - y[0]
    norm = np.sqrt(dx**2 + dy**2)
    dist = np.abs(dy * x - dx * y + x[-1] * y[0] - y[-1] * x[0]) / norm
    return int(np.argmax(dist))


def select_variable_junctions(avg_psi_per_tissue: pd.DataFrame,
                               n_override: int = 0,
                               outprefix: str = "") -> pd.Index:
    """
    Select tissue-discriminating junctions by elbow detection on the full
    ranked-variance curve of all junctions present in all reference matrices.

    n_override > 0 : skip elbow detection and use exactly that many junctions.
    """
    variance = avg_psi_per_tissue.dropna().var(axis=1).sort_values(ascending=False)
    vals     = variance.values

    if n_override > 0:
        n_select = min(n_override, len(vals))
    else:
        n_select = max(_detect_elbow(vals), 1)

    selected = variance.index[:n_select]
    print(f"  {n_select:,} junctions selected")

    if len(selected) < 10:
        print("  [WARNING] Fewer than 10 junctions selected. "
              "Consider --n-variable-junctions to override.")

    # Diagnostic scree plot
    if outprefix:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(np.arange(1, len(vals) + 1), vals, color="#4c8fca", lw=1.2)
        ax.axvline(n_select, color="#d95d5b", lw=1.2, ls="--",
                   label=f"selected N = {n_select:,}")
        ax.scatter([n_select], [vals[n_select - 1]], color="#d95d5b", s=50, zorder=4)
        ax.set_xlabel("Junctions ranked by variance", fontsize=9)
        ax.set_ylabel("Cross-tissue PSI variance", fontsize=9)
        ax.set_title("Ranked variance — elbow detection", fontsize=10)
        ax.legend(fontsize=8)
        ax.tick_params(labelsize=8)
        fig.tight_layout()
        path = f"{outprefix}_variable_junction_scree.pdf"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {path}\n")

    return selected


# ---------------------------------------------------------------------------
# Color helpers
# ---------------------------------------------------------------------------

def _tissue_palette(tissues: list, colors: list | None = None) -> dict:
    """Distinct color per reference tissue.

    If *colors* is provided it must be the same length as *tissues*; each
    entry may be any Matplotlib-compatible color string (hex, named, etc.).
    Missing or None entries fall back to the built-in defaults.
    """
    defaults=[
        "#912321","#002b58","#1e662a","#c23637","#0068a9","#3d892e",
        "#d95d5b","#4c8fca","#57aa3e","#ea9a9c","#91c4e9","#95c36e",
    ]
    palette = {}
    for i, tissue in enumerate(tissues):
        user_color = (colors[i] if colors and i < len(colors) else None)
        if user_color:
            user_color = user_color.strip("'\"")
        palette[tissue] = user_color if user_color else defaults[i % len(defaults)]
    return palette


def _query_palette(query_names: list, colors: list | None = None) -> dict:
    """Distinct color per query matrix (separate palette from tissue colors).

    If *colors* is provided it must be the same length as *query_names*; each
    entry may be any Matplotlib-compatible color string (hex, named, etc.).
    Missing or None entries fall back to the built-in defaults.
    """
    defaults = [
        "#ae450b","#6e2769","#005d6e","#926d17",
        "#ea6302","#9d4588","#009099","#c69528",
        "#f48f3e","#b271ab","#42b4b5","#e8c048",
        "#fbb875","#cca0ca","#8cccce","#f4d77e",
    ]
    palette = {}
    for i, name in enumerate(query_names):
        user_color = (colors[i] if colors and i < len(colors) else None)
        if user_color:
            user_color = user_color.strip("'\"")
        palette[name] = user_color if user_color else defaults[i % len(defaults)]
    return palette


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------

def run_pca(ref_psi: dict,
            query_psi: dict,
            variable_jxns: pd.Index,
            ref_names: list,
            query_names: list,
            outprefix: str,
            ref_colors: list | None = None,
            query_colors: list | None = None):
    """
    Joint PCA of all reference + query samples on variable junctions.
    Missing PSI values are imputed with the per-junction mean.

    Reference tissues → small semi-transparent circles.
    Query matrices    → larger diamonds, each in a distinct color.
    """
    frames       = []
    ref_idx      = {}   # tissue name → row indices into coords
    query_idx    = {}   # query name  → row indices into coords
    ref_counts   = {}   # tissue name → sample count (for legend)
    query_counts = {}   # query name  → sample count (for legend)
    cursor       = 0

    for name in ref_names:
        mat = ref_psi[name].reindex(variable_jxns)
        n   = mat.shape[1]
        frames.append(mat)
        ref_idx[name]    = np.arange(cursor, cursor + n)
        ref_counts[name] = n
        cursor += n

    for qname in query_names:
        mat = query_psi[qname].reindex(variable_jxns)
        n   = mat.shape[1]
        frames.append(mat)
        query_idx[qname]    = np.arange(cursor, cursor + n)
        query_counts[qname] = n
        cursor += n

    combined  = pd.concat(frames, axis=1)          # junctions × samples
    row_means = combined.mean(axis=1)
    combined  = combined.apply(lambda c: c.fillna(row_means), axis=0)  # fill before transpose
    combined  = combined.T                          # samples × junctions
    combined  = combined.dropna(axis=1)

    X       = StandardScaler().fit_transform(combined.values)
    pca     = PCA(n_components=min(10, X.shape[1], X.shape[0]))
    coords  = pca.fit_transform(X)
    var_exp = pca.explained_variance_ratio_ * 100

    tissue_pal = _tissue_palette(ref_names, ref_colors)
    query_pal  = _query_palette(query_names, query_colors)

    n_jxns = len(variable_jxns)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle(f"PCA of PSI values (N={n_jxns:,} discriminating junctions)", fontsize=12)

    for ax, (pcx, pcy) in zip(axes, [(0, 1), (0, 2)]):
        if coords.shape[1] <= pcy:
            ax.set_visible(False)
            continue
        for tissue in ref_names:
            idx = ref_idx[tissue]
            ax.scatter(coords[idx, pcx], coords[idx, pcy],
                       color=tissue_pal[tissue], alpha=0.30, s=14,
                       marker="o", linewidths=0, zorder=2)
        for qname in query_names:
            idx = query_idx[qname]
            ax.scatter(coords[idx, pcx], coords[idx, pcy],
                       color=query_pal[qname], alpha=0.90, s=80,
                       marker="D", edgecolors="white", linewidths=0.3,
                       zorder=4)
        ax.set_xlabel(f"PC{pcx+1} ({var_exp[pcx]:.1f}%)", fontsize=9)
        ax.set_ylabel(f"PC{pcy+1} ({var_exp[pcy]:.1f}%)", fontsize=9)
        ax.set_title(f"PC{pcx+1} vs PC{pcy+1}", fontsize=10)
        ax.tick_params(labelsize=8)

    ref_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=tissue_pal[t], markeredgecolor="none",
                   markersize=6, alpha=0.7,
                   label=f"{t} (n={ref_counts[t]})")
        for t in ref_names
    ]
    query_handles = [
        plt.Line2D([0], [0], marker="D", color="w",
                   markerfacecolor=query_pal[q], markeredgecolor="grey",
                   markeredgewidth=0.4, markersize=7,
                   label=f"{q} (n={query_counts[q]})")
        for q in query_names
    ]
    ncol = min(len(ref_handles) + len(query_handles), 6)
    fig.legend(handles=ref_handles + query_handles, loc="lower center", ncol=ncol,
               fontsize=8, bbox_to_anchor=(0.5, -0.02), frameon=False)
    fig.tight_layout(rect=[0, 0.06, 1, 1])

    path = f"{outprefix}_PCA.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")

    return pd.DataFrame(
        coords[:, :min(5, coords.shape[1])],
        index=combined.index,
        columns=[f"PC{i+1}" for i in range(min(5, coords.shape[1]))]
    ), var_exp


# ---------------------------------------------------------------------------
# Mean |ΔPSI| distance scoring
# ---------------------------------------------------------------------------

def compute_distance_scores(ref_avg_psi: pd.DataFrame,
                             query_psi_all: pd.DataFrame,
                             variable_jxns: pd.Index) -> pd.DataFrame:
    """
    For each query sample, compute mean |PSI_sample − avg_PSI_tissue| across
    all variable junctions with non-NaN PSI in that sample.
    Lower score = closer match.
    """
    jxns = variable_jxns.intersection(ref_avg_psi.index)
    ref  = ref_avg_psi.loc[jxns]
    q    = query_psi_all.reindex(jxns)

    scores = {}
    for sample in q.columns:
        s_psi = q[sample]
        valid = s_psi.notna()
        row   = {}
        for tissue in ref.columns:
            t_psi = ref[tissue]
            both  = valid & t_psi.notna()
            row[tissue] = (s_psi[both] - t_psi[both]).abs().mean() if both.sum() > 0 else np.nan
        scores[sample] = row

    return pd.DataFrame(scores).T   # samples × tissues


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_distance_heatmap(score_df: pd.DataFrame,
                          sample_to_query: dict,
                          query_names: list,
                          outprefix: str,
                          query_colors: list | None = None,
                          n_variable_jxns: int = 0):
    """
    Heatmap: rows = query samples (grouped + clustered by source matrix),
    cols = reference tissues.  A color sidebar identifies each group.
    """
    n_samples, n_tissues = score_df.shape
    query_pal = _query_palette(query_names, query_colors)

    # Cluster rows within each query group
    ordered_rows     = []
    group_boundaries = []
    for qname in query_names:
        members = [s for s in score_df.index if sample_to_query.get(s) == qname]
        if not members:
            continue
        sub      = score_df.loc[members].values.astype(float)
        sub_fill = np.where(np.isnan(sub), np.nanmean(sub, axis=0), sub)
        if len(members) > 2 and np.all(np.isfinite(sub_fill)):
            order = leaves_list(
                linkage(pdist(sub_fill, metric="euclidean"), method="average")
            )
        else:
            order = np.arange(len(members))
        ordered_rows.extend([members[i] for i in order])
        group_boundaries.append(len(ordered_rows))

    df_plot = score_df.loc[ordered_rows]

    fig_h = max(4, min(n_samples * 0.22 + 1.5, 40))
    fig_w = max(5, n_tissues * 1.2 + 3)

    fig, (ax_side, ax_heat) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [0.03, 1], "wspace": 0.01}
    )

    # Sidebar color strip
    side_arr   = np.zeros((n_samples, 1, 3))
    row_cursor = 0
    for qname in query_names:
        members = [s for s in ordered_rows if sample_to_query.get(s) == qname]
        n = len(members)
        if n == 0:
            continue
        side_arr[row_cursor:row_cursor + n, 0, :] = matplotlib.colors.to_rgb(query_pal[qname])
        row_cursor += n

    ax_side.imshow(side_arr, aspect="auto", interpolation="nearest")
    ax_side.set_xticks([])
    ax_side.set_yticks(range(n_samples))
    ax_side.set_yticklabels(df_plot.index, fontsize=max(5, min(8, 200 // n_samples)))
    ax_side.yaxis.set_tick_params(length=0)

    sns.heatmap(
        df_plot, ax=ax_heat,
        cmap="YlOrRd_r",
        linewidths=0.3, linecolor="#a0a0a0",
        annot=(n_samples <= 60), fmt=".2f", annot_kws={"size": 7},
        cbar_kws={"label": "Mean |ΔPSI|", "shrink": 0.6},
        yticklabels=False,
    )

    # White separator lines between groups
    for boundary in group_boundaries[:-1]:
        ax_heat.axhline(boundary, color="white", linewidth=2.5, zorder=5)
        ax_side.axhline(boundary - 0.5, color="white", linewidth=2.5, zorder=5)

    ax_heat.set_title(
        f"Tissue identity scores (mean |ΔPSI| for N={n_variable_jxns:,} discriminating junctions)",
        fontsize=11, pad=8,
    )
    ax_heat.set_xlabel("Reference tissue", fontsize=10)
    ax_heat.set_ylabel("")
    ax_heat.tick_params(axis="x", labelsize=9, rotation=30)

    # Legend placed on the right-hand side of the figure, below the colorbar.
    # We use the heatmap axes' inset colorbar axes as an anchor by positioning
    # relative to ax_heat in figure-fraction coordinates.
    legend_handles = [
        mpatches.Patch(color=query_pal[q], label=q)
        for q in query_names
        if any(sample_to_query.get(s) == q for s in ordered_rows)
    ]
    ax_heat.legend(
        handles=legend_handles,
        title="Query matrix",
        loc="upper left",
        bbox_to_anchor=(1.18, 0.44),   # right of ax_heat, below colorbar
        bbox_transform=ax_heat.transAxes,
        fontsize=8,
        title_fontsize=8,
        frameon=True,
        framealpha=0.85,
        borderpad=0.8,
    )

    fig.tight_layout()
    path = f"{outprefix}_distance_heatmap.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.outprefix) or ".", exist_ok=True)

    # ------------------------------------------------------------------
    # 0. Build BED lookup once — used as a filter after every matrix load
    # ------------------------------------------------------------------
    bed = read_bed(args.bed)
    print(f"Loaded {len(bed):,} regions from BED file\n")
    bed_dict = build_bed_dict(bed)
    del bed

    # ------------------------------------------------------------------
    # 1. Reference matrices: read → PSI → coverage filter → BED filter
    # ------------------------------------------------------------------
    ref_psi = {}
    for path, name in zip(args.matrix_refs, args.ref_names):
        ref_psi[name] = load_ref_matrix(
            path, name, bed_dict,
            min_coverage=args.min_coverage,
            min_ref_samples=args.min_ref_samples,
        )

    # ------------------------------------------------------------------
    # 2. Common junctions across all reference matrices
    # ------------------------------------------------------------------
    common_jxns = None
    for name in args.ref_names:
        idx = ref_psi[name].index
        common_jxns = idx if common_jxns is None else common_jxns.intersection(idx)
    print(f"\n{len(common_jxns):,} junctions present in all reference matrices")

    # ------------------------------------------------------------------
    # 3. Average PSI per reference tissue (over common junctions)
    # ------------------------------------------------------------------
    ref_avg = pd.DataFrame(
        {name: ref_psi[name].reindex(common_jxns).mean(axis=1)
         for name in args.ref_names}
    )   # junctions × tissues

    # ------------------------------------------------------------------
    # 4. Select tissue-discriminating junctions
    # ------------------------------------------------------------------
    print("\nSelecting tissue-discriminating junctions...")
    variable_jxns = select_variable_junctions(
        ref_avg,
        n_override=args.n_variable_junctions,
        outprefix=args.outprefix,
    )

    # Reference PSI matrices are only needed (downstream) on variable_jxns
    # for the PCA scatter plot. Shrink them now to drop the bulk of the
    # BED-filtered junctions (e.g. thousands -> ~100 rows per tissue).
    for name in args.ref_names:
        ref_psi[name] = ref_psi[name].reindex(variable_jxns)
    gc.collect()

    # ------------------------------------------------------------------
    # 5. Query matrices: read → PSI → coverage filter → BED filter
    # ------------------------------------------------------------------
    query_psi       = {}   # qname → PSI DataFrame
    sample_to_query = {}   # sample column name → query matrix name
    seen_columns    = set()

    for path, qname in zip(args.matrix_query, args.query_names):
        psi, s2q = load_query_matrix(
            path, qname, bed_dict,
            min_coverage=args.min_coverage,
            seen_columns=seen_columns,
            common_jxns=common_jxns,
            variable_jxns=variable_jxns,
        )
        query_psi[qname] = psi
        sample_to_query.update(s2q)

    # ------------------------------------------------------------------
    # 6. PCA
    # ------------------------------------------------------------------
    print("\nRunning PCA...")
    pca_coords, var_exp = run_pca(
        ref_psi, query_psi, variable_jxns,
        args.ref_names, args.query_names, args.outprefix,
        ref_colors=args.ref_colors,
        query_colors=args.query_colors,
    )

    # ------------------------------------------------------------------
    # 7. Mean |ΔPSI| distance scores
    # ------------------------------------------------------------------
    print("\nComputing mean |deltaPSI| distance scores...")
    query_psi_all = pd.concat(list(query_psi.values()), axis=1)
    ref_avg_var   = ref_avg.reindex(variable_jxns)
    score_df      = compute_distance_scores(ref_avg_var, query_psi_all, variable_jxns)

    # ------------------------------------------------------------------
    # 8. Plots
    # ------------------------------------------------------------------
    print("\nGenerating heatmap...")
    plot_distance_heatmap(
        score_df, sample_to_query, args.query_names, args.outprefix,
        query_colors=args.query_colors,
        n_variable_jxns=len(variable_jxns),
    )

    print("\nDone.")


if __name__ == "__main__":
    main()