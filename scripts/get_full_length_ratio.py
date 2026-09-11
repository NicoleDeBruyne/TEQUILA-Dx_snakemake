#!/usr/bin/env python3
"""
scripts/get_full_length_ratio.py
For each (gene, sample) pair, computes the "full-length ratio" (FLR) of each
read overlapping that gene in the sample's BAM: the fraction of the gene's
canonical transcript's exonic reference positions that the read's alignment
actually spans via a non-N (non-intron-skip) CIGAR operation (M/D/=/X).
avgFLR for a sample x gene pair is the mean FLR across every primary,
mapped read overlapping the transcript's span. A read that fully spans
every annotated exon of the canonical transcript (no dropped exons, no
truncation) scores close to 1; a read covering only part of the transcript
body (a truncated cDNA/library artifact, a partially-covered amplicon, an
alternate/partial isoform, etc.) scores lower.

Reads directly from each sample's own BAM (SAMPLES[s]["bam"], the same file
_8A/_8B/quantify_gene_expression.py use) with one open file handle per
sample, fetching each gene's canonical-transcript span with pysam's
region-indexed fetch() -- NOT phase_reads.py's per-gene bulk BAM files.
Those per-gene files are themselves nothing more than `samtools view
<region> <this same original BAM>` (see phase_reads.py's
filter_bam_by_region()), so reading them here would mean opening one small
file per (sample, gene) pair -- tens of thousands of file opens for a large
cohort x panel -- for data that's one indexed fetch() away in a file that's
opened once per sample anyway. This also means _8D no longer needs to wait
on phase_reads.py at all.

The "canonical transcript" per gene is picked from --gtf: the transcript
tagged "Ensembl_canonical" (GENCODE/Ensembl convention) if one exists for
that gene, otherwise the transcript with the largest total exonic length.
Its genomic span (min exon start to max exon end) is the fetch() region.

Invoked by rules/8_cohort_qc.smk (_8D_get_full_length_ratio).

Genes are split into two heatmaps by their median (across samples) read
count: genes with too few supporting reads have noisy/meaningless avgFLR,
so grouping them separately (rather than interleaving blank-looking low-
confidence cells throughout one big heatmap) keeps the well-supported genes
readable and makes clear which genes are low-confidence and why. Both live
as stacked subplots in one figure/PDF (high-read-count on top, low-read-
count below) so they're easy to compare side by side. Each panel carries
its own per-gene median-read-count bar (log-scaled, right), a per-sample
avgFLR boxplot (this panel's own genes, one jittered point per gene per
sample), and a per-sample total-read-count bar (log-scaled, bottom, same
value in both panels since it's summed across every gene, not just the
ones shown in that panel). The boxplot points and the total-reads bars are
both colored by sample_type (using the same colors as
_8C2_validate_sample_types -- resolved via the Snakefile's own
sample_type_color() and passed through --mapping-file's color column, not
recomputed here), with a shared legend on the figure -- no separate color
strip. Column width per sample shrinks (and the overall figure width is
capped) as sample count grows, so large cohorts (e.g. ~200 samples) stay a
print-friendly width instead of many feet wide.

Outputs:
  --outprefix + "_matrix.tsv"              -- genes x samples avgFLR matrix
                                              (every gene, unsplit).
  --outprefix + "_read_counts_matrix.tsv"  -- genes x samples read-count
                                              matrix backing the avgFLR
                                              values and the bar charts.
  --outprefix + "_heatmap.pdf"             -- one figure, two stacked
                                              subplots: avgFLR heatmap +
                                              per-sample avgFLR boxplot +
                                              read-count bar, for genes
                                              with median read count >=
                                              --min-reads (default 100) on
                                              top, and < --min-reads below.
"""

import argparse
import gzip
import os
import concurrent.futures

import pandas as pd
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches
import seaborn as sns
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-gene, per-sample full-length transcript coverage ratio (FLR) from each sample's BAM")
    parser.add_argument(
        "--mapping-file",
        required=True,
        help="TSV file without a header and with the columns: name, bam, sample_type, color -- same "
             "convention as scripts/quantify_gene_expression.py's mapping file, plus sample_type and "
             "a pre-resolved hex color (same colors _8C2_validate_sample_types uses, resolved by the "
             "calling rule via the Snakefile's sample_type_color() -- not recomputed here) for the "
             "heatmap's per-sample boxplot points and read-count bars.")
    parser.add_argument("--bed", required=True,
        help="BED file for the panel shared by every sample in --mapping-file. One row per gene "
             "(chrom, start, end, gene, ...); column 4 is the gene symbol. Only used to get the gene "
             "list -- the actual region fetched per gene is its canonical transcript's own span from --gtf.")
    parser.add_argument("--gtf", required=True, help="Path to the reference annotation GTF (plain or .gz).")
    parser.add_argument("--outprefix", required=True, help="Prefix for all output files.")
    parser.add_argument("--title", default="Full-length transcript coverage ratio (FLR)", help="Plot title prefix.")
    parser.add_argument("--min-reads", type=int, default=100,
        help="Genes whose median (across samples) read count is below this are split into their own "
             "low-confidence heatmap instead of the main one. Default: 100")
    parser.add_argument("--threads", type=int, default=1, help="Number of samples to process in parallel. Default: 1")
    return parser.parse_args()


def load_gene_list(bed):
    """Gene symbols from BED column 4, in file order, de-duplicated -- same
    convention as scripts/quantify_gene_expression.py's load_gene_regions()
    and scripts/phase_reads.py's extract_gene_regions()."""
    genes = []
    seen = set()
    with open(bed) as b:
        for line in b:
            if not line.strip() or line.startswith('#'):
                continue
            gene = line.strip().split('\t')[3]
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return genes


# Reference-consuming, non-N CIGAR ops: M=0, D=2, ==7, X=8. I=1/S=4/H=5/P=6
# don't consume reference; N=3 (intron skip) is explicitly excluded.
_REF_CONSUMING_NON_N = {0, 2, 7, 8}


def _open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _gtf_attr(attr_str, key):
    # GTF attribute fields look like: gene_name "FOO"; transcript_id "ENST...";
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) != 2:
            continue
        k, v = parts
        if k == key:
            return v.strip().strip('"')
    return None


def load_canonical_transcripts(gtf_path, genes):
    """Returns {gene: (chrom, span_start, span_end, [(exon_start, exon_end), ...], exonic_length)}
    for the canonical transcript of every gene in `genes` found in the GTF.
    All coordinates are 0-based half-open, matching pysam's reference
    coordinate convention. (chrom, span_start, span_end) is the transcript's
    overall genomic footprint (min exon start to max exon end) -- the region
    passed to AlignmentFile.fetch()."""
    genes = set(genes)

    # transcript_id -> dict(gene, chrom, strand, canonical, exons=[(start,end),...])
    transcripts = {}

    with _open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields[:9]
            if feature not in ("transcript", "exon"):
                continue
            gene_name = _gtf_attr(attrs, "gene_name") or _gtf_attr(attrs, "gene_id")
            if gene_name not in genes:
                continue
            transcript_id = _gtf_attr(attrs, "transcript_id")
            if transcript_id is None:
                continue

            if feature == "transcript":
                is_canonical = 'tag "Ensembl_canonical"' in attrs
                t = transcripts.setdefault(transcript_id, {
                    "gene": gene_name, "chrom": chrom, "strand": strand,
                    "canonical": False, "exons": [],
                })
                t["canonical"] = t["canonical"] or is_canonical
            else:  # exon
                t = transcripts.setdefault(transcript_id, {
                    "gene": gene_name, "chrom": chrom, "strand": strand,
                    "canonical": False, "exons": [],
                })
                # GTF is 1-based inclusive -> convert to 0-based half-open.
                t["exons"].append((int(start) - 1, int(end)))

    # Pick, per gene: the Ensembl_canonical-tagged transcript if any, else
    # the transcript with the largest total exonic length.
    by_gene = {}
    for tid, t in transcripts.items():
        if not t["exons"]:
            continue
        by_gene.setdefault(t["gene"], []).append(t)

    result = {}
    for gene, cands in by_gene.items():
        canonical_cands = [c for c in cands if c["canonical"]]
        pool = canonical_cands if canonical_cands else cands
        best = max(pool, key=lambda c: sum(e - s for s, e in c["exons"]))
        exons = sorted(best["exons"])
        # Merge any overlapping/adjacent exon records (defensive -- GENCODE
        # exons are normally already disjoint per transcript).
        merged = []
        for s, e in exons:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        exonic_length = sum(e - s for s, e in merged)
        span_start = merged[0][0]
        span_end = merged[-1][1]
        result[gene] = (best["chrom"], span_start, span_end, merged, exonic_length)

    return result


def _read_ref_blocks(read):
    """Reference-coordinate (0-based half-open) blocks covered by this
    read's non-N, reference-consuming CIGAR ops (M/D/=/X)."""
    blocks = []
    pos = read.reference_start
    block_start = None
    for op, length in read.cigartuples or []:
        if op in _REF_CONSUMING_NON_N:
            if block_start is None:
                block_start = pos
            pos += length
        else:
            if block_start is not None:
                blocks.append((block_start, pos))
                block_start = None
            if op == 3:  # N: intron skip, still consumes reference
                pos += length
    if block_start is not None:
        blocks.append((block_start, pos))
    return blocks


def _overlap_length(blocks_a, blocks_b):
    """Total overlap length between two lists of disjoint, sorted (start,end)
    interval tuples."""
    i = j = 0
    total = 0
    while i < len(blocks_a) and j < len(blocks_b):
        a_s, a_e = blocks_a[i]
        b_s, b_e = blocks_b[j]
        lo, hi = max(a_s, b_s), min(a_e, b_e)
        if lo < hi:
            total += hi - lo
        if a_e < b_e:
            i += 1
        else:
            j += 1
    return total


def compute_avg_flr(bam_handle, chrom, span_start, span_end, exons, exonic_length):
    """Returns (avg_flr, n_reads), fetching from an already-open
    AlignmentFile via an indexed region query -- no per-gene file open.
    avg_flr is NaN when there's no usable data; n_reads is always a real
    count (0 when there's nothing)."""
    if exonic_length <= 0:
        return np.nan, 0
    ratios = []
    for read in bam_handle.fetch(chrom, span_start, span_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        read_blocks = _read_ref_blocks(read)
        if not read_blocks:
            continue
        overlap = _overlap_length(read_blocks, exons)
        ratios.append(overlap / exonic_length)
    if not ratios:
        return np.nan, 0
    return float(np.mean(ratios)), len(ratios)


def _process_sample(sample, bam_path, transcripts):
    """Opens `bam_path` once and fetches every gene's region against that
    single handle -- one file open per sample, not one per (sample, gene)."""
    flr_out = {}
    count_out = {}
    if not bam_path or not os.path.exists(bam_path):
        for gene in transcripts:
            flr_out[gene] = np.nan
            count_out[gene] = 0
        return sample, flr_out, count_out

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for gene, (chrom, span_start, span_end, exons, exonic_length) in transcripts.items():
            avg_flr, n_reads = compute_avg_flr(bam, chrom, span_start, span_end, exons, exonic_length)
            flr_out[gene] = avg_flr
            count_out[gene] = n_reads
    return sample, flr_out, count_out


def _panel_height(flr_df):
    """Height (inches) _draw_heatmap_panel will need for this gene set --
    computed upfront (pure function of gene count) so the outer figure's
    gridspec can be given correct height_ratios before any subplot exists,
    rather than resizing after the fact (which would leave the actual
    per-panel gridspec proportions wrong)."""
    n_genes = len(flr_df.index)
    box_h = 1.6   # per-sample avgFLR boxplot row
    bar_h = 1.2   # per-sample total-read-count bar row
    if flr_df.empty or n_genes == 0:
        return 1.5
    row_height = 0.22 if n_genes <= 60 else max(0.05, 0.22 * (60 / n_genes))
    heat_h = min(max(4, row_height * n_genes), 36)
    return heat_h + box_h + bar_h + 1.2


def _draw_heatmap_panel(fig, outer_spec, flr_df, gene_median_reads, sample_total_reads,
                         sample_colors, panel_title, heat_w):
    """Draws one avgFLR heatmap, a per-sample avgFLR boxplot (this panel's
    own genes, points colored by sample_type), and a per-sample total-
    read-count bar (also colored by sample_type) into `outer_spec` (a
    SubplotSpec carved out of `fig`'s outer gridspec). `outer_spec`'s own
    height is assumed to already match _panel_height(flr_df) -- see
    make_combined_heatmap, which sizes the outer gridspec's row before
    calling this. `sample_total_reads` and `sample_colors` are expected to
    already cover every sample (not just the ones in `flr_df`'s columns),
    so they read the same across every panel."""
    n_genes = len(flr_df.index)
    bar_w = 1.6    # right (gene median reads) bar width, inches
    box_h = 1.6    # avgFLR boxplot row height, inches
    bar_h = 1.2    # total-reads bar row height, inches

    if flr_df.empty:
        # Still draw something -- this panel is part of a declared
        # Snakemake output, so the figure must exist even when a run's gene
        # panel happens to have nothing in this read-count bucket.
        ax = fig.add_subplot(outer_spec)
        ax.axis("off")
        ax.text(0.5, 0.5, "No genes in this group.", ha="center", va="center")
        ax.set_title(panel_title)
        return

    gene_order = flr_df.mean(axis=1).sort_values(ascending=False).index
    sample_order = flr_df.mean(axis=0).sort_values(ascending=False).index
    ordered = flr_df.loc[gene_order, sample_order]
    n_samples = len(sample_order)
    gmed = gene_median_reads.reindex(gene_order)
    stot = sample_total_reads.reindex(sample_order)
    scolors = [sample_colors.get(s, "#555555") for s in sample_order]

    row_height = 0.22 if n_genes <= 60 else max(0.05, 0.22 * (60 / n_genes))
    heat_h = min(max(4, row_height * n_genes), 36)

    # Only individually label genes when there's room to read them; past
    # that, still-visible sort order (best/worst avgFLR) carries the
    # information without needing every gene name printed.
    show_gene_labels = (heat_h / max(n_genes, 1)) * 72 >= 5  # ~5pt/row floor
    gene_fontsize = min(8, max(4, (heat_h * 72) / max(n_genes, 1) - 1))

    # Same idea on the sample axis: past a certain density, individual
    # sample-name labels just overlap into noise. Per-sample color (shared
    # by the boxplot points and the reads bar below) still shows grouping,
    # and exact identities are always in the companion _matrix.tsv.
    show_sample_labels = (heat_w / max(n_samples, 1)) * 72 >= 5
    sample_fontsize = min(7, max(3, (heat_w * 72) / max(n_samples, 1) - 1))

    inner = outer_spec.subgridspec(
        3, 3,
        width_ratios=[heat_w, bar_w, 0.25],
        height_ratios=[heat_h, box_h, bar_h],
        wspace=0.08, hspace=0.08,
    )
    ax_heat = fig.add_subplot(inner[0, 0])
    ax_gbar = fig.add_subplot(inner[0, 1], sharey=ax_heat)
    ax_cbar = fig.add_subplot(inner[0, 2])
    ax_box  = fig.add_subplot(inner[1, 0], sharex=ax_heat)
    ax_sbar = fig.add_subplot(inner[2, 0], sharex=ax_heat)

    sns.heatmap(ordered, ax=ax_heat, cmap="viridis", vmin=0, vmax=1,
                cbar=True, cbar_ax=ax_cbar, cbar_kws={"label": "avgFLR"},
                linewidths=0, yticklabels=show_gene_labels)
    if show_gene_labels:
        ax_heat.tick_params(axis="y", labelsize=gene_fontsize)
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("gene")
    ax_heat.set_title(panel_title, fontsize=10)
    ax_heat.tick_params(axis="x", labelbottom=False, bottom=False)

    # Per-gene median read count, right, log-scaled. Clipped to >=1 so
    # zero-read genes (common in the low-read-count panel) still show up as
    # a visible sliver rather than vanishing on a log axis.
    y = np.arange(n_genes) + 0.5
    ax_gbar.barh(y, gmed.clip(lower=1).to_numpy(), height=0.8, color="#555555")
    ax_gbar.set_xscale("log")
    ax_gbar.set_xlabel("median\nreads", fontsize=7)
    ax_gbar.tick_params(axis="y", left=False, labelleft=False)
    ax_gbar.set_ylim(ax_heat.get_ylim())

    # Per-sample avgFLR distribution across this panel's own genes -- box
    # outline in black/gray, one jittered point per gene colored by that
    # sample's sample_type. Column positions (x = i+0.5) match the heatmap
    # above it exactly, same convention as every other row in this panel.
    x = np.arange(n_samples) + 0.5
    box_data = [ordered[s].dropna().to_numpy() for s in sample_order]
    _black = dict(color="black", linewidth=0.7)
    bp = ax_box.boxplot(box_data, positions=x, widths=0.6, showfliers=False,
                         boxprops=_black, whiskerprops=_black, capprops=_black, medianprops=_black)
    rng = np.random.RandomState(0)
    for xi, vals, c in zip(x, box_data, scolors):
        if len(vals) == 0:
            continue
        jitter = (rng.rand(len(vals)) - 0.5) * 0.35
        ax_box.scatter(np.full(len(vals), xi) + jitter, vals, color=c, s=6, zorder=3, edgecolors="none")
    ax_box.set_ylim(-0.02, 1.02)
    ax_box.set_ylabel("avgFLR", fontsize=7)
    ax_box.tick_params(axis="x", labelbottom=False, bottom=False)

    # Per-sample total read count (summed across every gene, not just this
    # panel's subset), bottom, log-scaled, colored the same as the boxplot
    # points above (i.e. by sample_type).
    ax_sbar.bar(x, stot.clip(lower=1).to_numpy(), width=0.8, color=scolors)
    ax_sbar.set_yscale("log")
    ax_sbar.set_ylabel("total\nreads", fontsize=7)
    if show_sample_labels:
        ax_sbar.set_xticks(x)
        ax_sbar.set_xticklabels(sample_order, rotation=90, fontsize=sample_fontsize)
    else:
        ax_sbar.tick_params(axis="x", labelbottom=False, bottom=False)


def make_combined_heatmap(flr_high, flr_low, gene_median_reads, sample_total_reads,
                           sample_type_map, sample_colors, out_pdf, title_high, title_low,
                           suptitle, n_samples):
    """One figure, two stacked subplots (high-read-count panel on top,
    low-read-count panel below), each built by _draw_heatmap_panel. Panel
    heights are computed upfront so the outer gridspec's row proportions
    (and therefore the overall figure size) are correct on the first and
    only draw -- no resize-after-the-fact, which would leave the actual
    row split wrong relative to each panel's true gene count.

    Column width shrinks as sample count grows (same idea as the row-height
    shrink for genes) and the total width is capped, so a ~200-sample
    cohort gets a print-friendly figure instead of a many-foot-wide PDF."""
    col_width = 0.28 if n_samples <= 50 else max(0.05, 0.28 * (50 / n_samples))
    heat_w = min(max(5, col_width * n_samples + 1), 20)
    full_w = heat_w + 1.6 + 1.0   # heat_w + bar_w + margin, matches _draw_heatmap_panel's inner layout

    h_high = _panel_height(flr_high)
    h_low = _panel_height(flr_low)

    fig = plt.figure(figsize=(full_w, h_high + h_low + 1.0))
    outer = fig.add_gridspec(2, 1, height_ratios=[h_high, h_low], hspace=0.4)

    _draw_heatmap_panel(fig, outer[0, 0], flr_high, gene_median_reads, sample_total_reads,
                         sample_colors, title_high, heat_w)
    _draw_heatmap_panel(fig, outer[1, 0], flr_low, gene_median_reads, sample_total_reads,
                         sample_colors, title_low, heat_w)

    # One legend for the whole figure -- sample_type -> color is the same
    # in both panels (and matches rules/8_cohort_qc.smk's
    # _8C2_validate_sample_types, which resolves colors via the Snakefile's
    # own sample_type_color()).
    seen = {}
    for s, t in sample_type_map.items():
        seen.setdefault(t, sample_colors.get(s, "#555555"))
    legend_handles = [matplotlib.patches.Patch(color=c, label=t) for t, c in seen.items()]
    fig.legend(handles=legend_handles, loc="upper right", title="sample_type",
               fontsize=7, title_fontsize=7, frameon=False)

    fig.suptitle(suptitle)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()

    df = pd.read_csv(args.mapping_file, sep="\t", header=None, names=["sample", "bam", "sample_type", "color"])
    samples = sorted(df["sample"].dropna().unique())
    sample_to_bam = dict(zip(df["sample"], df["bam"]))
    sample_to_type = dict(zip(df["sample"], df["sample_type"]))
    sample_to_color = dict(zip(df["sample"], df["color"]))

    genes = load_gene_list(args.bed)
    if not genes:
        raise ValueError("No genes found in BED file: " + args.bed)

    print("Loading canonical transcripts for " + str(len(genes)) + " gene(s) from " + args.gtf + "...")
    transcripts = load_canonical_transcripts(args.gtf, genes)
    missing = [g for g in genes if g not in transcripts]
    if missing:
        print("[WARNING] No canonical transcript found in GTF for " + str(len(missing)) +
              " gene(s) (will be NaN in output): " + ", ".join(missing[:20]) +
              (" ..." if len(missing) > 20 else ""))

    print("Computing avgFLR across " + str(len(samples)) + " sample(s) x " + str(len(genes)) +
          " gene(s) using " + str(args.threads) + " thread(s)...")

    flr_results = {}
    count_results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(_process_sample, s, sample_to_bam.get(s), transcripts): s for s in samples}
        for f in concurrent.futures.as_completed(futures):
            sample, flr_vals, count_vals = f.result()
            flr_results[sample] = flr_vals
            count_results[sample] = count_vals
            print("Finished: " + sample)

    matrix = pd.DataFrame(flr_results).reindex(index=genes, columns=samples)
    matrix.index.name = "gene"

    count_matrix = pd.DataFrame(count_results).reindex(index=genes, columns=samples).fillna(0).astype(int)
    count_matrix.index.name = "gene"

    out_matrix = args.outprefix + "_matrix.tsv"
    matrix.to_csv(out_matrix, sep="\t")
    print("Saved avgFLR matrix: " + out_matrix)

    out_counts = args.outprefix + "_read_counts_matrix.tsv"
    count_matrix.to_csv(out_counts, sep="\t")
    print("Saved read-count matrix: " + out_counts)

    # Per-gene median read count (across samples) decides which heatmap a
    # gene lands in; per-sample total read count (across every gene, not
    # just the genes in one heatmap) backs the bottom bar in both.
    gene_median_reads = count_matrix.median(axis=1)
    sample_total_reads = count_matrix.sum(axis=0)

    high_genes = gene_median_reads[gene_median_reads >= args.min_reads].index
    low_genes = gene_median_reads.index.difference(high_genes)

    out_heatmap = args.outprefix + "_heatmap.pdf"
    make_combined_heatmap(
        matrix.loc[high_genes], matrix.loc[low_genes], gene_median_reads, sample_total_reads,
        sample_to_type, sample_to_color,
        out_heatmap,
        title_high="median reads >= " + str(args.min_reads) + " (" + str(len(high_genes)) + " gene(s))",
        title_low="median reads < " + str(args.min_reads) + " (" + str(len(low_genes)) + " gene(s))",
        suptitle=args.title + " -- heatmap",
        n_samples=len(samples),
    )
    print("Saved avgFLR heatmap: " + out_heatmap +
          " (" + str(len(high_genes)) + " high-read-count, " + str(len(low_genes)) + " low-read-count gene(s))")


if __name__ == "__main__":
    main()
