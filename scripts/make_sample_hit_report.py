"""
make_sample_hit_report.py

One PDF per sample (US Letter, landscape), in this order: an overview page; IGV-style views of every
candidate gene hit in merged_all_hits.tsv (default filter set, hit == TRUE, tier <= --max-tier,
ordered by `ranking`); then the gene-level sections -- gene expression (MOTR z-score), canonical-
transcript coverage, junction testability, full-length ratio.

QC pages. Every panel compares the sample with each sample type on the panel (the
sample is left out of its own type). Colours: the sample in a darker shade of its type's colour,
each sample type in its own colour; gene grids run red (0) - white (0.5) - blue (1), missing values
hatched; all 0-1 metric axes run 0-1.
  Overview page: read counts (stacked on-/off-target bars for the sample and each type's median
  mapped reads x median on-target rate; _9A), read length (on-/off-target boxes for the sample and,
  per type, a box built from the median Q1 / median / Q3 / min / max of its samples, whiskers
  1.5 IQR clipped to min / max; _9B), sample-type PCA (_9C2; the sample as a star), and a grid
  of the sample's candidate genes by combination of hit categories (variant, ASE, GTEx or cohort
  outlier junction; upset-style dots under the columns) x tier, with counts in the cells and
  marginal totals.
  After the gene pages: gene expression (MOTR z-score of the sample type, by splice-site
  assignment), then one section each for canonical-transcript coverage (fraction of bases at
  depth >= merge_min_dp_snv; _9D), junction testability (n_canonical_jxns_tested /
  n_canonical_jxns, from merged_all_hits.tsv) and full-length ratio (_9E): a summary strip (ranked-gene curves with a 25th-75th
  percentile ribbon and median per sample type and the sample as a bold line; a per-gene box plot
  of the sample and each type's median; a per-sample-mean box plot by type with the sample as a
  star) and a grid of every gene (name, a cell for the sample, one per type's median), sorted high
  to low by the sample's value, then the median of every other sample, then name. The label font
  shrinks with gene count (6.5 pt down to 4 pt); beyond that the grid continues on further pages.
  Gene names are text, so a PDF viewer's search finds them.
A panel whose input is not given is left out.

Every gene starts on a new page and is laid out as rows, 3 rows per page:
  1. Bulk alignments over the whole gene body (the gene's region in the BED), in the left 75%
     of the width; to its right, the gene's expression (by-assignment log2(MOTR + 1)) as a box
     per sample type with the sample as a star, and its MOTR z-score within its own type.
  2. hap1 and hap2 alignments over the gene body (same 75% width), if the gene was phased and
     has ASE or a hap1 / hap2 outlier junction.
  3. Variant zooms (15 bp either side of each candidate variant), up to 3 panels per row, bulk.
  4. Metric distributions (PSI_approx, PSI, IR, IPA) for each junction x metric that fired,
     across every sample on the panel split by sample type, up to 4 panels per row.
  5. Outlier splicing events (zooms), one per row, bulk: the event's junctions plus 10% of their span on
     either side. The IGV view (first two columns) has arcs for the event's junctions only; the
     third column is a sashimi plot of the same window (introns scaled by --intron-scale,
     unspliced non-exonic reads dropped for multi-exon genes, every junction with >= --sashimi-min-frac x the
     track's maximum coverage) for this sample and a representative cohort sample. An event called
     in hap1 or hap2 gets its own page with bulk, hap1 and hap2 rows (all three in the sashimi
     too). Title: event label(s) and each junction with its annotation (canonical / annotated /
     unannotated).
     Representative sample (one per gene): among the other samples of this sample's type (or, if
     there are none, of the cohort's most common sample type), leaving out samples with their own
     outlier junction call in the gene, the one whose bulk junction PSIs (cohort junction tables;
     junctions with coverage >= --coverage-threshold) are closest to the candidates' median PSI
     (mean |PSI - median| over at least 3 shared junctions, or all the gene has if fewer). Among
     candidates within 0.02 of the best, the one with the most junction reads.

IGV views (coverage and alignment tracks) show every read. In the sashimi plots of a multi-exon
gene, unspliced reads are left out unless they lie within an annotated exon (any transcript,
+/- 10 bp) -- so reads inside a long exon stay, pre-mRNA / intronic reads go. Sashimi arcs run from the last exon base before the
intron (ss1 - 1) to the first exon base after it (ss2 + 1).

Each alignment panel has a coverage track (all reads; positions where a non-reference base has
allele frequency > 0.2 at depth >= 10 are coloured by base, as in IGV), an alignment track
(primary and supplementary alignments, randomly downsampled to --max-reads reads; every mismatch
coloured by base, deletions >= 10 bp as black lines, insertions >= 10 bp as purple ticks,
smaller indels hidden unless they match a candidate variant -- same length, within 5 bp -- in the
gene view, or the zoomed variant in a variant zoom; splice gaps as thin grey lines), a reference-sequence track on
variant zooms, and the gene's canonical transcript (thick = CDS, thin = UTR), drawn once under the
bottom track of a bulk / hap1 / hap2 stack. Samples are always shown by name, never by alias.

Outlier splicing events come from the hit table: junctions are grouped by the event ID that _9L /
_8C assigned within each source (GTEx tissue or cohort) and phasing (older outputs without IDs fall
back to grouping same-label junctions that share a splice site); calls of the same junction set
from different sources, tissues or phasings are shown as one event.
"""
import argparse
import bisect
import gzip
import io
import itertools
import os
import random
import re
import textwrap
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection, LineCollection, PatchCollection
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle, FancyBboxPatch
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap, to_rgb, to_hex

from hits_io import read_hit_rows
from palette import RED_WHITE_BLUE
import tier_rules

rcParams["pdf.fonttype"] = 42
FS = 1.4                              # font scale for every text element in the report
rcParams["font.size"] = 6 * FS

# ----------------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------------
PAGE_W, PAGE_H = 11.0, 8.5            # US Letter, landscape
MARGIN_L, MARGIN_R = 0.75, 0.25
MARGIN_T, MARGIN_B = 0.70, 0.25
ROWS_PER_PAGE = 3
VARIANT_PANELS_PER_ROW = 3
DIST_PANELS_PER_ROW = 4
PANEL_GAP = 0.18                      # inches between panels in a multi-panel row

BASE_COLORS = {"A": "#00a600", "C": "#0000ff", "G": "#d17105", "T": "#ff0000", "N": "#888888"}
READ_COLOR = "#bdbdbd"
DEL_COLOR = "#000000"
INS_COLOR = "#8a2be2"
GAP_COLOR = "#9e9e9e"
COV_COLOR = "#a6a6a6"
ARC_COLOR = "#b2182b"                 # reset per sample (draw_qc_page) to the sample colour: a darker shade of its type colour
SPAN_COLOR = "#fff2a8"

RASTER_DPI = 300                      # resolution of the few other rasterised artists (coverage fill)
ALIGN_DPI = 1200                      # resolution of the alignment-track images
AF_THRESHOLD = 0.2
AF_MIN_DEPTH = 10
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400   # unmapped, secondary, QC-fail, duplicate (primary + supplementary kept)
MIN_INDEL = 10                             # insertions / deletions shorter than this are not drawn,
                                           # unless they are a candidate variant (see indel_specs)
GENE_IGV_FRAC = 0.75                       # share of the width for the gene-body IGV views (rest: expression box)
INDEL_POS_TOL = 5                          # bp an aligner may shift a variant indel within a repeat

GTEX_LABELS = {"exon_skipping_approx", "exon_inclusion_approx", "alt_5ss_approx", "alt_3ss_approx",
               "complex_approx"}
COHORT_LABEL_METRIC = {
    **{l: "junction_PSI_approx" for l in GTEX_LABELS},
    "exon_skipping": "junction_PSI", "exon_inclusion": "junction_PSI", "alt_5ss": "junction_PSI",
    "alt_3ss": "junction_PSI", "complex": "junction_PSI",
    "5ss_IR": "5ss_IR_ratio", "3ss_IR": "3ss_IR_ratio", "full_IR": "junction_full_IR_ratio",
    "IPA": "junction_IPA_ratio",
}
SINGLE_JUNCTION_LABELS = {"5ss_IR", "3ss_IR", "full_IR", "IPA"}
COHORT_DELTA_SUFFIX = {
    "junction_PSI": "deltaPSI", "junction_PSI_approx": "deltaPSIapprox", "5ss_IR_ratio": "delta5ssIR",
    "3ss_IR_ratio": "delta3ssIR", "junction_full_IR_ratio": "deltaFullIR", "junction_IPA_ratio": "deltaIPA",
}
# 8A per-gene table columns: the rescaled metric value, and the denominator that must reach
# --coverage-threshold for a sample to be plotted
_METRIC_RESCALED_COL = {
    "junction_PSI_approx":    "rescaled_junction_PSI_approx",
    "junction_PSI":           "rescaled_junction_PSI",
    "5ss_IR_ratio":           "rescaled_5ss_IR_ratio",
    "3ss_IR_ratio":           "rescaled_3ss_IR_ratio",
    "junction_full_IR_ratio": "rescaled_junction_full_IR_ratio",
    "junction_IPA_ratio":     "rescaled_junction_IPA_ratio",
}
_METRIC_DENOM_COL = {
    "junction_PSI_approx":    "junction_coverage_approx",
    "junction_PSI":           "junction_coverage",
    "5ss_IR_ratio":           "5ss_coverage",
    "3ss_IR_ratio":           "3ss_coverage",
    "junction_full_IR_ratio": "junction_coverage",
    "junction_IPA_ratio":     "5ss_coverage",
}
METRIC_SHORT = {
    "junction_PSI_approx": "PSI_approx", "junction_PSI": "PSI", "5ss_IR_ratio": "5'ss IR",
    "3ss_IR_ratio": "3'ss IR", "junction_full_IR_ratio": "full IR", "junction_IPA_ratio": "IPA",
}
_MISSING = {"", ".", "nan", "none", "n/a"}


def parse_args():
    p = argparse.ArgumentParser(description="Per-sample final report: a QC page, then IGV-style pages "
                                            "for the sample's candidate genes.")
    p.add_argument("--sample", required=True)
    p.add_argument("--hits-tsv", required=True, help="Bed-level merged_all_hits.tsv (default filter set).")
    p.add_argument("--bam", required=True, help="The sample's bulk BAM.")
    p.add_argument("--gene-bam-mapping", default=None,
                   help="The sample's phased_reads/{sample}_gene_bam_mapping_file.tsv (hap1/hap2 BAMs per gene).")
    p.add_argument("--bed", required=True, help="Panel BED; column 4 = gene, columns 1-3 = gene-body region.")
    p.add_argument("--gtf", required=True)
    p.add_argument("--genome", required=True, help="Indexed reference FASTA.")
    p.add_argument("--cohort-manifests", nargs="*", default=[], metavar="SAMPLE_TYPE:PATH",
                   help="8A gene manifests ('sample_type:path'), for the metric distribution panels.")
    p.add_argument("--samples", nargs="*", default=[], help="All samples on the panel (same order as --sample-types).")
    p.add_argument("--sample-types", nargs="*", default=[])
    p.add_argument("--colors", nargs="*", default=[],
                   help="Per-sample hex colour (same order as --samples), from the sample type.")
    p.add_argument("--coverage-threshold", type=int, default=50,
                   help="Minimum metric denominator for a sample to appear in a distribution panel.")
    p.add_argument("--top-n", type=int, default=0,
                   help="Number of top-ranked hit genes to show; 0 (default) shows every hit.")
    p.add_argument("--max-tier", type=int, default=5,
                   help="Only genes with a tier <= this are shown (default 5, i.e. tier 6 is left out).")
    p.add_argument("--max-reads", type=int, default=500, help="Reads drawn per alignment track.")
    p.add_argument("--variant-flank", type=int, default=15, help="Bases shown either side of a variant.")
    p.add_argument("--junction-flank", type=float, default=0.10,
                   help="Fraction of an event's span added on either side.")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--sample-bams", nargs="*", default=[],
                   help="Bulk BAM of every sample on the panel (same order as --samples), for the "
                        "representative sample in the sashimi plots.")
    p.add_argument("--amalgam-matrices", nargs="*", default=[], metavar="SAMPLE_TYPE:PATH",
                   help="Accepted for compatibility; not used (the representative sample is chosen from "
                        "junction PSI in the cohort junction tables).")
    p.add_argument("--sashimi-min-frac", type=float, default=0.10,
                   help="Sashimi: show a junction in a track if its read count is >= this fraction of "
                        "the track's maximum coverage (default 0.10).")
    p.add_argument("--intron-scale", type=float, default=0.25,
                   help="Sashimi: width factor for non-exonic bases (canonical transcript), default 0.25.")
    p.add_argument("--on-target-tsv", default=None,
                   help="QC page: _9A {bed_id}_on_target_rates.tsv (sample, mapped, on_target_rate, ...).")
    p.add_argument("--read-attributes-tsv", default=None,
                   help="QC page: _9B {bed_id}_read_attributes.tsv (sample, targetType, min, q1, median, q3, max).")
    p.add_argument("--flr-matrix", default=None,
                   help="QC page: _9E {bed_id}_full_length_ratio_matrix.tsv (gene x sample avgFLR).")
    p.add_argument("--coverage-matrix", default=None,
                   help="QC: _9D {bed_id}_breadth_coverage_matrix.tsv (gene x sample fraction of canonical-"
                        "transcript bases at depth >= merge_min_dp_snv).")
    p.add_argument("--coverage-depth", type=int, default=50,
                   help="Depth threshold behind --coverage-matrix (merge_min_dp_snv), for labels only.")
    p.add_argument("--zscore-cptm", default=None,
                   help="Accepted for compatibility; not drawn (the report shows MOTR z-scores only).")
    p.add_argument("--zscore-motr", default=None,
                   help="QC: the sample type's by-assignment gene_assignment_zscores_motr.tsv (gene x sample).")
    p.add_argument("--expression-matrices", nargs="*", default=[], metavar="SAMPLE_TYPE:PATH",
                   help="Each sample type's by-assignment gene_assignment_matrix_motr.tsv (gene x sample), "
                        "for the expression box plot beside each gene's bulk view.")
    p.add_argument("--zscore-threshold", type=float, default=-3.5,
                   help="Expression outlier cut-off; cells with |z| >= |this| are outlined (default -3.5).")
    p.add_argument("--pca-coords", default=None,
                   help="QC page: _9C2 {bed_id}_PCA_coords.tsv (sample-type validation PCA coordinates).")
    p.add_argument("--outfile", required=True)
    p.add_argument("--html-out", default=None,
                   help="Also write an offline HTML review page (the same pages as images, with a review "
                        "status and notes per gene) to this path.")
    p.add_argument("--html-dpi", type=int, default=200, help="Resolution of the page images in --html-out.")
    return p.parse_args()


# ----------------------------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------------------------
def _nonempty(v):
    return v is not None and not (isinstance(v, float) and np.isnan(v)) and str(v).strip().lower() not in _MISSING


def _items(cell):
    """';'-separated per-item list, keeping '.' placeholders so parallel columns stay aligned."""
    if not _nonempty(cell):
        return []
    return [x.strip() for x in str(cell).split(";")]


_TISSUE_RE = re.compile(r"\s*([^,()]+?)\s*\(([^)]+)\)")


def _tissue_values(cell):
    """'v1 (tissue1),v2 (tissue2)' -> [(v1, tissue1), (v2, tissue2)]."""
    return [(m.group(1).strip(), m.group(2).strip()) for m in _TISSUE_RE.finditer(str(cell or ""))]


def load_bed(path):
    regions = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or line.startswith(("#", "track", "browser")):
                continue
            regions[parts[3].strip()] = (parts[0], int(parts[1]), int(parts[2]))
    return regions


def load_gene_bams(path):
    out = {}
    if not path or not os.path.isfile(path):
        return out
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for _, r in df.iterrows():
        haps = [r.get("hap1_bam", ""), r.get("hap2_bam", "")]
        if all(h and os.path.isfile(h) for h in haps):
            out[r["gene"]] = haps
    return out


_ATTR_RE = {k: re.compile(k + r' "([^"]+)"') for k in ("gene_name", "transcript_id")}


def load_canonical_transcripts(gtf, genes):
    """gene -> dict(chrom, strand, tid, exons, cds): the Ensembl_canonical transcript, or the
    transcript with the most exonic bases if none is tagged. Coordinates are 0-based half-open."""
    genes = set(genes)
    txs = defaultdict(lambda: defaultdict(lambda: {"exons": [], "cds": [], "canonical": False}))
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9 or parts[2] not in ("exon", "CDS", "transcript"):
                continue
            m = _ATTR_RE["gene_name"].search(parts[8])
            if not m or m.group(1) not in genes:
                continue
            t = _ATTR_RE["transcript_id"].search(parts[8])
            if not t:
                continue
            rec = txs[m.group(1)][t.group(1)]
            rec["chrom"], rec["strand"] = parts[0], parts[6]
            if 'tag "Ensembl_canonical"' in parts[8]:
                rec["canonical"] = True
            if parts[2] == "exon":
                rec["exons"].append((int(parts[3]) - 1, int(parts[4])))
            elif parts[2] == "CDS":
                rec["cds"].append((int(parts[3]) - 1, int(parts[4])))
    out = {}
    for gene, tdict in txs.items():
        cands = {tid: r for tid, r in tdict.items() if r["exons"]}
        if not cands:
            continue
        canon = [tid for tid, r in cands.items() if r["canonical"]]
        tid = canon[0] if canon else max(cands, key=lambda k: sum(e - s for s, e in cands[k]["exons"]))
        r = cands[tid]
        canon_j, all_j = set(), set()
        merged = []                              # every annotated exon of the gene, merged
        for a, b in sorted(e for tr in cands.values() for e in tr["exons"]):
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        for t_id, tr in cands.items():           # junctions of every transcript (same IDs as the pipeline)
            ex = sorted(tr["exons"])
            for (_, e1), (s2, _) in zip(ex, ex[1:]):
                j = f"{tr['chrom']}_{e1 + 1}_{s2}"
                all_j.add(j)
                if tr["canonical"]:
                    canon_j.add(j)
        out[gene] = dict(chrom=r["chrom"], strand=r["strand"], tid=tid,
                         exons=sorted(r["exons"]), cds=sorted(r["cds"]),
                         canonical_junctions=canon_j, annotated_junctions=all_j, all_exons=merged)
    return out


def junction_annotation(tx, jxn):
    if tx is None:
        return "?"
    if jxn in tx["canonical_junctions"]:
        return "canonical"
    return "annotated" if jxn in tx["annotated_junctions"] else "unannotated"


# ----------------------------------------------------------------------------------------------
# Hits -> variants, events, distribution panels
# ----------------------------------------------------------------------------------------------
def parse_variants(row):
    ids = _items(row.get("variant_ID"))
    cols = {c: _items(row.get(c)) for c in row.index if c.startswith("variant_GT_")}
    clnsig = _items(row.get("variant_CLNSIG"))
    out = []
    for i, vid in enumerate(ids):
        parts = vid.rsplit("-", 3)
        if len(parts) != 4 or not parts[1].isdigit():
            continue
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
        gts = sorted({v[i] for v in cols.values() if i < len(v) and _nonempty(v[i])})
        out.append(dict(id=vid, chrom=chrom, pos=pos, ref=ref, alt=alt, gt="/".join(gts) if gts else ".",
                        clnsig=clnsig[i] if i < len(clnsig) and _nonempty(clnsig[i]) else ""))
    return out


def _jxn_coords(j):
    p = j.split("_")
    return "_".join(p[:-2]), int(p[-2]), int(p[-1])


_EID_RE = re.compile(r"\s*\[(E\d+)\]\s*$")


def _split_event_id(text):
    """'label(s) [E2]' -> ('label(s)', 'E2'); no ID -> (text, None)."""
    m = _EID_RE.search(str(text))
    return (str(text)[:m.start()], m.group(1)) if m else (str(text), None)


def junction_calls(row):
    """Every (source, tissue, phasing, junction, label, metric, delta, event_id) call for this gene row."""
    calls = []
    for ph in ("bulk", "hap1", "hap2"):
        jxns = _items(row.get(f"{ph}_jxns"))
        events = _items(row.get(f"{ph}_jxn_event"))
        deltas = _items(row.get(f"{ph}_deltaPSI"))
        for i, j in enumerate(jxns):
            if not _nonempty(j):
                continue
            dmap = dict((t, v) for v, t in _tissue_values(deltas[i] if i < len(deltas) else ""))
            for lab, tissue in _tissue_values(events[i] if i < len(events) else ""):
                lab, eid = _split_event_id(lab)
                for l in lab.split(","):
                    l = l.strip()
                    if l in GTEX_LABELS:
                        calls.append(dict(source="GTEx", tissue=tissue, phasing=ph, junction=j, label=l,
                                          metric="junction_PSI_approx", delta=dmap.get(tissue, ""), eid=eid))
        cj = _items(row.get(f"cohort_{ph}_jxns"))
        cev = _items(row.get(f"cohort_{ph}_jxn_event"))
        cdel = {m: _items(row.get(f"cohort_{ph}_{s}")) for m, s in COHORT_DELTA_SUFFIX.items()}
        for i, j in enumerate(cj):
            if not _nonempty(j) or i >= len(cev):
                continue
            labs, eid = _split_event_id(cev[i])
            for l in labs.split(","):
                l = l.strip()
                metric = COHORT_LABEL_METRIC.get(l)
                if metric is None:
                    continue
                dl = cdel.get(metric, [])
                calls.append(dict(source="cohort", tissue=None, phasing=ph, junction=j, label=l,
                                  metric=metric, delta=dl[i] if i < len(dl) else "", eid=eid))
    return calls


def build_events(calls):
    """Group calls into events. Calls carrying an event ID (written by _9L / _8C) are grouped by
    (source, tissue, phasing, event ID) -- exactly the junctions called together. Calls without
    one (older outputs) fall back to grouping same-label junctions that share a splice site.
    Calls of the same junction set from different sources / tissues / phasings become one event.
    Returns a list of dicts sorted by position."""
    events = {}

    def _add_event(g, labels, src, phasing, cs):
        key = frozenset(g)
        ev = events.setdefault(key, dict(junctions=sorted(g, key=lambda j: _jxn_coords(j)[1:]),
                                         labels=[], calls=[], deltas=defaultdict(list)))
        for label in labels:
            if label not in ev["labels"]:
                ev["labels"].append(label)
        ev["calls"].append(f"{'/'.join(dict.fromkeys(labels))} [{src}, {phasing}]")
        for c in cs:
            if c["junction"] in key and _nonempty(c["delta"]):
                ev["deltas"][c["junction"]].append(c["delta"])

    by_eid = defaultdict(list)
    for c in calls:
        if c.get("eid"):
            by_eid[(c["source"], c["tissue"], c["phasing"], c["eid"])].append(c)
    for (source, tissue, phasing, _eid), cs in by_eid.items():
        src = f"GTEx {tissue}" if source == "GTEx" else "cohort"
        _add_event(sorted({c["junction"] for c in cs}), [c["label"] for c in cs], src, phasing, cs)

    by_key = defaultdict(list)
    for c in calls:
        if not c.get("eid"):
            by_key[(c["source"], c["tissue"], c["phasing"], c["label"])].append(c)
    for (source, tissue, phasing, label), cs in by_key.items():
        jxns = sorted({c["junction"] for c in cs})
        if label in SINGLE_JUNCTION_LABELS:
            groups = [[j] for j in jxns]
        else:
            parent = {j: j for j in jxns}

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            site = defaultdict(list)
            for j in jxns:
                _, s1, s2 = _jxn_coords(j)
                site[("s1", s1)].append(j)
                site[("s2", s2)].append(j)
            for members in site.values():
                for other in members[1:]:
                    parent[find(other)] = find(members[0])
            comp = defaultdict(list)
            for j in jxns:
                comp[find(j)].append(j)
            groups = list(comp.values())
        src = f"GTEx {tissue}" if source == "GTEx" else "cohort"
        for g in groups:
            _add_event(g, [label], src, phasing, cs)
    out = list(events.values())
    out.sort(key=lambda e: min(_jxn_coords(j)[1] for j in e["junctions"]))
    return out


def metric_junction_pairs(calls):
    seen, out = set(), []
    for c in sorted(calls, key=lambda c: (_jxn_coords(c["junction"])[1], c["metric"])):
        key = (c["metric"], c["junction"])
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def load_metric_values(manifests, gene, pairs, coverage_threshold):
    """{(metric, junction): {sample: value}} from this gene's 8A raw TSVs (bulk rows only)."""
    vals = defaultdict(dict)
    needed = {j for _, j in pairs}
    for stype, mpath in manifests.items():
        if not os.path.isfile(mpath):
            continue
        man = pd.read_csv(mpath, sep="\t", dtype=str, keep_default_na=False)
        paths = man.loc[man["gene"] == gene, "result_path"].tolist()
        for rp in paths:
            rp = str(rp).strip()
            if not rp or rp in ("None", "nan", "NA") or not os.path.isfile(rp):
                continue
            df = pd.read_csv(rp, sep="\t", dtype=str)
            if not {"sample", "phasing", "junction"} <= set(df.columns):
                continue
            df = df[(df["phasing"] == "bulk") & df["junction"].isin(needed)]
            for metric, jxn in pairs:
                rc, dc = _METRIC_RESCALED_COL.get(metric), _METRIC_DENOM_COL.get(metric)
                if rc not in df.columns or dc not in df.columns:
                    continue
                sub = df[df["junction"] == jxn]
                v = pd.to_numeric(sub[rc], errors="coerce")
                d = pd.to_numeric(sub[dc], errors="coerce")
                ok = (d >= coverage_threshold) & v.notna()
                for s, x in zip(sub["sample"][ok], v[ok]):
                    vals[(metric, jxn)][str(s)] = float(x)
    return vals


# ----------------------------------------------------------------------------------------------
# Reads
# ----------------------------------------------------------------------------------------------
def _keep(read):
    return not (read.flag & SKIP_FLAGS)


def _is_spliced(read):
    return any(op == 3 for op, _ in (read.cigartuples or ()))


def _keep_spliced(read):
    return _keep(read) and _is_spliced(read)


EXON_TOL = 10    # bp an unspliced read may extend past an annotated exon and still count as exonic


def read_filter(tx):
    """Read filter for a gene's sashimi tracks. A multi-exon gene keeps
    spliced reads plus unspliced reads lying within an annotated exon (any transcript, +/- EXON_TOL),
    so reads inside a long exon are kept while unspliced pre-mRNA / intronic reads are dropped.
    A single-exon gene (or no transcript) keeps every read."""
    if tx is None or len(tx["exons"]) <= 1:
        return _keep
    exons = tx.get("all_exons") or sorted(tx["exons"])
    starts = [a for a, _ in exons]

    def keep(read):
        if not _keep(read):
            return False
        if _is_spliced(read):
            return True
        rs, re_ = read.reference_start, read.reference_end or read.reference_start
        i = bisect.bisect_right(starts, rs + EXON_TOL) - 1
        return i >= 0 and exons[i][0] - EXON_TOL <= rs and re_ <= exons[i][1] + EXON_TOL
    return keep


def sample_reads(bam, chrom, start, end, max_reads, seed, keep=_keep):
    """Reservoir-sample up to max_reads reads overlapping the region; returns (reads, n_total)."""
    rng = random.Random(seed)
    res, n = [], 0
    for r in bam.fetch(chrom, max(0, start), end):
        if not keep(r):
            continue
        n += 1
        if len(res) < max_reads:
            res.append(r)
        else:
            k = rng.randrange(n)
            if k < max_reads:
                res[k] = r
    res.sort(key=lambda r: (r.reference_start, r.reference_end or 0))
    return res, n


def coverage_counts(bam, chrom, start, end, keep=_keep):
    """(4, L) array of A/C/G/T counts per position from all kept reads."""
    acgt = bam.count_coverage(chrom, max(0, start), end, quality_threshold=0, read_callback=keep)
    return np.array(acgt, dtype=np.int64)


def junction_read_counts(bam, chrom, start, end, junctions):
    want = {}
    for j in junctions:
        _, s1, s2 = _jxn_coords(j)
        want[(s1, s2)] = j
    counts = {j: 0 for j in junctions}
    for r in bam.fetch(chrom, max(0, start), end):
        if not _keep(r) or r.cigartuples is None:
            continue
        pos = r.reference_start
        for op, ln in r.cigartuples:
            if op in (0, 2, 7, 8):
                pos += ln
            elif op == 3:
                j = want.get((pos + 1, pos + ln))
                if j is not None:
                    counts[j] += 1
                pos += ln
    return counts


# ----------------------------------------------------------------------------------------------
# Drawing
# ----------------------------------------------------------------------------------------------
def _fmt_pos(x, _=None):
    return f"{int(x):,}"


def draw_coverage(ax, acgt, ref_seq, start, end, label, arcs=None, show_xticks=True):
    depth = acgt.sum(axis=0)
    L = end - start
    x = np.arange(start, end)
    if L > 4000:                                  # bin long regions to keep the PDF light
        nb = 2000
        edges = np.linspace(0, L, nb + 1).astype(int)
        mx = np.maximum.reduceat(depth, edges[:-1]) if L else depth
        xs = np.repeat(start + edges, 2)[1:-1]
        ys = np.repeat(mx, 2)
        ax.fill_between(xs, ys, color=COV_COLOR, linewidth=0, rasterized=True)
    else:
        ax.bar(x + 0.5, depth, width=1.0, color=COV_COLOR, linewidth=0)
    top = max(int(depth.max()) if depth.size else 0, 1)

    # positions where a non-reference base exceeds the AF threshold: stacked base-coloured bars
    if ref_seq:
        ref_idx = np.array(["ACGT".find(b) for b in ref_seq.upper()])
        refcount = np.where(ref_idx >= 0, acgt[np.clip(ref_idx, 0, 3), np.arange(L)], 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            af = np.where(depth > 0, (depth - refcount) / depth, 0)
        hits = np.where((af > AF_THRESHOLD) & (depth >= AF_MIN_DEPTH))[0]
        w = max(1.0, L / 500.0)
        for i in hits:
            bottom = 0
            for b, base in enumerate("ACGT"):
                c = acgt[b, i]
                if c:
                    ax.add_patch(Rectangle((start + i + 0.5 - w / 2, bottom), w, c,
                                           color=BASE_COLORS[base], linewidth=0))
                    bottom += c

    ymax = top
    if arcs:
        ymax = top * (1.25 + 0.3 * len(arcs))
        maxc = max([c for _, _, c in arcs] + [1])
        for k, (a, b, c) in enumerate(sorted(arcs, key=lambda t: t[1] - t[0])):
            h = top * (1.1 + 0.3 * (k + 1))
            ya = depth[min(max(a - start, 0), L - 1)] if L else 0
            yb = depth[min(max(b - start, 0), L - 1)] if L else 0
            xa, xb = a + 0.5, b + 0.5          # centres of the flanking exon bases
            path = Path([(xa, ya), (xa, h), (xb, h), (xb, yb)],
                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
            ax.add_patch(PathPatch(path, facecolor="none", edgecolor=ARC_COLOR,
                                   linewidth=0.6 + 2.4 * c / maxc, clip_on=False))
            ax.text((xa + xb) / 2, h * 0.78 + 0.22 * max(ya, yb), str(c), ha="center", va="bottom",
                    fontsize=5 * FS, color=ARC_COLOR,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
    ax.set_xlim(start, end)
    ax.set_ylim(0, ymax)
    ax.set_yticks([])
    ax.text(0.003, 0.97, f"[0-{top:,}]", transform=ax.transAxes, va="top", ha="left", fontsize=5 * FS)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)
    if show_xticks:
        ax.xaxis.tick_top()
        ax.xaxis.set_major_locator(plt.MaxNLocator(4, integer=True))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(_fmt_pos))
        ax.tick_params(axis="x", labelsize=5 * FS, length=2, pad=1)
    else:
        ax.set_xticks([])
    if label:
        _track_label(ax, label)


def _track_label(ax, label):
    """Left-hand track label; for 'name\nrest' (e.g. 'bulk\ncoverage') the track name is bold."""
    if "\n" not in label:
        ax.set_ylabel(label, fontsize=6 * FS, rotation=0, ha="right", va="center", labelpad=4)
        return
    name, rest = label.split("\n", 1)
    tr = matplotlib.transforms.offset_copy(ax.transAxes, fig=ax.figure, x=-4, y=0, units="points")
    ax.text(0, 0.5, name, transform=tr, fontsize=6 * FS, fontweight="bold", ha="right", va="bottom")
    ax.text(0, 0.5, rest, transform=tr, fontsize=6 * FS, ha="right", va="top")


def indel_specs(variants):
    """(kind, 0-based reference position, length) of every candidate indel, so a read's matching
    indel is drawn even when shorter than MIN_INDEL. Deletion: the deleted bases start at VCF POS
    (0-based), length len(REF) - len(ALT). Insertion: before 0-based POS, length len(ALT) - len(REF)."""
    out = []
    for v in variants:
        d = len(v["alt"]) - len(v["ref"])
        if d < 0:
            out.append(("D", v["pos"], -d))
        elif d > 0:
            out.append(("I", v["pos"], d))
    return out


def _is_variant_indel(kind, pos, ln, specs):
    return any(k == kind and l == ln and abs(p - pos) <= INDEL_POS_TOL for k, p, l in specs)


def _rgb(c):
    return tuple(int(round(v * 255)) for v in matplotlib.colors.to_rgb(c))


def draw_alignments(ax, reads, n_total, ref_seq, start, end, label, keep_indels=()):
    """Paint the alignment track straight into an RGB image sized to the axes at ALIGN_DPI and
    embed it unresampled -- far sharper (and lighter on memory) than rasterising thousands of
    matplotlib patches onto a full-page canvas."""
    L = max(end - start, 1)
    gap = max(1, int(L * 0.003))
    row_ends, placed = [], []
    for r in reads:                                    # pack reads into rows first
        rs, re_ = r.reference_start, r.reference_end or r.reference_start
        row = next((i for i, e in enumerate(row_ends) if e + gap <= rs), None)
        if row is None:
            row = len(row_ends)
            row_ends.append(re_)
        else:
            row_ends[row] = re_
        placed.append((r, row))
    nrows = max(len(row_ends), 8)

    pos_ax = ax.get_position()
    W = max(int(round(pos_ax.width * PAGE_W * ALIGN_DPI)), 1)
    H = max(int(round(pos_ax.height * PAGE_H * ALIGN_DPI)), 1)
    img = np.full((H, W, 3), 255, dtype=np.uint8)
    row_px = H / nrows
    # a gap between rows only when it is at least one pixel; sub-pixel gaps cause moire stripes
    padpx = int(row_px * 0.1) if row_px >= 10 else 0
    px_per_base = W / L
    read_c, gap_c, del_c, ins_c = _rgb(READ_COLOR), _rgb(GAP_COLOR), _rgb(DEL_COLOR), _rgb(INS_COLOR)
    base_c = {b: _rgb(c) for b, c in BASE_COLORS.items()}

    def xcols(a, b):
        c0 = int((a - start) * px_per_base)
        c1 = max(int(round((b - start) * px_per_base)), c0 + 1)
        return max(c0, 0), min(c1, W)

    line_h = max(1, int(round(row_px * 0.15)))
    ref_up = ref_seq.upper() if ref_seq else ""
    mm_w = max(1, int(round(px_per_base)), int(W / 1500))   # mismatches stay visible when zoomed out
    for r, row in placed:
        r0 = int(round(row * row_px)) + padpx
        r1 = max(int(round((row + 1) * row_px)) - padpx, r0 + 1)
        mid = (r0 + r1) // 2
        l0, l1 = max(mid - line_h // 2, r0), min(mid - line_h // 2 + line_h, r1)
        seq = r.query_sequence
        pos, qpos = r.reference_start, 0
        mism = []
        for op, ln in r.cigartuples or []:
            if op in (0, 7, 8):
                a, b = max(pos, start), min(pos + ln, end)
                if a < b:
                    c0, c1 = xcols(a, b)
                    img[r0:r1, c0:c1] = read_c
                    if seq and ref_up:
                        q = np.frombuffer(seq[qpos + (a - pos): qpos + (b - pos)].upper().encode(), dtype="S1")
                        rf = np.frombuffer(ref_up[a - start: b - start].encode(), dtype="S1")
                        for k in np.nonzero(q != rf)[0]:
                            mism.append((a + k, q[k].decode()))
                pos += ln
                qpos += ln
            elif op == 1:
                if start <= pos <= end and (ln >= MIN_INDEL or _is_variant_indel("I", pos, ln, keep_indels)):
                    c = min(int((pos - start) * px_per_base), W - 1)
                    img[r0:r1, max(c - 1, 0):c + 1] = ins_c
                qpos += ln
            elif op == 2:
                a, b = max(pos, start), min(pos + ln, end)
                if a < b:
                    c0, c1 = xcols(a, b)
                    if ln >= MIN_INDEL or _is_variant_indel("D", pos, ln, keep_indels):
                        img[l0:l1, c0:c1] = del_c
                    else:                       # small deletion: hidden (drawn as part of the read)
                        img[r0:r1, c0:c1] = read_c
                pos += ln
            elif op == 3:
                a, b = max(pos, start), min(pos + ln, end)
                if a < b:
                    c0, c1 = xcols(a, b)
                    img[l0:l1, c0:c1] = gap_c
                pos += ln
            elif op == 4:
                qpos += ln
        for x, base in mism:                    # drawn last so neighbouring reads don't cover them
            c = int((x + 0.5 - start) * px_per_base)
            c0 = max(c - mm_w // 2, 0)
            img[r0:r1, c0:min(c0 + mm_w, W)] = base_c.get(base, base_c["N"])

    ax.imshow(img, extent=(start, end, nrows, 0), aspect="auto", interpolation="none")
    ax.set_xlim(start, end)
    ax.set_ylim(nrows, 0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if label:
        ax.set_ylabel(label, fontsize=6 * FS, rotation=0, ha="right", va="center", labelpad=4)


def draw_reference(ax, ref_seq, start, end, highlight=None):
    for i, b in enumerate(ref_seq.upper()):
        ax.text(start + i + 0.5, 0.5, b, ha="center", va="center", fontsize=5.5 * FS,
                color=BASE_COLORS.get(b, "#555"), fontweight="bold", family="monospace")
    if highlight:
        a, b = highlight
        ax.add_patch(Rectangle((a, 0.02), b - a, 0.96, fill=False, edgecolor="black", linewidth=0.6))
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.axis("off")


def draw_transcript(ax, tx, start, end, label=None, marks=None, spans=None, xform=None):
    """Canonical transcript track. xform: optional genomic -> display coordinate map (sashimi
    intron scaling); positions are tested in genomic space and drawn in display space."""
    X = (lambda v: float(xform(v))) if xform is not None else (lambda v: v)
    ax.set_xlim(X(start), X(end))
    ax.set_ylim(0, 1)
    ax.axis("off")
    for a, b in spans or []:
        ax.axvspan(X(a), X(b), color=SPAN_COLOR, zorder=0)
    if tx is None:
        ax.text(0.5, 0.5, "no transcript in GTF", transform=ax.transAxes, ha="center", va="center", fontsize=5 * FS)
        return
    exons, cds = tx["exons"], tx["cds"]
    t0, t1 = exons[0][0], exons[-1][1]
    a, b = max(t0, start), min(t1, end)
    if a < b:
        ax.plot([X(a), X(b)], [0.5, 0.5], color="#1f3a93", linewidth=0.6, zorder=1)
        # strand arrows evenly along the transcript: blue on the intron line, white inside exons
        n_arrows = max(3, int(round(ax.get_position().width * PAGE_W * 1.3)))
        marker = ">" if tx["strand"] == "+" else "<"
        if xform is not None:      # evenly spaced on the display axis
            disp = np.linspace(X(a), X(b), n_arrows + 2)[1:-1]
            gpos = np.interp(disp, xform(np.arange(a, b + 1)), np.arange(a, b + 1))
        else:
            gpos = np.linspace(a, b, n_arrows + 2)[1:-1]
        for xa in gpos:
            in_exon = any(s <= xa < e for s, e in exons)
            ax.plot(X(xa), 0.5, marker=marker, markersize=2.2,
                    color="white" if in_exon else "#1f3a93", zorder=5 if in_exon else 2)
    for s, e in exons:
        if e <= start or s >= end:
            continue
        ax.add_patch(Rectangle((X(s), 0.32), X(e) - X(s), 0.36, color="#1f3a93", linewidth=0, zorder=3))
    for s, e in cds:
        if e <= start or s >= end:
            continue
        ax.add_patch(Rectangle((X(s), 0.15), X(e) - X(s), 0.70, color="#1f3a93", linewidth=0, zorder=3))
    for x in marks or []:
        ax.plot(X(x), 0.95, marker="v", markersize=3, color="#d62728", zorder=4, clip_on=False)
    if label:
        ax.text(-0.003, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=5.5 * FS)


def draw_distribution(ax, metric, junction, values, sample_to_type, type_colors, outliers, sample):
    if not values:
        ax.text(0.5, 0.5, "no samples with\nenough coverage", ha="center", va="center", fontsize=5 * FS,
                transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        types = sorted({sample_to_type.get(s, "unknown") for s in values})
        rng = np.random.default_rng(42)
        for i, t in enumerate(types, 1):
            ss = [s for s in values if sample_to_type.get(s, "unknown") == t]
            vs = [values[s] for s in ss]
            black = dict(color="black", linewidth=0.6)
            ax.boxplot([vs], positions=[i], widths=0.5, showfliers=False, medianprops=dict(color="black", linewidth=1),
                       boxprops=black, whiskerprops=black, capprops=black)
            col = type_colors.get(t, "#555555")
            xs = rng.normal(i, 0.06, len(vs))
            for x, y, s in zip(xs, vs, ss):
                if s == sample:
                    ax.scatter(x, y, marker="*", s=50, color=col, edgecolors="black", linewidths=0.5, zorder=4)
                    ax.annotate(s, (x, y), xytext=(4, 0), textcoords="offset points",
                                fontsize=5 * FS, va="center", fontweight="bold")
                elif s in outliers:         # other samples' outlier calls: outlined, not labelled
                    ax.scatter(x, y, s=9, color=col, edgecolors="black", linewidths=0.4, zorder=3)
                else:
                    ax.scatter(x, y, s=6, color=col, alpha=0.8, linewidths=0, zorder=2)
        ax.set_xticks(range(1, len(types) + 1))
        ax.set_xticklabels(types, fontsize=5 * FS)
        ax.set_xlim(0.4, len(types) + 0.6)
        ax.set_ylim(-0.02, 1.05)
        ax.tick_params(axis="y", labelsize=5 * FS, length=2)
    ax.set_title(f"{METRIC_SHORT.get(metric, metric)} | {junction}", fontsize=5.5 * FS, pad=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ----------------------------------------------------------------------------------------------
# Layout
# ----------------------------------------------------------------------------------------------
def _flush_pdf_images(pdf):
    """Write the images of the pages saved so far into the PDF now and drop them from memory.

    matplotlib's PDF backend keeps every embedded image (each alignment track is a ~90 MB RGBA
    array at ALIGN_DPI) until the whole file is closed, so a report with many genes ran out of
    memory. Between pages no content stream is open, so the pending image objects can be written
    immediately; the image registry keeps only each object's name and reference (under a key that
    cannot collide with a later image), which is all the file's final bookkeeping needs. If this
    matplotlib version doesn't have the expected internals, nothing is done (old behaviour)."""
    f = getattr(pdf, "_file", None)
    needed = ("_images", "_unpack", "_writeImg", "reserveObject", "writeImages")
    if f is None or not all(hasattr(f, a) for a in needed) or getattr(f, "currentstream", None) is not None:
        return
    if not getattr(f, "_hit_report_patched", False):
        orig = f.writeImages

        def write_remaining_images(_f=f, _orig=orig):
            pending = {k: v for k, v in _f._images.items() if v[0] is not None}
            done = {k: v for k, v in _f._images.items() if v[0] is None}
            _f._images = pending
            try:
                _orig()
            finally:
                _f._images = {**done, **pending}
        f.writeImages = write_remaining_images
        f._hit_report_patched = True
    for key, (img, name, ob) in list(f._images.items()):
        if img is None:
            continue
        data, adata = f._unpack(img)
        smask = None
        if adata is not None:
            smask = f.reserveObject("smask")
            f._writeImg(adata, smask.id)
        f._writeImg(data, ob.id, smask)
        del f._images[key]
        f._images[("written", name)] = (None, name, ob)


REVIEW_CTX = ("overview", None)          # what the page being drawn belongs to (for --html-out)


def link_text(fig, x_in, y_in, text, target, ha="right"):
    """Clickable header text at (x, y) inches from the bottom-left; `target` is resolved to a page
    after the PDF is written (see add_pdf_navigation)."""
    t = fig.text(x_in / PAGE_W, y_in / PAGE_H, text, fontsize=6 * FS, color="#0068a9", ha=ha, va="top")
    t._pdf_link = target
    return t


def _link_rects(fig):
    """(rect in PDF points, target) of every link_text on the figure."""
    texts = [t for t in fig.texts if getattr(t, "_pdf_link", None) is not None]
    if not texts:
        return []
    r = fig.canvas.get_renderer()
    k = 72.0 / fig.dpi
    out = []
    for t in texts:
        bb = t.get_window_extent(renderer=r)
        out.append(((bb.x0 * k - 1, bb.y0 * k - 1, bb.x1 * k + 1, bb.y1 * k + 1), t._pdf_link))
    return out


def nav_links(fig, prev_gene=None, next_gene=None):
    """'« previous gene: X    contents    next gene: Y »' at the top right of a page."""
    y, x = PAGE_H - 0.16, PAGE_W - MARGIN_R
    items = ([(f"next gene: {next_gene} \u00bb", ("gene", next_gene))] if next_gene else []) \
        + [("contents", ("contents",))] \
        + ([(f"\u00ab previous gene: {prev_gene}", ("gene", prev_gene))] if prev_gene else [])
    r = fig.canvas.get_renderer()
    for text, target in items:                       # right to left, each placed by the width of the last
        t = link_text(fig, x, y, text, target)
        x -= t.get_window_extent(renderer=r).width / fig.dpi + 0.25


def draw_contents(pdf, shown, genes):
    """Contents pages, drawn last (when every page number is known) and moved to just after the
    overview by add_pdf_navigation. Every row links to its first page."""
    global REVIEW_CTX
    fp_ = pdf.first_page
    qc_sections = sorted(((n, i) for (k, n), i in fp_.items() if k == "qc"), key=lambda x: x[1])
    rows = [("Overview", None, None, None, ("overview", None))]
    rows += [(g["gene"], g["rank"], g["tier"], g["hits"], ("gene", g["gene"])) for g in genes
             if ("gene", g["gene"]) in fp_]
    rows += [(n, None, None, None, ("qc", n)) for n, _ in qc_sections]
    per_page = 30
    n_cont = max(1, -(-len(rows) // per_page))
    REVIEW_CTX = ("contents", None)
    for k in range(n_cont):
        fig = plt.figure(figsize=(PAGE_W, PAGE_H))
        fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H, f"{shown}  |  Contents" + (f"  (page {k + 1})" if k else ""),
                 fontsize=11 * FS, fontweight="bold", va="top")
        fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H, "click a name to jump to its pages; every page links "
                 "back here (top right)", fontsize=6.5 * FS, va="top", color="#333")
        y = PAGE_H - 0.95
        for x, h in ((MARGIN_L, "Rank"), (MARGIN_L + 0.5, "Gene / section"), (MARGIN_L + 2.3, "Tier"),
                     (MARGIN_L + 2.8, "Hits"), (PAGE_W - MARGIN_R - 0.4, "Page")):
            fig.text(x / PAGE_W, y / PAGE_H, h, fontsize=6.5 * FS, fontweight="bold", va="top")
        y -= 0.24
        prev_kind = "overview"
        for name, rank, tier, hits, target in rows[k * per_page:(k + 1) * per_page]:
            if target[0] != prev_kind:
                y -= 0.06
                prev_kind = target[0]
            i = fp_.get(target)
            if rank:
                fig.text(MARGIN_L / PAGE_W, y / PAGE_H, rank, fontsize=6.5 * FS, va="top")
            link_text(fig, MARGIN_L + 0.5, y, name, target, ha="left").set_fontsize(6.5 * FS)
            if tier:
                fig.text((MARGIN_L + 2.3) / PAGE_W, y / PAGE_H, tier, fontsize=6.5 * FS, va="top")
            if hits:
                h = hits if len(hits) <= 95 else hits[:92] + "..."
                fig.text((MARGIN_L + 2.8) / PAGE_W, y / PAGE_H, h, fontsize=5.5 * FS, va="top", color="#333")
            if i is not None:                        # final page number: contents pages sit after page 1
                fig.text((PAGE_W - MARGIN_R - 0.4) / PAGE_W, y / PAGE_H, str(i + 1 + (n_cont if i > 0 else 0)),
                         fontsize=6.5 * FS, va="top")
            y -= 0.205
        pdf.savefig(fig)
        plt.close(fig)
    return rows


def add_pdf_navigation(path, cap, rows):
    """Finish the PDF in place with an incremental update (no extra libraries): move the contents
    pages to just after the overview, turn every link_text into a link to its target's first page,
    and add bookmarks. Relies on matplotlib's plain PDF layout (catalog, page tree with /Kids, an
    /Annots object per page); if that is not what it finds, the PDF is left as written."""
    with open(path, "rb") as fh:
        data = fh.read()

    def obj(num):
        m = re.search(rb"(?<![0-9])%d 0 obj\s*(.*?)\s*endobj" % num, data, re.S)
        return m.group(1) if m else None

    try:
        tr = data[data.rindex(b"trailer"):]
        size = int(re.search(rb"/Size (\d+)", tr).group(1))
        root = int(re.search(rb"/Root (\d+) 0 R", tr).group(1))
        info = re.search(rb"/Info (\d+) 0 R", tr)
        prev = int(re.search(rb"startxref\s+(\d+)", data[data.rindex(b"startxref"):]).group(1))
        pages_num = int(re.search(rb"/Pages (\d+) 0 R", obj(root)).group(1))
        kids = [int(x) for x in re.findall(rb"(\d+) 0 R", re.search(rb"/Kids \[(.*?)\]", obj(pages_num), re.S).group(1))]
        annots = [int(re.search(rb"/Annots (\d+) 0 R", obj(k)).group(1)) for k in kids]
    except Exception as e:                           # unexpected layout: leave the file alone
        print(f"  [WARNING] PDF navigation not added ({e!r}).")
        return
    if len(kids) != cap.n:
        print(f"  [WARNING] PDF navigation not added (expected {cap.n} pages, found {len(kids)}).")
        return

    cont = [i for i in range(cap.n) if cap.page_ctx[i][0] == "contents"]
    order = [0] + cont + [i for i in range(1, cap.n) if i not in cont]
    pos = {orig: new for new, orig in enumerate(order)}           # original index -> final index

    def dest(target):
        i = 0 if target[0] == "overview" else (cont[0] if target[0] == "contents" and cont else
                                                cap.first_page.get(target))
        return None if i is None else kids[i]

    new_objs = {}
    nxt = [size]

    def add(body):
        n = nxt[0]
        nxt[0] += 1
        new_objs[n] = body
        return n

    per_page = defaultdict(list)
    for i, (x0, y0, x1, y1), target in cap.links:
        d = dest(target)
        if d is not None:
            a = add(b"<< /Type /Annot /Subtype /Link /Rect [ %.2f %.2f %.2f %.2f ] /Border [ 0 0 0 ] "
                    b"/Dest [ %d 0 R /Fit ] >>" % (x0, y0, x1, y1, d))
            per_page[i].append(a)
    for i, refs in per_page.items():
        new_objs[annots[i]] = b"[ " + b" ".join(b"%d 0 R" % r for r in refs) + b" ]"
    old_pages = obj(pages_num)
    new_objs[pages_num] = re.sub(rb"/Kids \[.*?\]", b"/Kids [ " + b" ".join(b"%d 0 R" % kids[i] for i in order)
                                 + b" ]", old_pages, flags=re.S)

    # bookmarks: Overview, Contents, Candidate genes > each gene, QC and expression > each section
    def esc(t):
        return t.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)").encode("latin-1", "replace")
    tops = [("Overview", ("overview", None), [])]
    if cont:
        tops.append(("Contents", ("contents",), []))
    genes_ = [(f"{r[1]}. {r[0]} (tier {r[2]})", r[4]) for r in rows if r[4][0] == "gene" and dest(r[4]) is not None]
    qcs = [(r[0], r[4]) for r in rows if r[4][0] == "qc" and dest(r[4]) is not None]
    if genes_:
        tops.append(("Candidate genes", genes_[0][1], genes_))
    if qcs:
        tops.append(("QC and expression", qcs[0][1], qcs))
    outl = nxt[0]; nxt[0] += 1
    top_nums = [nxt[0] + k for k in range(len(tops))]; nxt[0] += len(tops)
    for k, (title, target, children) in enumerate(tops):
        ch_nums = [nxt[0] + c for c in range(len(children))]; nxt[0] += len(children)
        for c, (ct, ctarget) in enumerate(children):
            body = b"<< /Title (%s) /Parent %d 0 R /Dest [ %d 0 R /Fit ]" % (esc(ct), top_nums[k], dest(ctarget))
            if c > 0:
                body += b" /Prev %d 0 R" % ch_nums[c - 1]
            if c + 1 < len(children):
                body += b" /Next %d 0 R" % ch_nums[c + 1]
            new_objs[ch_nums[c]] = body + b" >>"
        body = b"<< /Title (%s) /Parent %d 0 R /Dest [ %d 0 R /Fit ]" % (esc(title), outl, dest(target))
        if k > 0:
            body += b" /Prev %d 0 R" % top_nums[k - 1]
        if k + 1 < len(tops):
            body += b" /Next %d 0 R" % top_nums[k + 1]
        if children:
            body += b" /First %d 0 R /Last %d 0 R /Count %d" % (ch_nums[0], ch_nums[-1], len(children))
        new_objs[top_nums[k]] = body + b" >>"
    new_objs[outl] = b"<< /Type /Outlines /First %d 0 R /Last %d 0 R /Count %d >>" % (
        top_nums[0], top_nums[-1], len(tops))
    cat = obj(root).rstrip()
    assert cat.endswith(b">>")
    new_objs[root] = cat[:-2].rstrip() + b" /Outlines %d 0 R /PageMode /UseOutlines >>" % outl

    out = bytearray(data if data.endswith(b"\n") else data + b"\n")
    offsets = {}
    for n in sorted(new_objs):
        offsets[n] = len(out)
        out += b"%d 0 obj\n" % n + new_objs[n] + b"\nendobj\n"
    xref = len(out)
    out += b"xref\n0 1\n0000000000 65535 f \n"
    for n in sorted(offsets):
        out += b"%d 1\n%010d 00000 n \n" % (n, offsets[n])
    out += b"trailer\n<< /Size %d /Root %d 0 R" % (nxt[0], root)
    if info:
        out += b" /Info %s 0 R" % info.group(1)
    out += b" /Prev %d >>\nstartxref\n%d\n%%%%EOF\n" % (prev, xref)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(out)
    os.replace(tmp, path)
    print(f"  contents page, bookmarks and {sum(len(v) for v in per_page.values())} page link(s) added")


def _compact_png(png):
    """The page PNG reduced to a 256-colour palette (about a third of the size, no visible change
    for these plots); unchanged if Pillow (a matplotlib dependency) is unavailable."""
    try:
        from PIL import Image
    except ImportError:
        return png
    im = Image.open(io.BytesIO(png)).convert("RGB")
    q = im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    out = io.BytesIO()
    q.save(out, format="PNG", optimize=True)
    return out.getvalue() if out.tell() < len(png) else png


class CapturePdf:
    """PdfPages stand-in that also keeps each saved page as a PNG (for the HTML review page),
    tagged with what it belongs to: the figure's _review_ctx (set when a gene page is created), or
    REVIEW_CTX at the time it is saved."""

    def __init__(self, pdf, dpi):
        self._pdf, self._dpi, self.pages = pdf, dpi, []
        self.n = 0                    # pages saved so far
        self.first_page = {}          # review context -> its first page index
        self.links = []               # (page index, (x0, y0, x1, y1) in PDF points, target)
        self.page_ctx = []            # review context of every page, in order

    def savefig(self, fig, *a, **k):
        ctx = getattr(fig, "_review_ctx", REVIEW_CTX)
        self.first_page.setdefault(ctx, self.n)
        self.page_ctx.append(ctx)
        self.links += [(self.n, r, t) for r, t in _link_rects(fig)]
        self._pdf.savefig(fig, *a, **k)
        self.n += 1
        if self._dpi and ctx[0] != "contents":
            nav = [t for t in fig.texts if getattr(t, "_pdf_link", None) is not None]
            for t in nav:                 # PDF-only links: not clickable in the HTML page images
                t.set_visible(False)
            # the alignment tracks are 1200-dpi images: drawn nearest-neighbour at screen resolution
            # they alias (reads drop out, splice lines turn into stripes), so smooth them for the PNG
            ims = [im for ax in fig.axes for im in ax.images]
            old = [im.get_interpolation() for im in ims]
            for im in ims:
                im.set_interpolation("antialiased")
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=self._dpi)
            for im, o in zip(ims, old):
                im.set_interpolation(o)
            for t in nav:
                t.set_visible(True)
            self.pages.append((ctx, _compact_png(buf.getvalue())))

    def __getattr__(self, name):
        return getattr(self._pdf, name)


class Pager:
    """Hands out row slots, 4 per page, with a page header naming the gene."""

    def __init__(self, pdf, header_fn):
        self.pdf, self.header_fn = pdf, header_fn
        self.fig, self.row = None, ROWS_PER_PAGE
        self.part = 0

    def _new_page(self):
        self.flush()
        self.fig = plt.figure(figsize=(PAGE_W, PAGE_H))
        self.fig._review_ctx = REVIEW_CTX
        self.part += 1
        self.header_fn(self.fig, self.part)
        self.row = 0

    def fresh_page(self):
        """Start the next row on a new page (unless the current page is still empty)."""
        if self.fig is not None and self.row > 0:
            self.row = ROWS_PER_PAGE

    def new_gene(self, header_fn):
        self.flush()
        self.header_fn, self.part, self.row = header_fn, 0, ROWS_PER_PAGE

    def next_row(self):
        if self.row >= ROWS_PER_PAGE:
            self._new_page()
        h = (PAGE_H - MARGIN_T - MARGIN_B) / ROWS_PER_PAGE
        top = PAGE_H - MARGIN_T - self.row * h
        self.row += 1
        return self.fig, top, h

    def flush(self):
        if self.fig is not None:
            self.pdf.savefig(self.fig, dpi=RASTER_DPI)
            _flush_pdf_images(self.pdf)
            plt.close(self.fig)
            self.fig = None


def _ax(fig, left, bottom, width, height):
    return fig.add_axes([left / PAGE_W, bottom / PAGE_H, width / PAGE_W, height / PAGE_H])


def panel_columns(n, per_row):
    """(left, width) in inches for n panels, each 1/per_row of the content width."""
    total = PAGE_W - MARGIN_L - MARGIN_R
    w = (total - PANEL_GAP * (per_row - 1)) / per_row
    return [(MARGIN_L + i * (w + PANEL_GAP), w) for i in range(n)]


def alignment_panel(fig, left, width, top, h, title, bam, chrom, start, end, fasta, tx, args, seed,
                    label_tracks=True, arcs_junctions=None, ref_track=False, highlight=None,
                    marks=None, spans=None, track_label="", show_tx=True, title_width=None, keep_indels=()):
    max_chars = int((title_width or width) * 19 / FS)
    if len(title) > max_chars:
        title = title[:max_chars - 3] + "..."
    title_h, tick_h, cov_h, ref_h, pad = 0.19, 0.13, 0.36, (0.15 if ref_track else 0.0), 0.03
    ann_h = 0.20 if show_tx else 0.0
    aln_h = h - title_h - tick_h - cov_h - ann_h - ref_h - (4 if show_tx else 3) * pad - 0.06
    y = top - title_h
    if title:
        fig.text(left / PAGE_W, (y + 0.02) / PAGE_H, title, fontsize=6.5 * FS, fontweight="bold", va="bottom", ha="left")
    y -= tick_h
    ref_seq = fasta.fetch(chrom, max(0, start), end) if fasta else ""
    if len(ref_seq) < end - max(0, start):
        ref_seq = ref_seq + "N" * (end - max(0, start) - len(ref_seq))
    start = max(0, start)
    # IGV views show every read (unmapped, secondary, QC-fail and duplicate reads aside)
    keep = _keep
    acgt = coverage_counts(bam, chrom, start, end, keep)
    arcs = None
    if arcs_junctions:
        counts = junction_read_counts(bam, chrom, start, end, arcs_junctions)
        # arc from the last exon base before the intron (ss1 - 1, 1-based) to the first exon base
        # after it (ss2 + 1, 1-based), as 0-based base indices
        arcs = [(_jxn_coords(j)[1] - 2, _jxn_coords(j)[2], counts[j]) for j in arcs_junctions]
    ax = _ax(fig, left, y - cov_h, width, cov_h)
    draw_coverage(ax, acgt, ref_seq, start, end, f"{track_label}\ncoverage" if label_tracks else "", arcs=arcs)
    y -= cov_h + pad
    reads, n = sample_reads(bam, chrom, start, end, args.max_reads, seed, keep)
    ax = _ax(fig, left, y - aln_h, width, aln_h)
    draw_alignments(ax, reads, n, ref_seq, start, end, "alignments" if label_tracks else "",
                    keep_indels=keep_indels)
    y -= aln_h + pad
    if ref_track:
        ax = _ax(fig, left, y - ref_h, width, ref_h)
        draw_reference(ax, ref_seq, start, end, highlight=highlight)
        y -= ref_h + pad
    if show_tx:
        ax = _ax(fig, left, y - ann_h, width, ann_h)
        draw_transcript(ax, tx, start, end, label=None,
                        marks=marks, spans=spans)


# ----------------------------------------------------------------------------------------------
# Sashimi plots (outlier splicing events) and the representative cohort sample
# ----------------------------------------------------------------------------------------------
REP_MIN_JUNCTIONS = 3        # a candidate must share this many covered junctions with the median
REP_DEPTH_MARGIN = 0.02      # candidates this close to the best distance are ranked by depth


def load_gene_psi(manifests, gene, coverage_threshold):
    """This gene's bulk junction PSI from the cohort junction tables (8A raw TSVs), as
    (psi: sample x junction, depth: sample -> total junction reads). A PSI is kept only where the
    junction's coverage in that sample is >= coverage_threshold."""
    rows = []
    for stype, mpath in manifests.items():
        if not os.path.isfile(mpath):
            continue
        man = pd.read_csv(mpath, sep="\t", dtype=str, keep_default_na=False)
        for rp in man.loc[man["gene"] == gene, "result_path"].tolist():
            rp = str(rp).strip()
            if not rp or rp in ("None", "nan", "NA") or not os.path.isfile(rp):
                continue
            df = pd.read_csv(rp, sep="\t", dtype=str)
            need = {"sample", "phasing", "junction", "rescaled_junction_PSI", "junction_coverage"}
            if not need <= set(df.columns):
                continue
            df = df[df["phasing"] == "bulk"]
            rows.append(pd.DataFrame({"sample": df["sample"].astype(str), "junction": df["junction"],
                                      "psi": pd.to_numeric(df["rescaled_junction_PSI"], errors="coerce"),
                                      "cov": pd.to_numeric(df["junction_coverage"], errors="coerce")}))
    if not rows:
        return pd.DataFrame(), pd.Series(dtype=float)
    long = pd.concat(rows, ignore_index=True).drop_duplicates(["sample", "junction"])
    depth = long.groupby("sample")["cov"].sum(min_count=1)
    ok = long[(long["cov"] >= coverage_threshold) & long["psi"].notna()]
    psi = ok.pivot(index="sample", columns="junction", values="psi")
    return psi, depth


def choose_representative(gene, sample, samples, sample_to_type, sample_bams, psi, depth, exclude):
    """The most typical cohort sample for this gene's splicing: among the candidates (the other
    samples of this sample's type, or, if there are none, the samples of the cohort's most common
    sample type), leaving out any sample with its own outlier junction call in the gene, the one
    whose junction PSIs are closest (mean |PSI - median|) to the candidates' median PSI, over
    junctions covered in both (at least REP_MIN_JUNCTIONS of them, or all the gene has if fewer).
    Among candidates within REP_DEPTH_MARGIN of the best distance, the one with the most junction
    reads, so the sashimi is well populated. Returns dict(sample, sample_type, dist, n) or None."""
    stype = sample_to_type.get(sample)
    usable = [s for s in samples if s != sample and s in sample_bams and s not in exclude]
    cands = [s for s in usable if sample_to_type.get(s) == stype]
    rtype = stype
    if not cands:
        counts = defaultdict(int)
        for s in usable:
            counts[sample_to_type.get(s)] += 1
        if not counts:
            return None
        order = {t: i for i, t in enumerate(dict.fromkeys(sample_to_type.get(s) for s in samples))}
        rtype = max(counts, key=lambda t: (counts[t], -order.get(t, 0)))
        cands = [s for s in usable if sample_to_type.get(s) == rtype]
    if psi is None or psi.empty:
        return None
    cands = [s for s in cands if s in psi.index and psi.loc[s].notna().any()]
    if not cands:
        return None
    sub = psi.loc[cands].dropna(axis=1, how="all")
    median = sub.median(axis=0, skipna=True)
    min_shared = min(REP_MIN_JUNCTIONS, sub.shape[1])
    scored = []
    for c in cands:
        diff = (sub.loc[c] - median).abs().dropna()
        if len(diff) >= min_shared:
            scored.append((float(diff.mean()), c))
    if not scored:
        return None
    best_d = min(d for d, _ in scored)
    close = [(d, c) for d, c in scored if d <= best_d + REP_DEPTH_MARGIN]
    d, best = max(close, key=lambda dc: (float(depth.get(dc[1], 0) or 0), -dc[0]))
    return dict(sample=best, sample_type=rtype, dist=d, n=len(cands))


def sashimi_data(bam, chrom, start, end, keep, min_frac):
    """Coverage (matches + deletions, 0-based start..end) and the junctions inside the window with
    >= min_frac x max coverage reads. Returns (coverage, {junction: count}, n_reads)."""
    L = end - start
    cov = np.zeros(L, dtype=np.int64)
    jx = defaultdict(int)
    n = 0
    for r in bam.fetch(chrom, start, end):
        if not keep(r) or r.cigartuples is None:
            continue
        n += 1
        pos = r.reference_start
        for op, ln in r.cigartuples:
            if op in (0, 2, 7, 8):
                lo, hi = max(pos, start) - start, min(pos + ln, end) - start
                if lo < hi:
                    cov[lo:hi] += 1
                pos += ln
            elif op == 3:
                ss1, ss2 = pos + 1, pos + ln             # 1-based first / last intron base
                if ss1 > start and ss2 <= end:
                    jx[f"{chrom}_{ss1}_{ss2}"] += 1
                pos += ln
    thr = max(1, int(np.floor(min_frac * cov.max()))) if L else 1
    return cov, {j: c for j, c in jx.items() if c >= thr}, n


def intron_transform(start, end, exons, scale):
    """Monotonic genomic -> display map (0-based): exonic bases of the canonical transcript keep
    width 1, every other base gets width `scale`."""
    if scale == 1.0 or not exons:
        return lambda pos: np.asarray(pos, dtype=np.float64)
    merged = []
    for a, b in sorted((max(a, start), min(b, end)) for a, b in exons):
        if a >= b:
            continue
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    gx, dx, cur, d = [start], [float(start)], start, float(start)
    for a, b in merged:
        if cur < a:
            d += (a - cur) * scale
            gx.append(a); dx.append(d)
        d += b - a
        gx.append(b); dx.append(d)
        cur = b
    if cur < end:
        d += (end - cur) * scale
        gx.append(end); dx.append(d)
    gx, dx = np.array(gx, dtype=np.float64), np.array(dx, dtype=np.float64)
    return lambda pos: np.interp(np.asarray(pos, dtype=np.float64), gx, dx)


def _arc_sides(tracks):
    """Junction -> 1 (above the axis) / 0 (below): by mean usage across tracks, a junction goes
    below if it overlaps one already placed above."""
    usage = defaultdict(float)
    for t in tracks:
        for j, c in t["junctions"].items():
            usage[j] += c / max(t["n"], 1) / len(tracks)
    sides = {}
    for j in sorted(usage, key=usage.get, reverse=True):
        _, a, b = _jxn_coords(j)
        above = [(x, y) for k, (x, y) in ((k, _jxn_coords(k)[1:]) for k in sides) if sides[k] == 1]
        sides[j] = 0 if any(a < y and b > x for x, y in above) else 1
    return sides


def draw_sashimi(fig, left, width, top, bottom, tracks, tx, chrom, start, end, args):
    """tracks: list of dict(label, bam, color). Sashimi tracks stacked from `top` down to `bottom`
    (inches), the canonical transcript (intron-scaled) under the last one."""
    keep = read_filter(tx)
    start = max(0, start)
    for t in tracks:
        t["cov"], t["junctions"], t["n"] = sashimi_data(t["bam"], chrom, start, end, keep, args.sashimi_min_frac)
    sides = _arc_sides(tracks)
    xf = intron_transform(start, end, tx["exons"] if tx else [], args.intron_scale)
    x = xf(np.arange(start, end) + 0.5)
    lab_h, ann_h, pad = 0.15, 0.20, 0.04
    tr_h = (top - bottom - ann_h - len(tracks) * (lab_h + pad)) / len(tracks)
    ax_left, ax_w = left + 0.30, width - 0.30
    y = top
    for t in tracks:
        fig.text(ax_left / PAGE_W, (y - lab_h + 0.02) / PAGE_H, t["label"], fontsize=5.5 * FS,
                 ha="left", va="bottom", color="#222")
        y -= lab_h
        ax = _ax(fig, ax_left, y - tr_h, ax_w, tr_h)
        y -= tr_h + pad
        cov = t["cov"].astype(float)
        ymax = max(cov.max(), 1.0)
        ax.fill_between(x, 0, cov, color=t["color"], alpha=0.7, linewidth=0, rasterized=True)
        maxj = max(t["junctions"].values()) if t["junctions"] else 1
        for j, c in t["junctions"].items():
            _, s1, s2 = _jxn_coords(j)
            i1, i2 = s1 - 2 - start, s2 - start              # last exon base before / first after
            y1 = cov[i1] if 0 <= i1 < len(cov) else 0.0
            y2 = cov[i2] if 0 <= i2 < len(cov) else 0.0
            x1, x2 = float(xf(s1 - 2 + 0.5)), float(xf(s2 + 0.5))
            lw = 2.0 * c / maxj + 0.3
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
            top_y = max(y1, y2)
            if sides.get(j, 1) == 1:
                piv = top_y * 1.25
                pts = [(x1, y1), (x1, piv + (y2 - y1) * 0.1), (x2, piv + (y1 - y2) * 0.1), (x2, y2)]
                va, off = "bottom", 1.5
            else:
                piv = -ymax * 0.35
                pts = [(x1, 0), (x1, piv), (x2, piv), (x2, 0)]
                va, off = "top", -1.5
            ax.add_patch(PathPatch(Path(pts, codes), facecolor="none", edgecolor=t["color"], lw=lw, alpha=0.85))
            peak = 0.125 * pts[0][1] + 0.375 * pts[1][1] + 0.375 * pts[2][1] + 0.125 * pts[3][1]   # curve at t = 0.5
            ax.annotate(str(c), ((x1 + x2) / 2, peak), xytext=(0, off + (lw / 2 if off > 0 else -lw / 2)),
                        textcoords="offset points", ha="center", va=va, fontsize=4.5 * FS, annotation_clip=False)
        ax.set_xlim(float(xf(start)), float(xf(end)))
        ax.set_ylim(-ymax * 0.4, ymax * 1.55)
        ax.set_xticks([])
        ax.set_yticks([ymax])
        ax.set_yticklabels([f"{int(ymax):,}"], fontsize=4.5 * FS)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.plot([float(xf(start))] * 2, [0, ymax], color="black", lw=0.6, clip_on=False)
        ax.axhline(0, color="black", lw=0.4)
    ax = _ax(fig, ax_left, bottom, ax_w, ann_h)
    draw_transcript(ax, tx, start, end, xform=xf)


# ----------------------------------------------------------------------------------------------
# QC page
# ----------------------------------------------------------------------------------------------
QC_SAMPLE_COLOR = "#333333"   # reset per sample to a darker shade of its type colour
QC_NAN_COLOR = "#d9d9d9"
# every box plot in the pipeline: plain black outline, no fill
_BLACK_LINE = dict(color="black", linewidth=0.8)
BOX_STYLE = dict(boxprops=_BLACK_LINE, whiskerprops=_BLACK_LINE, capprops=_BLACK_LINE, medianprops=_BLACK_LINE)
# both heatmaps: red (0) - white (0.5) - blue (1), the pipeline's continuous scale (palette.py)
QC_HEAT_CMAP = RED_WHITE_BLUE


def darker(color, f=0.5):
    """`color` mixed with black (f = fraction of black)."""
    r, g, b = to_rgb(color)
    return to_hex((r * (1 - f), g * (1 - f), b * (1 - f)))


def _read_tsv(path, **kw):
    if not path or not os.path.isfile(path):
        if path:
            print(f"  [WARNING] QC input not found: {path}; panel skipped.")
        return None
    return pd.read_csv(path, sep="\t", **kw)


def _qc_spines(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _qc_unavailable(ax, title, msg="not available"):
    ax.set_title(title, fontsize=7 * FS, fontweight="bold")
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, color="#777")
    ax.set_xticks([]); ax.set_yticks([])
    _qc_spines(ax)


def _xlabels(ax, labels, n_rot=3):
    """Tick labels; rotated when there are many."""
    if len(labels) > n_rot:
        ax.set_xticklabels(labels, fontsize=5 * FS, rotation=35, ha="right")
    else:
        ax.set_xticklabels(labels, fontsize=5.5 * FS)


def qc_read_counts(ax, df, sample, shown, tgroups):
    """Stacked on-/off-target bar for the sample and the median of each sample type
    (median mapped reads x median on-target rate), in the types' colours."""
    title = "Read counts"
    if df is None:
        return _qc_unavailable(ax, title)
    df = df.copy()
    df["sample"] = df["sample"].astype(str)
    for c in ("mapped", "on_target_rate"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.set_index("sample")
    if sample not in df.index:
        return _qc_unavailable(ax, title, f"{shown} not in the on-target table")
    bars = [(shown, df.at[sample, "mapped"], df.at[sample, "on_target_rate"] / 100, QC_SAMPLE_COLOR)]
    for t, members, col in tgroups:
        sub = df.reindex([s for s in members if s in df.index])
        if len(sub):
            bars.append((f"{t}\nmedian", sub["mapped"].median(), sub["on_target_rate"].median() / 100, col))
    x = np.arange(len(bars))
    mapped = np.array([b[1] for b in bars], dtype=float)
    rate = np.array([b[2] for b in bars], dtype=float)
    scale, unit = (1e6, "M") if np.nanmax(mapped) >= 1e6 else (1e3, "k")
    on, off = mapped * rate / scale, mapped * (1 - rate) / scale
    cols = [b[3] for b in bars]
    # on-target solid, off-target white, both outlined in the bar's colour
    ax.bar(x, on, color=cols, width=0.65, edgecolor=cols, linewidth=0.8)
    ax.bar(x, off, bottom=on, color="white", width=0.65, edgecolor=cols, linewidth=0.8)
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor="#777777", edgecolor="#777777", linewidth=0.8, label="on-target"),
                      Patch(facecolor="white", edgecolor="#777777", linewidth=0.8, label="off-target")]
    for i in range(len(bars)):
        ax.text(x[i], mapped[i] / scale * 1.02, f"{mapped[i] / scale:.1f}{unit}\n{rate[i]:.0%}",
                ha="center", va="bottom", fontsize=(5.5 if len(bars) <= 4 else 4.5) * FS)
    ax.set_xticks(x)
    _xlabels(ax, [b[0] for b in bars])
    ax.set_ylabel(f"mapped reads ({unit})")
    ax.set_ylim(0, np.nanmax(mapped) / scale * 1.4)
    ax.set_title(title, fontsize=7 * FS, fontweight="bold")
    ax.legend(handles=legend_handles, fontsize=5 * FS, frameon=False, loc="upper right", ncol=2,
              bbox_to_anchor=(1, 1.03))
    _qc_spines(ax)

def qc_read_lengths(ax, df, sample, shown, tgroups):
    """On-/off-target read-length boxes: the sample, and a box per sample type built from the
    median Q1 / median / Q3 (and min / max) of its samples, filled in the types' colours."""
    title = "Read length"
    if df is None:
        return _qc_unavailable(ax, title)
    df = df.copy()
    df["sample"] = df["sample"].astype(str)
    for c in ("min", "q1", "median", "q3", "max"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    n_box = 1 + len(tgroups)
    stats, pos, cols, ticks, groups = [], [], [], [], []
    for g, tt in enumerate(("on_target", "off_target")):
        sub = df[df["targetType"] == tt].set_index("sample")
        rows = [(shown, sub.loc[[sample]] if sample in sub.index else sub.iloc[0:0], QC_SAMPLE_COLOR)]
        rows += [(t, sub.reindex([s for s in members if s in sub.index]), col) for t, members, col in tgroups]
        x0 = g * (n_box + 1)
        for k, (lab, part, col) in enumerate(rows):
            ticks.append((x0 + k, lab))
            if part.empty:
                continue
            q = part[["min", "q1", "median", "q3", "max"]].median()
            iqr = q["q3"] - q["q1"]
            stats.append(dict(med=q["median"], q1=q["q1"], q3=q["q3"], fliers=[],
                              whislo=max(q["min"], q["q1"] - 1.5 * iqr),
                              whishi=min(q["max"], q["q3"] + 1.5 * iqr)))
            pos.append(x0 + k)
            cols.append(col)
        groups.append((x0 + (n_box - 1) / 2, tt.replace("_", "-")))
    if not stats:
        return _qc_unavailable(ax, title, f"{shown} not in the read-length table")
    # read-length boxes keep a colour fill (sample / sample types), black outline
    b = ax.bxp(stats, positions=pos, widths=0.6, showfliers=False, patch_artist=True,
               boxprops=dict(edgecolor="black", linewidth=0.8), whiskerprops=_BLACK_LINE, capprops=_BLACK_LINE,
               medianprops=_BLACK_LINE)
    for patch, c in zip(b["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("black")
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], fontsize=(5 if n_box <= 3 else 4.2) * FS, rotation=35, ha="right")
    ax.set_xlim(-0.6, 2 * n_box + 0.6)
    for xc, lab in groups:
        ax.text(xc, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=6 * FS, fontweight="bold")
    ax.set_ylabel("read length (bp)")
    ax.set_title(title, fontsize=7 * FS, fontweight="bold", pad=14 * FS)
    _qc_spines(ax)

def qc_pca(ax, df, sample, shown):
    """Sample-type validation PCA (PC1 vs PC2) with this sample as a star."""
    title = "Sample type"
    if df is None or not {"PC1", "PC2"} <= set(df.columns):
        return _qc_unavailable(ax, title)
    df = df.copy()
    df["sample"] = df["sample"].astype(str)
    handles = []
    for kind, kw in (("reference", dict(s=6, alpha=0.25, linewidths=0, marker="o", zorder=2)),
                     ("query", dict(s=34, alpha=0.9, marker="D", edgecolors="white", linewidths=0.3, zorder=3))):
        sub = df[df["kind"] == kind]
        for grp, gsub in sub.groupby("group", sort=False):
            col = gsub["color"].iloc[0]
            name = textwrap.fill(str(grp), 18) + f" (n={len(gsub)})"
            h = ax.scatter(gsub["PC1"], gsub["PC2"], color=col, label=name, **kw)
            handles.append(h)
    me = df[(df["kind"] == "query") & (df["sample"] == sample)]
    if not me.empty:
        ax.scatter(me["PC1"], me["PC2"], s=140, marker="*", color=QC_SAMPLE_COLOR,
                   edgecolors="black", linewidths=0.8, zorder=5)
    else:
        ax.text(0.02, 0.98, f"{shown} not in the PCA", transform=ax.transAxes, va="top", color="#b00",
                fontsize=5.5 * FS)
    v1 = df["PC1_var_pct"].iloc[0] if "PC1_var_pct" in df else np.nan
    v2 = df["PC2_var_pct"].iloc[0] if "PC2_var_pct" in df else np.nan
    ax.set_xlabel(f"PC1 ({v1:.1f}%)"); ax.set_ylabel(f"PC2 ({v2:.1f}%)")
    ax.set_title(title, fontsize=7 * FS, fontweight="bold")
    # shape key only (colours are in the page's colour key)
    shape_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color="#888888", markersize=4, alpha=0.6, label="GTEx"),
        plt.Line2D([0], [0], marker="D", linestyle="", color="#888888", markersize=5, label="TEQUILA-Dx"),
    ]
    ax.legend(handles=shape_handles, fontsize=4.5 * FS, frameon=False, loc="best", handletextpad=0.3,
              borderaxespad=0.2, labelspacing=0.3)
    _qc_spines(ax)
    refs = df[df["kind"] == "reference"]
    return [(str(g), gs["color"].iloc[0], len(gs)) for g, gs in refs.groupby("group", sort=False)]


def _qc_box(ax, groups, title, ylabel=None, star=None, ylim=(0, 1), points=False):
    """groups: list of (label, values, color). Tukey boxes (whiskers 1.5 IQR, fliers hidden);
    star = (x index, value) marks this sample; points = also draw each value as a dot;
    ylim = None scales the axis to the data."""
    data = [(lab, np.asarray(v, dtype=float), c) for lab, v, c in groups]
    data = [(lab, v[~np.isnan(v)], c) for lab, v, c in data]
    pos = [i for i, (_, v, _) in enumerate(data) if len(v)]
    if points:      # every value as a dot, under the (unfilled) box
        rng = np.random.default_rng(0)
        for i in pos:
            v = data[i][1]
            many = len(v) > 50
            ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=3 if many else 6, color=data[i][2],
                       alpha=0.35 if many else 0.8, linewidths=0, zorder=1)
    if pos:
        ax.boxplot([data[i][1] for i in pos], positions=pos, widths=0.6, showfliers=False, zorder=2,
                   **BOX_STYLE)
    if star is not None and not np.isnan(star[1]):
        ax.scatter([star[0]], [star[1]], s=110, marker="*", color=QC_SAMPLE_COLOR, edgecolors="black",
                   linewidths=0.5, zorder=5)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels([d[0] for d in data], fontsize=5 * FS, rotation=40, ha="right")
    ax.set_xlim(-0.6, len(data) - 0.4)
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        vals = np.concatenate([d[1] for d in data] + ([np.array([star[1]])] if star is not None else []))
        vals = vals[~np.isnan(vals)]
        if len(vals):
            lo, hi = vals.min(), vals.max()
            padv = max((hi - lo) * 0.15, 0.02)
            ax.set_ylim(max(0, lo - padv), min(1.0, hi + padv) if hi <= 1 else hi + padv)
    ax.tick_params(axis="y", labelsize=5 * FS)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=5 * FS)
    ax.set_title(title, fontsize=6 * FS, fontweight="bold")
    _qc_spines(ax)


def draw_expression_box(ax, gene, sample, stype, tgroups, expr, zmat, threshold):
    """This gene's by-assignment MOTR expression, log2(MOTR + 1): one box per sample type (dots
    behind; the sample left out of its own type) with the sample as a star, and the sample's
    z-score within its type."""
    title = "Expression"
    vals_by_type = {t: m.loc[gene] for t, m in expr.items() if m is not None and gene in m.index}
    own_val = vals_by_type.get(stype, pd.Series(dtype=float)).get(sample, np.nan)
    if not vals_by_type or np.isnan(own_val):
        return _qc_unavailable(ax, title, "no expression data")
    groups = []
    for t, members, col in tgroups:
        v = vals_by_type.get(t)
        if v is not None:
            groups.append((t, np.log2(v.reindex([m_ for m_ in members if m_ in v.index]).astype(float) + 1), col))
    own = next((k for k, (t, *_r) in enumerate(groups) if t == stype), None)
    if own is None:
        groups.insert(0, (stype, pd.Series(dtype=float), QC_SAMPLE_COLOR))
        own = 0
    _qc_box(ax, groups, title, ylabel="log2(MOTR + 1)", star=(own, float(np.log2(own_val + 1))),
            ylim=None, points=True)
    z = zmat.loc[gene, sample] if zmat is not None and gene in zmat.index and sample in zmat.columns else np.nan
    if not np.isnan(z):
        flag = abs(z) >= abs(threshold)
        ax.text(0.5, 1.0, f"z = {z:+.1f} vs {stype}", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=5 * FS, color=QC_SAMPLE_COLOR if flag else "#555",
                fontweight="bold" if flag else "normal")
        ax.set_title(title, fontsize=6 * FS, fontweight="bold", pad=10 * FS)


def _ranked(values, grid):
    """Gene values sorted high -> low, as a curve over % of genes (0-100), sampled on `grid`."""
    v = np.sort(np.asarray(values, dtype=float)[~np.isnan(np.asarray(values, dtype=float))])[::-1]
    if len(v) == 0:
        return np.full(len(grid), np.nan)
    if len(v) == 1:
        return np.full(len(grid), v[0])
    return np.interp(grid, np.linspace(0, 100, len(v)), v)


def qc_curve(ax, matrix, sample, shown, tgroups, ylim=(-0.02, 1.02), ylabel=None):
    """Ranked-gene curves: for each sample type (own type first) the 25th-75th percentile ribbon
    and median line of its samples' curves (thin individual lines if it has < 3 samples); this
    sample as a bold line. The sample is not part of its type's ribbon."""
    grid = np.linspace(0, 100, 201)
    handles = []
    for t, members, col in tgroups:
        members = [s for s in members if s in matrix.columns]
        if not members:
            continue
        curves = np.vstack([_ranked(matrix[s].to_numpy(), grid) for s in members])
        if len(members) >= 3:
            lo, med, hi = np.nanpercentile(curves, [25, 50, 75], axis=0)
            ax.fill_between(grid, lo, hi, color=col, alpha=0.25, linewidth=0, zorder=1)
            ax.plot(grid, med, color=col, linewidth=0.8, zorder=2)
        else:
            for c in curves:
                ax.plot(grid, c, color=col, linewidth=0.6, alpha=0.7, zorder=2)
    ax.plot(grid, _ranked(matrix[sample].to_numpy(), grid), color=QC_SAMPLE_COLOR, linewidth=1.6, zorder=3)
    ax.set_xlim(0, 100)
    ax.set_ylim(*ylim)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=5 * FS)
    ax.set_xlabel("% of genes (ranked)", fontsize=5 * FS)
    ax.tick_params(labelsize=5 * FS)
    ax.set_title("ranked genes", fontsize=6 * FS, fontweight="bold")
    _qc_spines(ax)


GRID_FONT_MAX, GRID_FONT_MIN = 6.5, 4.0     # gene-label font (pt): shrinks with gene count, then pages are added


def _grid_layout(n_genes, n_cells, max_label, width, h_first, h_next):
    """Font size and grid geometry for the labelled gene grid: the largest font (<= GRID_FONT_MAX)
    at which every gene fits on the first page, else GRID_FONT_MIN with extra pages."""
    def geom(f):
        row_h = f * 1.25 / 72                            # inches
        cell_w = f * 1.7 / 72
        label_w = max_label * 0.62 * f / 72
        col_w = label_w + 0.05 + n_cells * cell_w + 0.16
        ncols = max(1, int(width // col_w))
        return row_h, cell_w, label_w, col_w, ncols, max(1, int(h_first // row_h)), max(1, int(h_next // row_h))
    f = GRID_FONT_MAX
    while f > GRID_FONT_MIN:
        g = geom(f)
        if g[4] * g[5] >= n_genes:
            return (f,) + g
        f = round(f - 0.25, 2)
    return (GRID_FONT_MIN,) + geom(GRID_FONT_MIN)


def _metric_header(fig, shown, title, sub):
    nav_links(fig)
    fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H, f"{shown}  |  {title}", fontsize=11 * FS, fontweight="bold", va="top")
    fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H, sub, fontsize=6.5 * FS, va="top", color="#333")


def _grid_pages(pdf, names, values, cell_cols, cmap, norm, shown, title, sub, summary_fn, cbar_label,
                outline=None, gap_after_first=False, cbar_ticks=None):
    """Labelled gene grid over as many pages as needed: each gene is a row (its name, then one cell
    per column of `values`), rows stacked into page columns. cell_cols: header tick colour per cell;
    outline: genes x cells booleans for cells to outline; gap_after_first: a small gap after the
    first cell. The first page also gets `summary_fn(fig)` above the grid."""
    n_cells = values.shape[1]
    gap = 0.45 if gap_after_first else 0.0          # in cell widths
    cbar_space = 0.6                                 # room right of the grid for the colour bar
    grid_left, grid_w = MARGIN_L, PAGE_W - MARGIN_L - MARGIN_R - cbar_space
    first_top, next_top, bottom = 5.3, PAGE_H - MARGIN_T - 0.05, 0.3
    f, row_h, cell_w, label_w, col_w, ncols, rows1, rows2 = _grid_layout(
        len(names), n_cells + gap, max(len(n) for n in names), grid_w, first_top - bottom, next_top - bottom)

    def cell_x(x0, j):
        return x0 + label_w + 0.05 + (j + (gap if (gap_after_first and j > 0) else 0)) * cell_w

    def grid(fig, start, top, nrows):
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, PAGE_W); ax.set_ylim(0, PAGE_H); ax.axis("off")
        end = min(len(names), start + nrows * ncols)
        rects, colors, nan_rects, out_rects = [], [], [], []
        for k in range(start, end):
            c, r = divmod(k - start, nrows)
            x0 = grid_left + c * col_w
            y = top - (r + 1) * row_h
            ax.text(x0 + label_w, y + row_h / 2, names[k], ha="right", va="center", fontsize=f)
            for j, v in enumerate(values[k]):
                rx, ry, rw, rh = cell_x(x0, j), y + row_h * 0.06, cell_w * 0.92, row_h * 0.88
                if np.isnan(v):
                    nan_rects.append(Rectangle((rx, ry), rw, rh))
                else:
                    rects.append(Rectangle((rx, ry), rw, rh)); colors.append(cmap(norm(float(v))))
                    if outline is not None and outline[k, j]:
                        out_rects.append(Rectangle((rx, ry), rw, rh))
        if rects:
            ax.add_collection(PatchCollection(rects, facecolors=colors, edgecolors="none"))
        if nan_rects:
            ax.add_collection(PatchCollection(nan_rects, facecolors=QC_NAN_COLOR, edgecolors="white",
                                              hatch="//////", linewidths=0))
        if out_rects:
            ax.add_collection(PatchCollection(out_rects, facecolors="none", edgecolors="black", linewidths=0.4))
        # coloured ticks above each page column's cells: which cell is which
        used_cols = (end - start + nrows - 1) // nrows
        for c in range(used_cols):
            x0 = grid_left + c * col_w
            for j, col in enumerate(cell_cols):
                ax.add_patch(Rectangle((cell_x(x0, j), top + 0.02), cell_w * 0.92, 0.035, facecolor=col,
                                       edgecolor="none"))
        cax = _ax(fig, PAGE_W - MARGIN_R - cbar_space + 0.2, top - 1.6, 0.09, 1.6)
        cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
        cb.outline.set_visible(False)
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)
        cb.ax.tick_params(labelsize=5 * FS)
        cb.set_label(cbar_label, fontsize=5 * FS)
        return end

    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    _metric_header(fig, shown, title, sub)
    summary_fn(fig)
    done = grid(fig, 0, first_top, rows1)
    page = 1
    while True:
        pdf.savefig(fig, dpi=RASTER_DPI)
        plt.close(fig)
        if done >= len(names):
            break
        page += 1
        fig = plt.figure(figsize=(PAGE_W, PAGE_H))
        _metric_header(fig, shown, f"{title} (continued, page {page})", sub)
        done = grid(fig, done, next_top, rows2)


def qc_mean_panel(ax, matrix, sample, stype, tgroups, title):
    """Box plot of every sample's mean over genes, one box per sample type (dots behind), with this
    sample as a star among its own type's dots; 0-1 axis."""
    if matrix is None or matrix.empty or sample not in matrix.columns:
        return _qc_unavailable(ax, title)
    per_sample = matrix.mean(axis=0, skipna=True)
    groups = [(t, per_sample.reindex([c for c in m_ if c in matrix.columns]), col) for t, m_, col in tgroups]
    groups = [g for g in groups if len(g[1])]
    own = next((k for k, (t, *_r) in enumerate(groups) if t == stype), None)
    if own is None:
        groups.insert(0, (stype, pd.Series(dtype=float), QC_SAMPLE_COLOR))
        own = 0
    _qc_box(ax, groups, title, star=(own, float(per_sample.get(sample, np.nan))), points=True)


def draw_metric_pages(pdf, matrix, sample, shown, stype, tgroups, title, desc, value_label):
    """One QC metric (gene x sample matrix of 0-1 values): a summary strip (ranked-gene curves,
    per-gene box plot, per-sample-mean box plot) and a grid of every gene -- its name, then a cell
    for the sample and one per sample type's median -- continued on further pages if needed.
    Genes are sorted high to low by the sample's value, then the median of every other sample, then name."""
    if matrix is None or matrix.empty or sample not in matrix.columns:
        print(f"  {title}: no data for {shown}; page skipped.")
        return
    others = [s_ for s_ in matrix.columns if s_ != sample]
    rest_med = matrix[others].median(axis=1, skipna=True) if others else pd.Series(np.nan, index=matrix.index)
    # high to low: the sample's value, then the median of every other sample (NaN last), then name
    order = (pd.DataFrame({"s": matrix[sample], "c": rest_med, "g": matrix.index.astype(str)}, index=matrix.index)
             .sort_values(["s", "c", "g"], ascending=[False, False, True], na_position="last",
                          kind="mergesort").index)
    meds = []
    for t, members, col in tgroups:
        cols_ = [s_ for s_ in members if s_ in matrix.columns]
        if cols_:
            meds.append((t, cols_, col, matrix[cols_].median(axis=1, skipna=True)))
    cells = [(shown, QC_SAMPLE_COLOR, matrix[sample])] + [(f"{t} median", col, m_) for t, _, col, m_ in meds]
    values = np.vstack([v.reindex(order).to_numpy(dtype=float) for _, _, v in cells]).T   # genes x cells

    def summary(fig):
        kax = _ax(fig, 0.95, 6.15, 2.4, 1.05)
        qc_curve(kax, matrix, sample, shown, tgroups)
        bax1 = _ax(fig, 3.85, 6.15, 1.45, 1.05)
        _qc_box(bax1, [(shown, matrix[sample], QC_SAMPLE_COLOR)] + [(t, m_, col) for t, _, col, m_ in meds],
                "per gene (type = median)", points=True)
        qc_mean_panel(_ax(fig, 5.95, 6.15, 1.45, 1.05), matrix, sample, stype, tgroups, "per-sample mean")
        draw_page_legend(fig, 8.2, 7.45, shown, tgroups)

    _grid_pages(pdf, [str(g) for g in order], values, [c for _, c, _ in cells], QC_HEAT_CMAP,
                plt.Normalize(0, 1), shown, title, desc, summary, value_label)


Z_CAP = 5.0          # z-score colour scale runs -Z_CAP .. +Z_CAP


def _cluster_order(df):
    """Column order from average-linkage clustering on the columns' values (NaN -> 0)."""
    if df.shape[1] < 3:
        return list(df.columns)
    from scipy.cluster.hierarchy import linkage, leaves_list
    X = df.fillna(0).to_numpy(dtype=float).T
    return [df.columns[i] for i in leaves_list(linkage(X, method="average", metric="euclidean"))]


def draw_zscore_pages(pdf, zmat, sample, shown, stype, tgroups, norm_name, threshold):
    """Gene-expression z-scores (by splice-site assignment, one normalisation) for this sample and
    every other sample of its type: a summary strip and a grid of every gene -- its name, this
    sample's cell, then the other samples clustered by similarity. Genes sorted by this sample's
    z-score, lowest first. Cells at or beyond +/- |threshold| are outlined."""
    title = f"Gene expression ({norm_name} z-score)"
    if zmat is None or zmat.empty or sample not in zmat.columns:
        print(f"  {title}: no data for {shown}; page skipped.")
        return
    others = _cluster_order(zmat[[c for c in zmat.columns if c != sample]])
    cols = [sample] + others
    order = (pd.DataFrame({"s": zmat[sample], "g": zmat.index.astype(str)}, index=zmat.index)
             .sort_values(["s", "g"], na_position="last", kind="mergesort").index)
    values = zmat.loc[order, cols].to_numpy(dtype=float)
    thr = abs(float(threshold))
    outline = np.abs(np.nan_to_num(values, nan=0.0)) >= thr
    tcol = next((c for t, _, c in tgroups if t == stype), "#999999")

    def summary(fig):
        kax = _ax(fig, 0.95, 6.15, 3.0, 1.05)
        own = [(t, [c for c in m_ if c in zmat.columns], col) for t, m_, col in tgroups if t == stype]
        qc_curve(kax, zmat.clip(-Z_CAP, Z_CAP), sample, shown, own, ylim=(-Z_CAP - 0.3, Z_CAP + 0.3),
                 ylabel="z-score")
        kax.axhline(-thr, color="#777", linewidth=0.6, linestyle=":")
        kax.axhline(thr, color="#777", linewidth=0.6, linestyle=":")
        draw_page_legend(fig, 8.2, 7.45, shown, [g for g in tgroups if g[0] == stype])

    _grid_pages(pdf, [str(g) for g in order], values, [QC_SAMPLE_COLOR] + [tcol] * len(others),
                QC_HEAT_CMAP, plt.Normalize(-Z_CAP, Z_CAP), shown, title, "", summary,
                f"z-score (capped at +/-{Z_CAP:g})", outline=outline, gap_after_first=True,
                cbar_ticks=[-Z_CAP, -thr, 0, thr, Z_CAP])


def load_flr_matrix(path):
    df = _read_tsv(path, index_col=0)
    if df is None:
        return None
    df.columns = df.columns.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def load_testability_matrix(hits_tsv):
    """gene x sample n_canonical_jxns_tested / n_canonical_jxns (NaN when no canonical junctions)."""
    if not hits_tsv or not os.path.isfile(hits_tsv):
        return None
    cols = {"sample", "gene", "n_canonical_jxns", "n_canonical_jxns_tested"}
    df = pd.read_csv(hits_tsv, sep="\t", dtype=str, keep_default_na=False, usecols=lambda c: c in cols)
    if not cols <= set(df.columns):
        print("  [WARNING] merged_all_hits.tsv has no n_canonical_jxns / n_canonical_jxns_tested columns; "
              "testability panel skipped.")
        return None
    n = pd.to_numeric(df["n_canonical_jxns"], errors="coerce")
    t = pd.to_numeric(df["n_canonical_jxns_tested"], errors="coerce")
    df["frac"] = (t / n).where(n > 0)
    df["sample"] = df["sample"].astype(str)
    return df.pivot_table(index="gene", columns="sample", values="frac", aggfunc="first", dropna=False)


def _hit_flag(v):
    return _nonempty(v) and str(v).strip().lower() not in ("false", "0")


# Combination x tier cells the tier rules can produce (tier_rules.py, shared with merge_hits.py);
# every other cell of the overview grid is drawn grey with a dash
_POSSIBLE_CELLS = tier_rules.possible_cells()


def draw_hit_tier_grid(fig, rows, left=1.3, bottom=0.35):
    """Tier x hit-category grid: columns are the combinations of hit categories (variant, ASE,
    GTEx or cohort outlier junction; upset-style dots underneath), rows are tiers 1-6, each cell
    the number of this sample's candidate genes (printed; zero left blank), shaded white -> the
    sample colour (square-root scale). Bars above / to the right give the column / row totals."""
    cats = ["variant", "ASE", "junction"]          # dot-matrix rows, top to bottom
    order = ["variant", "junction", "ASE"]         # combination order, as in the cohort upset plot
    combos = [frozenset(c) for k in range(1, 4) for c in itertools.combinations(order, k)]
    tiers = [1, 2, 3, 4, 5, 6]
    grid = np.zeros((len(tiers), len(combos)), dtype=int)
    no_cat = 0
    for _, r in rows.iterrows():
        have = frozenset(c for c, ok in (
            ("variant", _hit_flag(r.get("variant"))),
            ("ASE", _hit_flag(r.get("ASE"))),
            ("junction", _hit_flag(r.get("outlier_junction")) or _hit_flag(r.get("cohort_outlier_junction"))),
        ) if ok)
        tier = pd.to_numeric(pd.Series([r.get("tier")]), errors="coerce").iloc[0]
        if pd.isna(tier) or int(tier) not in tiers:
            continue
        if not have:
            no_cat += 1
            continue
        grid[tiers.index(int(tier)), combos.index(have)] += 1
    if no_cat:
        print(f"  [WARNING] {no_cat} candidate gene(s) have a tier but no variant / ASE / junction hit; "
              "left out of the hit-category grid.")
    n_genes = int(grid.sum())

    cw, rh = 0.62, 0.32                             # cell width / height (inches)
    gw, gh = cw * len(combos), rh * len(tiers)
    grid_bottom = bottom + 1.0
    # main grid
    ax = _ax(fig, left, grid_bottom, gw, gh)
    vmax = max(int(grid.max()), 1)
    cmap = LinearSegmentedColormap.from_list("hits", ["#ffffff", QC_SAMPLE_COLOR])
    for i in range(len(tiers)):
        for j in range(len(combos)):
            v = grid[i, j]
            impossible = (combos[j], tiers[i]) not in _POSSIBLE_CELLS
            frac = np.sqrt(v / vmax)
            if impossible and not v:     # cannot happen under the tier rules: light grey with a dash
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="#eeeeee", edgecolor="#dddddd", linewidth=0.5))
                ax.text(j + 0.5, i + 0.5, "\u2013", ha="center", va="center", fontsize=6.5 * FS, color="#aaaaaa")
                continue
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(frac), edgecolor="#dddddd", linewidth=0.5))
            if v:
                ax.text(j + 0.5, i + 0.5, str(v), ha="center", va="center", fontsize=6.5 * FS,
                        color="white" if frac > 0.6 else "black")
            if impossible and v:         # should never happen: the tier rules and the hit flags disagree
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="none", edgecolor="#c23637", linewidth=1.2))
                print(f"  [WARNING] {v} gene(s) in tier {tiers[i]} with hit categories "
                      f"{'+'.join(sorted(combos[j]))}, which the tier rules should not produce.")
    ax.set_xlim(0, len(combos)); ax.set_ylim(len(tiers), 0)
    ax.set_xticks([]); ax.set_yticks(np.arange(len(tiers)) + 0.5)
    ax.set_yticklabels([f"Tier {t}" for t in tiers], fontsize=6 * FS)
    ax.tick_params(axis="y", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # column totals above, row totals to the right
    col_tot, row_tot = grid.sum(axis=0), grid.sum(axis=1)
    axt = _ax(fig, left, grid_bottom + gh + 0.08, gw, 0.7)
    axt.bar(np.arange(len(combos)) + 0.5, col_tot, width=0.6, color=QC_SAMPLE_COLOR)
    for j, v in enumerate(col_tot):
        axt.text(j + 0.5, v, str(v), ha="center", va="bottom", fontsize=5.5 * FS)
    axt.set_xlim(0, len(combos)); axt.set_ylim(0, max(col_tot.max(), 1) * 1.3)
    axt.axis("off")
    axt.set_title(f"Candidate genes by hit category and tier ({n_genes} genes)", fontsize=7 * FS,
                  fontweight="bold", pad=4)
    axr = _ax(fig, left + gw + 0.08, grid_bottom, 0.9, gh)
    axr.barh(np.arange(len(tiers)) + 0.5, row_tot, height=0.6, color=QC_SAMPLE_COLOR)
    for i, v in enumerate(row_tot):
        axr.text(v, i + 0.5, f" {v}", ha="left", va="center", fontsize=5.5 * FS)
    axr.set_ylim(len(tiers), 0); axr.set_xlim(0, max(row_tot.max(), 1) * 1.35)
    axr.axis("off")
    # upset dots under the columns
    axm = _ax(fig, left, bottom, gw, 0.85)
    cat_y = {c: len(cats) - k for k, c in enumerate(cats)}
    for j, c in enumerate(combos):
        ys = [cat_y[k] for k in c]
        if len(ys) > 1:
            axm.plot([j + 0.5, j + 0.5], [min(ys), max(ys)], color="black", linewidth=1.5, zorder=2)
        for k in cats:
            axm.scatter([j + 0.5], [cat_y[k]], s=40, zorder=3, color="black" if k in c else "#dddddd")
    axm.set_xlim(0, len(combos)); axm.set_ylim(0.5, len(cats) + 0.5)
    axm.set_yticks(list(cat_y.values())); axm.set_yticklabels(list(cat_y.keys()), fontsize=6 * FS)
    axm.tick_params(axis="y", length=0)
    axm.set_xticks([])
    for sp in axm.spines.values():
        sp.set_visible(False)


def draw_page_legend(fig, left, top, shown, tgroups, refs=None):
    """Page-wide colour key: the sample (star), each sample type (square) and, if given, the GTEx
    reference tissues of the sample-type PCA (faint dots)."""
    width, height = PAGE_W - MARGIN_R - left, 3.0
    ax = _ax(fig, left, top - height, width, height)
    ax.set_xlim(0, width); ax.set_ylim(top - height, top); ax.axis("off")   # data units = inches
    mx, tx, step = 0.15, 0.28, 0.13
    y = top
    x0 = 0.08                             # inset from the box edge
    ax.text(x0, y - 0.06, "TEQUILA-Dx", fontsize=6 * FS, fontweight="bold", va="top")
    y -= 0.06
    y -= step * 1.5
    entries = [("star", QC_SAMPLE_COLOR, shown)] + [("square", col, f"{t} (n={len(m_)})") for t, m_, col in tgroups]
    if refs:
        entries.append(("head", None, "GTEx"))
        entries += [("dot", col, textwrap.shorten(str(g), 24, placeholder="...")) for g, col, _ in refs]
    for kind, col, lab in entries:
        if kind == "head":
            y -= step * 0.3
            ax.text(0.08, y, lab, fontsize=5 * FS, va="center", color="#555")
        else:
            if kind == "star":
                ax.scatter([mx], [y], marker="*", s=80, color=col, edgecolors="black", linewidths=0.5, clip_on=False)
            elif kind == "square":
                ax.scatter([mx], [y], marker="s", s=35, color=col, clip_on=False)
            else:
                ax.scatter([mx], [y], marker="o", s=14, color=col, alpha=0.5, linewidths=0, clip_on=False)
            ax.text(tx, y, lab, fontsize=5 * FS, va="center")
        y -= step
    # little box around the key, sized to the entries and the longest label
    longest = max(len(e[2]) for e in entries + [("", None, "TEQUILA-Dx")])
    box_w = min(width, tx + longest * 0.62 * 5 * FS / 72 + 0.1)
    ax.add_patch(FancyBboxPatch((0, y + step * 0.4), box_w, top - (y + step * 0.4),
                                boxstyle="round,pad=0,rounding_size=0.05", facecolor="none",
                                edgecolor="#999999", linewidth=0.6, clip_on=False))


def draw_qc_page(pdf, args, sample, shown, sample_to_type, type_colors, hit_rows):
    """Overview page (read counts, read length, sample-type PCA, hit categories, per-sample QC means).
    Returns what draw_qc_sections needs to draw the QC sections after the gene pages."""
    stype = sample_to_type.get(sample, "unknown")
    global QC_SAMPLE_COLOR
    QC_SAMPLE_COLOR = darker(type_colors.get(stype, "#4c72b0"))
    global ARC_COLOR
    ARC_COLOR = QC_SAMPLE_COLOR       # IGV sashimi arcs in the sample's colour
    # comparison groups (one per sample type, own type first, then by size) exclude the sample itself
    all_samples = [str(s) for s in args.samples if str(s) != sample]
    counts = {}
    for s_ in all_samples:
        t = sample_to_type.get(s_)
        if t:
            counts[t] = counts.get(t, 0) + 1
    type_order = [stype] + sorted((t for t in counts if t != stype), key=lambda t: -counts[t])
    tgroups = [(t, [s for s in all_samples if sample_to_type.get(s) == t], type_colors.get(t, "#7f7f7f"))
               for t in type_order]
    tgroups = [g for g in tgroups if g[1]]

    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    nav_links(fig)
    fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H, f"{shown}  |  Overview  (sample type: {stype})",
             fontsize=11 * FS, fontweight="bold", va="top")
    fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H,
             ", ".join(f"{t} (n = {len(m_)}{' others' if t == stype else ''})" for t, m_, _ in tgroups),
             fontsize=6.5 * FS, va="top", color="#333")
    qc_read_counts(_ax(fig, 0.95, 5.0, 2.3, 2.6), _read_tsv(args.on_target_tsv), sample, shown, tgroups)
    qc_read_lengths(_ax(fig, 3.95, 5.0, 2.6, 2.4), _read_tsv(args.read_attributes_tsv), sample, shown, tgroups)
    refs = qc_pca(_ax(fig, 7.0, 5.0, 2.3, 2.6), _read_tsv(args.pca_coords), sample, shown)
    draw_page_legend(fig, 9.55, 7.6, shown, tgroups, refs)
    flr = load_flr_matrix(args.flr_matrix)
    testability = load_testability_matrix(args.hits_tsv)
    coverage = load_flr_matrix(args.coverage_matrix)
    draw_hit_tier_grid(fig, hit_rows, left=1.1, bottom=0.35)
    # sample quality: per-sample means of the three gene-level QC metrics, each on its own 0-1 axis
    # with a y label saying what the metric measures (spaced so each label clears the panel before)
    for k, (mat, name, ylab) in enumerate((
            (coverage, "Coverage",
             f"Fraction of canonical transcript at >= {args.coverage_depth}x\n(average across panel)"),
            (testability, "Junction testability", "Fraction of canonical junctions tested\n(average across panel)"),
            (flr, "Full-length ratio", "Fraction of canonical transcript covered by a single read\n(average across panel)"))):
        ax = _ax(fig, 7.2 + k * 1.4, 1.2, 0.75, 2.9)
        qc_mean_panel(ax, mat, sample, stype, tgroups, name)
        ax.set_title(name, fontsize=5.5 * FS, fontweight="bold")
        ax.set_ylabel(ylab, fontsize=4.8 * FS, labelpad=2, linespacing=1.1)
    pdf.savefig(fig, dpi=RASTER_DPI)
    plt.close(fig)

    return dict(stype=stype, tgroups=tgroups, flr=flr, testability=testability, coverage=coverage)


def draw_qc_sections(pdf, args, sample, shown, ctx):
    """After the gene pages: gene expression (MOTR z-score), canonical-transcript coverage,
    junction testability, full-length ratio."""
    global REVIEW_CTX
    stype, tgroups = ctx["stype"], ctx["tgroups"]
    REVIEW_CTX = ("qc", "Gene expression (MOTR z-score)")
    draw_zscore_pages(pdf, load_flr_matrix(args.zscore_motr), sample, shown, stype, tgroups, "MOTR",
                      args.zscore_threshold)
    for mat, title, desc, label in (
            (ctx["coverage"], "Canonical-transcript coverage",
             f"fraction of each gene's canonical-transcript bases at depth >= {args.coverage_depth}",
             f"fraction >= {args.coverage_depth}x"),
            (ctx["testability"], "Junction testability",
             "fraction of each gene's canonical junctions tested", "fraction tested"),
            (ctx["flr"], "Full-length ratio",
             "fraction of each gene's canonical-transcript bases covered by a single read, on average",
             "full-length ratio")):
        REVIEW_CTX = ("qc", title)
        draw_metric_pages(pdf, mat, sample, shown, stype, tgroups, title, desc, label)
    _flush_pdf_images(pdf)


def main():
    args = parse_args()
    sample = args.sample
    shown_name = sample                # the report always uses sample names, never aliases
    sample_to_type = dict(zip(args.samples, args.sample_types))
    type_colors = {}
    for t, c in zip(args.sample_types, args.colors):
        type_colors.setdefault(t, c)
    manifests = {}
    for item in args.cohort_manifests:
        st, _, path = item.partition(":")
        if path:
            manifests[st] = path

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    hits = read_hit_rows(args.hits_tsv, sep="\t", dtype=str, keep_default_na=False) \
        if os.path.isfile(args.hits_tsv) else pd.DataFrame()
    # (metric, junction) -> every sample on the panel with an outlier call for it (labelled in the
    # distribution panels)
    outliers_by_pair = defaultdict(set)
    for _, hrow in hits.iterrows():
        for c in junction_calls(hrow):
            outliers_by_pair[(c["metric"], c["junction"])].add(str(hrow["sample"]))
    mine = hits[hits["sample"].astype(str) == str(sample)].copy() if not hits.empty else hits
    if not mine.empty:
        mine["_rank"] = pd.to_numeric(mine["ranking"], errors="coerce")
        mine = mine[pd.to_numeric(mine["tier"], errors="coerce") <= args.max_tier]
        mine = mine.sort_values("_rank")
        if args.top_n > 0:
            mine = mine.head(args.top_n)
    genes = list(mine["gene"]) if not mine.empty else []
    print(f"{sample}: {len(genes)} top hit gene(s): {', '.join(genes) if genes else 'none'}")

    regions = load_bed(args.bed)
    tx_by_gene = load_canonical_transcripts(args.gtf, genes) if genes else {}
    hap_bams = load_gene_bams(args.gene_bam_mapping)
    sample_bams = dict(zip(map(str, args.samples), args.sample_bams)) if args.sample_bams else {}
    expr = {}
    for item in args.expression_matrices:
        st, _, path = item.partition(":")
        if path and os.path.isfile(path):
            expr[st] = load_flr_matrix(path)
    zmat_motr = load_flr_matrix(args.zscore_motr) if args.zscore_motr else None
    # gene -> samples with their own outlier junction call there (never used as the representative)
    jxn_outlier_samples = defaultdict(set)
    for _, hrow in hits.iterrows():
        if junction_calls(hrow):
            jxn_outlier_samples[str(hrow["gene"])].add(str(hrow["sample"]))
    type_color = lambda st: type_colors.get(st, "#7f7f7f")
    fasta = pysam.FastaFile(args.genome)
    bulk = pysam.AlignmentFile(args.bam, "rb")

    global REVIEW_CTX
    review_genes = []
    with PdfPages(args.outfile) as real_pdf:
        pdf = CapturePdf(real_pdf, args.html_dpi if args.html_out else 0)
        REVIEW_CTX = ("overview", None)
        qc_ctx = draw_qc_page(pdf, args, str(sample), shown_name, sample_to_type, type_colors,
                              hits[hits["sample"].astype(str) == str(sample)] if not hits.empty else hits)
        if not genes:
            fig = plt.figure(figsize=(PAGE_W, PAGE_H))
            fig.text(0.5, 0.5, f"{shown_name}: no candidate gene hits.", ha="center", va="center", fontsize=14 * FS)
            pdf.savefig(fig)
            plt.close(fig)
            draw_qc_sections(pdf, args, str(sample), shown_name, qc_ctx)
            pages = pdf.pages
            genes_done = True
        else:
            genes_done = False
        if not genes_done:
            pages = _draw_gene_pages(pdf, args, sample, shown_name, mine, regions, tx_by_gene, hap_bams,
                                     sample_bams, expr, zmat_motr, jxn_outlier_samples, type_color, fasta,
                                     bulk, sample_to_type, type_colors, manifests, outliers_by_pair, qc_ctx,
                                     review_genes)
        contents_rows = draw_contents(pdf, shown_name, review_genes)
    bulk.close()
    fasta.close()
    add_pdf_navigation(args.outfile, pdf, contents_rows)
    print(f"Saved: {args.outfile}")
    if args.html_out:
        import review_html
        review_html.write_review_page(args.html_out, shown_name, qc_ctx["stype"], review_genes, pages)
        print(f"Saved: {args.html_out}")


def _draw_gene_pages(pdf, args, sample, shown_name, mine, regions, tx_by_gene, hap_bams, sample_bams, expr,
                     zmat_motr, jxn_outlier_samples, type_color, fasta, bulk, sample_to_type, type_colors,
                     manifests, outliers_by_pair, qc_ctx, review_genes):
    """Every candidate gene's pages, then the QC sections. Returns the captured pages."""
    global REVIEW_CTX
    gene_order = [g for g in mine["gene"] if g in regions]
    if True:

        pager = Pager(pdf, None)
        for _, row in mine.iterrows():
            gene = row["gene"]
            if gene not in regions:
                print(f"  [WARNING] {gene} not in BED; skipped.")
                continue
            chrom, gstart, gend = regions[gene]
            tx = tx_by_gene.get(gene)
            variants = parse_variants(row)
            calls = junction_calls(row)
            events = build_events(calls)
            pairs = metric_junction_pairs(calls)
            cats = [c for c in ("pathogenic_variant", "variant", "ASE", "outlier_junction",
                                "cohort_outlier_junction")
                    if _nonempty(row.get(c)) and str(row.get(c)) not in ("False", "FALSE")]
            gi = gene_order.index(gene)
            prev_g = gene_order[gi - 1] if gi > 0 else None
            next_g = gene_order[gi + 1] if gi + 1 < len(gene_order) else None

            def header(fig, part, gene=gene, row=row, cats=cats, prev_g=prev_g, next_g=next_g):
                nav_links(fig, prev_g, next_g)
                fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H,
                         f"{shown_name}  |  {gene}  |  rank {row.get('ranking', '.')}, tier {row.get('tier', '.')}"
                         + (f"  (page {part})" if part > 1 else ""),
                         fontsize=11 * FS, fontweight="bold", va="top")
                sub = "hits: " + (", ".join(f"{c}={row.get(c)}" for c in cats) or ".")
                if _nonempty(row.get("inheritance_patterns")):
                    inh = str(row.get("inheritance_patterns"))
                    if tier_rules.is_dominant(inh) and str(row.get("haploinsufficient", "")).strip().upper() == "TRUE":
                        inh += " (haploinsufficient)"
                    sub += f"   |   inheritance: {inh}"
                fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H, sub, fontsize=6.5 * FS, va="top", color="#333")

            REVIEW_CTX = ("gene", gene)
            review_genes.append(dict(
                gene=gene, rank=str(row.get("ranking", ".")), tier=str(row.get("tier", ".")),
                hits=", ".join(f"{c}={row.get(c)}" for c in cats), cats=cats,
                inheritance=str(row.get("inheritance_patterns", "")) if _nonempty(row.get("inheritance_patterns")) else "",
                haploinsufficient=str(row.get("haploinsufficient", "")).strip().upper() == "TRUE"))
            pager.new_gene(header)
            print(f"  {gene}: {len(variants)} variant(s), {len(events)} event(s), {len(pairs)} distribution(s)")
            full_left, full_w = MARGIN_L, PAGE_W - MARGIN_L - MARGIN_R
            var_marks = [v["pos"] - 1 for v in variants if v["chrom"] == chrom]
            ev_spans = [(min(_jxn_coords(j)[1] for j in e["junctions"]) - 1,
                         max(_jxn_coords(j)[2] for j in e["junctions"])) for e in events]

            # 1-3. gene body: bulk, then hap1 / hap2 -- haplotypes only when the gene has ASE or a
            # hap1 / hap2 outlier junction (GTEx or cohort)
            tracks = [("bulk", bulk)]
            hap_handles = []
            has_ase = str(row.get("ASE", "")).strip().upper() == "TRUE"
            has_hap_jxn = any(c["phasing"] in ("hap1", "hap2") for c in calls)
            if gene in hap_bams and (has_ase or has_hap_jxn):
                for name, path in zip(("hap1", "hap2"), hap_bams[gene]):
                    hh = pysam.AlignmentFile(path, "rb")
                    hap_handles.append(hh)
                    tracks.append((name, hh))
            # the canonical transcript is drawn once, under the last (bottom) track. The IGV views
            # take GENE_IGV_FRAC of the width; beside the bulk view, the gene's expression box plot
            igv_gene_w = full_w * GENE_IGV_FRAC
            for k, (name, handle) in enumerate(tracks):
                fig, top, h = pager.next_row()
                if k == 0:
                    bx_left = full_left + igv_gene_w + 0.75
                    bx = _ax(fig, bx_left, top - h + 0.75, full_left + full_w - bx_left - 0.05, h - 1.15)
                    draw_expression_box(bx, gene, str(sample), qc_ctx["stype"], qc_ctx["tgroups"], expr,
                                        zmat_motr, args.zscore_threshold)
                alignment_panel(fig, full_left, igv_gene_w, top, h,
                                "",
                                handle, chrom, gstart, gend, fasta, tx, args, args.seed,
                                marks=var_marks, spans=ev_spans, track_label=name,
                                show_tx=(k == len(tracks) - 1), keep_indels=indel_specs(variants))
            for hh in hap_handles:
                hh.close()

            # 4. variant zooms, up to 3 per row
            for i in range(0, len(variants), VARIANT_PANELS_PER_ROW):
                chunk = variants[i:i + VARIANT_PANELS_PER_ROW]
                fig, top, h = pager.next_row()
                for k, ((left, w), v) in enumerate(zip(panel_columns(len(chunk), VARIANT_PANELS_PER_ROW), chunk)):
                    vs = v["pos"] - 1
                    ve = vs + max(len(v["ref"]), 1)
                    title = v["id"] + (f"  {v['clnsig']}" if v["clnsig"] else "")
                    alignment_panel(fig, left, w, top, h, title, bulk, v["chrom"],
                                    vs - args.variant_flank, ve + args.variant_flank, fasta, tx, args,
                                    args.seed, label_tracks=(k == 0), ref_track=True, highlight=(vs, ve),
                                    track_label="bulk", keep_indels=indel_specs([v]))

            # 5. distributions of junctions outside every event (rare), up to 4 per row; each event's
            # own distributions are drawn under its title, just before its zoom (6)
            values = load_metric_values(manifests, gene, pairs, args.coverage_threshold) if (pairs and manifests) else {}
            if pairs and not manifests:
                print(f"  {gene}: no cohort manifests given; distribution panels skipped.")
            ev_jxns = {j for e in events for j in e["junctions"]}

            def dist_rows(chunk_pairs, title=None):
                """Distribution panels, 4 per row; `title` above the first row. Returns the last figure."""
                fig = None
                for i in range(0, len(chunk_pairs), DIST_PANELS_PER_ROW):
                    chunk = chunk_pairs[i:i + DIST_PANELS_PER_ROW]
                    fig, top, h = pager.next_row()
                    t_off = 0.0
                    if title and i == 0:
                        fig.text(full_left / PAGE_W, (top - 0.17) / PAGE_H, title, fontsize=6.5 * FS,
                                 fontweight="bold", va="bottom", ha="left")
                        t_off = 0.2
                    for (left, w), (metric, jxn) in zip(panel_columns(len(chunk), DIST_PANELS_PER_ROW), chunk):
                        ax = _ax(fig, left + 0.25, top - h + 0.3, w - 0.3, h - 0.55 - t_off)
                        draw_distribution(ax, metric, jxn, values.get((metric, jxn), {}), sample_to_type, type_colors,
                                          outliers_by_pair.get((metric, jxn), set()), sample)
                return fig

            if manifests:
                dist_rows([p_ for p_ in pairs if p_[1] not in ev_jxns])

            # 6. splicing events: one row (bulk); a haplotype event gets its own page with bulk,
            # hap1 and hap2 rows. One title, on the bulk row. The IGV view takes the first two
            # columns; the third holds a sashimi plot of this sample (bulk, and hap1 / hap2 for a
            # haplotype event) and a representative cohort sample (chosen once per gene).
            rep, rep_bam = None, None
            if events:
                psi, depth = load_gene_psi(manifests, gene, args.coverage_threshold)
                rep = choose_representative(gene, str(sample), [str(x) for x in args.samples], sample_to_type,
                                            sample_bams, psi, depth, jxn_outlier_samples.get(gene, set()))
                if rep:
                    rep_bam = pysam.AlignmentFile(sample_bams[rep["sample"]], "rb")
                    print(f"  {gene}: representative {rep['sample']} ({rep['sample_type']}, mean |dPSI| from "
                          f"median = {rep['dist']:.3f}, {rep['n']} candidate(s))")
                else:
                    print(f"  {gene}: no representative sample (no covered junction PSI in other samples)")
            cols3 = panel_columns(3, 3)
            igv_w = cols3[1][0] + cols3[1][1] - full_left
            sash_left, sash_w = cols3[2]
            for e in events:
                s = min(_jxn_coords(j)[1] for j in e["junctions"]) - 1
                t = max(_jxn_coords(j)[2] for j in e["junctions"])
                flank = max(int(round((t - s) * args.junction_flank)), 5)
                title = (f"{' / '.join(e['labels'])}  —  "
                         + "; ".join(f"{j} ({junction_annotation(tx, j)})" for j in e["junctions"]))
                max_chars = int(full_w * 19 / FS)
                if len(title) > max_chars:
                    title = title[:max_chars - 3] + "..."
                # the event's distributions first, under its title
                ev_pairs = [p_ for p_ in pairs if p_[1] in e["junctions"]] if manifests else []
                dist_fig = dist_rows(ev_pairs, title) if ev_pairs else None
                is_hap_event = any(f", hap1]" in c or f", hap2]" in c for c in e["calls"])
                ev_tracks = [("bulk", None)]
                if is_hap_event and gene in hap_bams:
                    ev_tracks += list(zip(("hap1", "hap2"), hap_bams[gene]))
                    pager.fresh_page()
                handles, slots = [], []
                for k, (name, path) in enumerate(ev_tracks):
                    handle = bulk if path is None else pysam.AlignmentFile(path, "rb")
                    handles.append((name, handle))
                    fig, top, h = pager.next_row()
                    slots.append((fig, top, h))
                    # the title sits over the distributions; repeated on the zoom only if that starts a new page
                    show_title = k == 0 and fig is not dist_fig
                    alignment_panel(fig, full_left, igv_w, top, h, title if show_title else "", handle, chrom,
                                    s - flank, t + flank, fasta, tx, args, args.seed,
                                    arcs_junctions=e["junctions"], track_label=name,
                                    show_tx=(k == len(ev_tracks) - 1), title_width=full_w)
                # sashimi in the third column, over the same rows (all on one page)
                fig0, top0, _ = slots[0]
                same = [sl for sl in slots if sl[0] is fig0]
                sash_top = top0 - 0.19 - 0.06
                sash_bottom = same[-1][1] - same[-1][2] + 0.10
                sash_tracks = [dict(label=f"{sample} {name}", bam=handle, color=darker(type_color(sample_to_type.get(str(sample)))))
                               for name, handle in handles[:len(same)]]
                if rep_bam is not None:
                    sash_tracks.append(dict(label=f"{rep['sample']} (representative {rep['sample_type']})",
                                            bam=rep_bam, color=type_color(rep["sample_type"])))
                else:
                    sash_tracks[0]["label"] += "  (no representative sample)"
                draw_sashimi(fig0, sash_left, sash_w, sash_top, sash_bottom, sash_tracks, tx, chrom,
                             s - flank, t + flank, args)
                for name, handle in handles:
                    if handle is not bulk:
                        handle.close()
                if len(ev_tracks) > 1:
                    pager.fresh_page()
            if rep_bam is not None:
                rep_bam.close()

        pager.flush()
        draw_qc_sections(pdf, args, str(sample), shown_name, qc_ctx)
        return pdf.pages


if __name__ == "__main__":
    main()