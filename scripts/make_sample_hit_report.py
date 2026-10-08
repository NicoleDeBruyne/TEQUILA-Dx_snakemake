"""
make_sample_hit_report.py

One PDF per sample (US Letter, landscape): a QC page, then IGV-style views of every candidate gene
hit in merged_all_hits.tsv (default filter set, hit == TRUE, tier <= --max-tier, ordered by `ranking`).

QC page (first page). Every panel compares the sample with the median of the samples of its own
sample type and the median of the whole cohort on the panel (both medians include the sample):
  1. Read counts: stacked on-target / off-target bar; the medians are median mapped reads x median
     on-target rate (_9A).
  2. Sample-type PCA (_9C2): GTEx references faint, cohort samples coloured by sample type, this
     sample as a star.
  3. Read length (_9B): on-target and off-target boxes for the sample, and boxes built from the
     median Q1 / median / Q3 (and min / max) of its sample type and of the cohort. Whiskers run
     to Q1 - 1.5 IQR and Q3 + 1.5 IQR, clipped to the min / max (the per-sample summaries hold no
     individual read lengths).
  4. Full-length ratio per gene (_9E) and 5. junction testability per gene
     (n_canonical_jxns_tested / n_canonical_jxns, from merged_all_hits.tsv): heatmaps with rows
     sample / sample-type median / cohort median, genes sorted by the cohort median; grey = no
     value (e.g. single-exon genes, or too few reads).
A panel whose input is not given is left out.

Every gene starts on a new page and is laid out as rows, 3 rows per page:
  1. Bulk alignments over the whole gene body (the gene's region in the BED).
  2. hap1 and 3. hap2 alignments over the gene body, if the gene was phased and has ASE or a
     hap1 / hap2 outlier junction.
  4. Variant zooms (15 bp either side of each candidate variant), up to 3 panels per row, bulk.
  5. Outlier splicing events, one per row, bulk: the event's junctions plus 10% of their span on
     either side. The IGV view (first two columns) has arcs for the event's junctions only; the
     third column is a sashimi plot of the same window (introns scaled by --intron-scale,
     unspliced reads dropped for multi-exon genes, every junction with >= --sashimi-min-frac x the
     track's maximum coverage) for this sample and a representative cohort sample. An event called
     in hap1 or hap2 gets its own page with bulk, hap1 and hap2 rows (all three in the sashimi
     too). Title: event label(s) and each junction with its annotation (canonical / annotated /
     unannotated).
     Representative sample (one per gene): among the other samples of this sample's type (or, if
     there are none, of the cohort's most common sample type), the one whose AMALGAM transcript
     counts for the gene have the highest Pearson correlation with the candidates' mean profile.
     If the correlation is undefined (e.g. a single transcript), the sample with the median gene
     total is used (r = n/a).
  6. Metric distributions (PSI_approx, PSI, IR, IPA) for each junction x metric that fired,
     across every sample on the panel split by sample type, up to 4 panels per row.

If the gene's canonical transcript has more than one exon, unspliced (monoexonic) reads are left
out of the coverage and alignment tracks. Sashimi arcs run from the last exon base before the
intron (ss1 - 1) to the first exon base after it (ss2 + 1).

Each alignment panel has a coverage track (all reads; positions where a non-reference base has
allele frequency > 0.2 at depth >= 10 are coloured by base, as in IGV), an alignment track
(primary and supplementary alignments, randomly downsampled to --max-reads reads; every mismatch
coloured by base, deletions >= 10 bp as black lines, insertions >= 10 bp as purple ticks, smaller
indels hidden, splice gaps as thin grey lines), a reference-sequence track on
variant zooms, and the gene's canonical transcript (thick = CDS, thin = UTR), drawn once under the
bottom track of a bulk / hap1 / hap2 stack. Samples are always shown by name, never by alias.

Outlier splicing events come from the hit table: junctions are grouped by the event ID that _9L /
_8C assigned within each source (GTEx tissue or cohort) and phasing (older outputs without IDs fall
back to grouping same-label junctions that share a splice site); calls of the same junction set
from different sources, tissues or phasings are shown as one event.
"""
import argparse
import gzip
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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from matplotlib import rcParams

from hits_io import read_hit_rows

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
ARC_COLOR = "#b2182b"
SPAN_COLOR = "#fff2a8"

RASTER_DPI = 300                      # resolution of the few other rasterised artists (coverage fill)
ALIGN_DPI = 1200                      # resolution of the alignment-track images
AF_THRESHOLD = 0.2
AF_MIN_DEPTH = 10
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400   # unmapped, secondary, QC-fail, duplicate (primary + supplementary kept)
MIN_INDEL = 10                             # insertions / deletions shorter than this are not drawn

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
                   help="AMALGAM transcript_amalgam_matrix_raw.tsv per sample type (transcript_id, "
                        "gene_id = gene symbol, one column per sample), for choosing the representative.")
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
    p.add_argument("--pca-coords", default=None,
                   help="QC page: _9C2 {bed_id}_PCA_coords.tsv (sample-type validation PCA coordinates).")
    p.add_argument("--outfile", required=True)
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
        for t_id, tr in cands.items():           # junctions of every transcript (same IDs as the pipeline)
            ex = sorted(tr["exons"])
            for (_, e1), (s2, _) in zip(ex, ex[1:]):
                j = f"{tr['chrom']}_{e1 + 1}_{s2}"
                all_j.add(j)
                if tr["canonical"]:
                    canon_j.add(j)
        out[gene] = dict(chrom=r["chrom"], strand=r["strand"], tid=tid,
                         exons=sorted(r["exons"]), cds=sorted(r["cds"]),
                         canonical_junctions=canon_j, annotated_junctions=all_j)
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
        ax.set_ylabel(label, fontsize=6 * FS, rotation=0, ha="right", va="center", labelpad=4)


def _rgb(c):
    return tuple(int(round(v * 255)) for v in matplotlib.colors.to_rgb(c))


def draw_alignments(ax, reads, n_total, ref_seq, start, end, label):
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
                if ln >= MIN_INDEL and start <= pos <= end:
                    c = min(int((pos - start) * px_per_base), W - 1)
                    img[r0:r1, max(c - 1, 0):c + 1] = ins_c
                qpos += ln
            elif op == 2:
                a, b = max(pos, start), min(pos + ln, end)
                if a < b:
                    c0, c1 = xcols(a, b)
                    if ln >= MIN_INDEL:
                        img[l0:l1, c0:c1] = del_c
                    else:                       # small deletion: drawn as part of the read
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
                elif s in outliers:
                    ax.scatter(x, y, s=9, color=col, edgecolors="black", linewidths=0.4, zorder=3)
                    ax.annotate(s, (x, y), xytext=(3, 0), textcoords="offset points",
                                fontsize=4 * FS, va="center")
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


class Pager:
    """Hands out row slots, 4 per page, with a page header naming the gene."""

    def __init__(self, pdf, header_fn):
        self.pdf, self.header_fn = pdf, header_fn
        self.fig, self.row = None, ROWS_PER_PAGE
        self.part = 0

    def _new_page(self):
        self.flush()
        self.fig = plt.figure(figsize=(PAGE_W, PAGE_H))
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
                    marks=None, spans=None, track_label="", show_tx=True, title_width=None):
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
    # for a spliced gene (multi-exon canonical transcript), unspliced reads are left out of both tracks
    keep = _keep_spliced if (tx is not None and len(tx["exons"]) > 1) else _keep
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
    draw_alignments(ax, reads, n, ref_seq, start, end, "alignments" if label_tracks else "")
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
def load_amalgam_matrix(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.astype(str)
    df["gene_id"] = df["gene_id"].astype(str)
    return df


def choose_representative(gene, sample, samples, sample_to_type, sample_bams, amalgam_paths, cache):
    """The cohort sample whose AMALGAM transcript counts for `gene` correlate best (Pearson) with
    the mean profile of the candidates: the other samples of this sample's type, or, if there are
    none, the samples of the cohort's most common sample type. Returns
    dict(sample, sample_type, r, n) or None."""
    stype = sample_to_type.get(sample)
    usable = [s for s in samples if s != sample and s in sample_bams]
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
    path = amalgam_paths.get(rtype)
    if not path or not os.path.isfile(path):
        print(f"  [WARNING] no AMALGAM transcript matrix for {rtype}; no representative sample.")
        return None
    if path not in cache:
        cache[path] = load_amalgam_matrix(path)
    mat = cache[path]
    cands = [s for s in cands if s in mat.columns]
    rows = mat[mat["gene_id"] == gene]
    if not cands or rows.empty:
        return None
    expr = rows[cands].apply(pd.to_numeric, errors="coerce").fillna(0)
    mean_profile = expr.mean(axis=1)
    best, best_r = None, -np.inf
    for c in cands:
        col = expr[c]
        if len(col) < 2 or col.std() == 0 or mean_profile.std() == 0:
            continue
        r = float(np.corrcoef(col.to_numpy(float), mean_profile.to_numpy(float))[0, 1])
        if r > best_r:
            best, best_r = c, r
    if best is None:
        # correlation undefined (e.g. a single transcript): the sample with the median gene total
        tot = expr.sum(axis=0)
        best = (tot - tot.median()).abs().idxmin()
        best_r = np.nan
    return dict(sample=best, sample_type=rtype, r=best_r, n=len(cands))


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
    keep = _keep_spliced if (tx is not None and len(tx["exons"]) > 1) else _keep
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
QC_ON_COLOR, QC_OFF_COLOR = "#4c72b0", "#c8c8c8"
QC_SAMPLE_COLOR, QC_COHORT_COLOR = "#333333", "#999999"
QC_NAN_COLOR = "#bbbbbb"


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


def qc_read_counts(ax, df, sample, shown, stype, type_samples, all_samples):
    """Stacked on-/off-target bar for the sample, its sample-type median and the cohort median."""
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
    bars = [(shown, df.at[sample, "mapped"], df.at[sample, "on_target_rate"] / 100)]
    for lab, group in ((f"{stype}\nmedian", type_samples), ("cohort\nmedian", all_samples)):
        sub = df.reindex([s for s in group if s in df.index])
        if len(sub):
            bars.append((lab, sub["mapped"].median(), sub["on_target_rate"].median() / 100))
    x = np.arange(len(bars))
    mapped = np.array([b[1] for b in bars], dtype=float)
    rate = np.array([b[2] for b in bars], dtype=float)
    scale, unit = (1e6, "M") if np.nanmax(mapped) >= 1e6 else (1e3, "k")
    on, off = mapped * rate / scale, mapped * (1 - rate) / scale
    ax.bar(x, on, color=QC_ON_COLOR, width=0.65, label="on-target")
    ax.bar(x, off, bottom=on, color=QC_OFF_COLOR, width=0.65, label="off-target")
    for i in range(len(bars)):
        ax.text(x[i], mapped[i] / scale * 1.02, f"{mapped[i] / scale:.1f}{unit}\n{rate[i]:.0%} on",
                ha="center", va="bottom", fontsize=5.5 * FS)
    ax.set_xticks(x)
    ax.set_xticklabels([b[0] for b in bars], fontsize=5.5 * FS)
    ax.set_ylabel(f"mapped reads ({unit})")
    ax.set_ylim(0, np.nanmax(mapped) / scale * 1.4)
    ax.set_title(title, fontsize=7 * FS, fontweight="bold")
    ax.legend(fontsize=5 * FS, frameon=False, loc="upper right", ncol=2, bbox_to_anchor=(1, 1.03))
    _qc_spines(ax)


def qc_read_lengths(ax, df, sample, shown, stype, type_samples, all_samples, type_color):
    """On-/off-target read-length boxes: sample, sample-type median quartiles, cohort median quartiles."""
    title = "Read length\n(whiskers: 1.5 IQR, clipped to min/max)"
    if df is None:
        return _qc_unavailable(ax, title)
    df = df.copy()
    df["sample"] = df["sample"].astype(str)
    for c in ("min", "q1", "median", "q3", "max"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    stats, pos, cols, ticks, groups = [], [], [], [], []
    for g, tt in enumerate(("on_target", "off_target")):
        sub = df[df["targetType"] == tt].set_index("sample")
        rows = [(shown, sub.loc[[sample]] if sample in sub.index else sub.iloc[0:0], QC_SAMPLE_COLOR),
                (stype, sub.reindex([s for s in type_samples if s in sub.index]), type_color),
                ("cohort", sub.reindex([s for s in all_samples if s in sub.index]), QC_COHORT_COLOR)]
        x0 = g * 4
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
        groups.append((x0 + 1, tt.replace("_", "-")))
    if not stats:
        return _qc_unavailable(ax, title, f"{shown} not in the read-length table")
    b = ax.bxp(stats, positions=pos, widths=0.6, patch_artist=True, showfliers=False)
    for patch, c in zip(b["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("black"); patch.set_linewidth(0.6)
    for m in b["medians"]:
        m.set_color("black")
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], fontsize=5 * FS, rotation=30, ha="right")
    ax.set_xlim(-0.6, 6.6)
    for xc, lab in groups:
        ax.text(xc, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=6 * FS, fontweight="bold")
    ax.set_ylabel("read length (bp)")
    ax.set_title(title, fontsize=7 * FS, fontweight="bold", pad=14 * FS)
    _qc_spines(ax)


def qc_pca(ax, df, sample, shown):
    """Sample-type validation PCA (PC1 vs PC2) with this sample as a star."""
    title = "Sample type (junction PSI PCA)"
    if df is None or not {"PC1", "PC2"} <= set(df.columns):
        return _qc_unavailable(ax, title)
    df = df.copy()
    df["sample"] = df["sample"].astype(str)
    handles = []
    for kind, kw in (("reference", dict(s=6, alpha=0.25, linewidths=0, marker="o", zorder=2)),
                     ("query", dict(s=22, alpha=0.9, marker="D", edgecolors="white", linewidths=0.3, zorder=3))):
        sub = df[df["kind"] == kind]
        for grp, gsub in sub.groupby("group", sort=False):
            col = gsub["color"].iloc[0]
            name = textwrap.fill(str(grp), 18) + f" (n={len(gsub)})"
            h = ax.scatter(gsub["PC1"], gsub["PC2"], color=col, label=name, **kw)
            handles.append(h)
    me = df[(df["kind"] == "query") & (df["sample"] == sample)]
    if not me.empty:
        ax.scatter(me["PC1"], me["PC2"], s=140, marker="*", color=me["color"].iloc[0],
                   edgecolors="black", linewidths=0.8, zorder=5)
        ax.annotate(shown, (me["PC1"].iloc[0], me["PC2"].iloc[0]), xytext=(6, 4), textcoords="offset points",
                    fontsize=6 * FS, fontweight="bold", zorder=6)
    else:
        ax.text(0.02, 0.98, f"{shown} not in the PCA", transform=ax.transAxes, va="top", color="#b00",
                fontsize=5.5 * FS)
    v1 = df["PC1_var_pct"].iloc[0] if "PC1_var_pct" in df else np.nan
    v2 = df["PC2_var_pct"].iloc[0] if "PC2_var_pct" in df else np.nan
    ax.set_xlabel(f"PC1 ({v1:.1f}%)"); ax.set_ylabel(f"PC2 ({v2:.1f}%)")
    ax.set_title(title, fontsize=7 * FS, fontweight="bold")
    leg = ax.legend(handles=handles, fontsize=4.5 * FS, frameon=False, loc="center left",
                    bbox_to_anchor=(1.0, 0.5), markerscale=0.8)
    for lh in leg.legend_handles if hasattr(leg, "legend_handles") else leg.legendHandles:
        lh.set_alpha(0.9)
    _qc_spines(ax)


def qc_heatmap(fig, ax, cax, matrix, sample, shown, stype, type_samples, all_samples, title, cmap, cbar_label):
    """Three-row gene heatmap: sample, sample-type median, cohort median; sorted by cohort median."""
    if matrix is None or matrix.empty or sample not in matrix.columns:
        cax.set_visible(False)
        return _qc_unavailable(ax, title, "not available" if matrix is None or matrix.empty
                               else f"{shown} not in the table")
    type_cols = [s for s in type_samples if s in matrix.columns]
    all_cols = [s for s in all_samples if s in matrix.columns] or list(matrix.columns)
    cohort_med = matrix[all_cols].median(axis=1, skipna=True)
    rows = [matrix[sample], matrix[type_cols].median(axis=1, skipna=True), cohort_med]
    order = cohort_med.sort_values(na_position="last", kind="mergesort").index
    m = np.vstack([r.reindex(order).to_numpy(dtype=float) for r in rows])
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(QC_NAN_COLOR)
    im = ax.imshow(np.ma.masked_invalid(m), aspect="auto", cmap=cm, vmin=0, vmax=1, interpolation="nearest")
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels([shown, f"{stype}\nmedian (n={len(type_cols)})", f"cohort\nmedian (n={len(all_cols)})"],
                       fontsize=5.5 * FS, linespacing=0.9)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(order)} panel genes, sorted by cohort median (grey = no value)", fontsize=5.5 * FS)
    ax.set_title(title, fontsize=7 * FS, fontweight="bold", loc="left")
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=5 * FS)
    cb.set_label(cbar_label, fontsize=5 * FS)


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


def draw_qc_page(pdf, args, sample, shown, sample_to_type, type_colors):
    stype = sample_to_type.get(sample, "unknown")
    all_samples = [str(s) for s in args.samples] or [sample]
    type_samples = [s for s in all_samples if sample_to_type.get(s) == stype] or [sample]
    testability = load_testability_matrix(args.hits_tsv)
    if testability is None and not any((args.on_target_tsv, args.read_attributes_tsv,
                                        args.flr_matrix, args.pca_coords)):
        print("  No QC inputs; QC page skipped.")
        return
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H, f"{shown}  |  QC  (sample type: {stype})",
             fontsize=11 * FS, fontweight="bold", va="top")
    fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H,
             f"compared with {stype} (n = {len(type_samples)}) and the whole cohort (n = {len(all_samples)}); "
             "medians include this sample",
             fontsize=6.5 * FS, va="top", color="#333")
    qc_read_counts(_ax(fig, 0.95, 5.1, 2.3, 2.4), _read_tsv(args.on_target_tsv),
                   sample, shown, stype, type_samples, all_samples)
    qc_read_lengths(_ax(fig, 3.95, 5.1, 2.6, 2.2), _read_tsv(args.read_attributes_tsv),
                    sample, shown, stype, type_samples, all_samples, type_colors.get(stype, "#4c72b0"))
    qc_pca(_ax(fig, 7.15, 5.1, 2.2, 2.4), _read_tsv(args.pca_coords), sample, shown)
    qc_heatmap(fig, _ax(fig, 1.75, 2.85, 8.2, 1.05), _ax(fig, 10.05, 2.85, 0.12, 1.05),
               load_flr_matrix(args.flr_matrix), sample, shown, stype, type_samples, all_samples,
               "Full-length ratio per gene", "viridis", "full-length ratio")
    qc_heatmap(fig, _ax(fig, 1.75, 0.85, 8.2, 1.05), _ax(fig, 10.05, 0.85, 0.12, 1.05),
               testability, sample, shown, stype, type_samples, all_samples,
               "Junction testability per gene (fraction of canonical junctions tested)", "magma",
               "fraction tested")
    pdf.savefig(fig, dpi=RASTER_DPI)
    _flush_pdf_images(pdf)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------
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
    amalgam_paths = {}
    for item in args.amalgam_matrices:
        st, _, path = item.partition(":")
        if path:
            amalgam_paths[st] = path
    amalgam_cache = {}
    type_color = lambda st: type_colors.get(st, "#7f7f7f")
    fasta = pysam.FastaFile(args.genome)
    bulk = pysam.AlignmentFile(args.bam, "rb")

    with PdfPages(args.outfile) as pdf:
        draw_qc_page(pdf, args, str(sample), shown_name, sample_to_type, type_colors)
        if not genes:
            fig = plt.figure(figsize=(PAGE_W, PAGE_H))
            fig.text(0.5, 0.5, f"{shown_name}: no candidate gene hits.", ha="center", va="center", fontsize=14 * FS)
            pdf.savefig(fig)
            plt.close(fig)
            return

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
            def header(fig, part, gene=gene, row=row, cats=cats):
                fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H,
                         f"{shown_name}  |  {gene}  |  rank {row.get('ranking', '.')}, tier {row.get('tier', '.')}"
                         + (f"  (page {part})" if part > 1 else ""),
                         fontsize=11 * FS, fontweight="bold", va="top")
                sub = "hits: " + (", ".join(f"{c}={row.get(c)}" for c in cats) or ".")
                if _nonempty(row.get("inheritance_patterns")):
                    sub += f"   |   inheritance: {row.get('inheritance_patterns')}"
                fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H, sub, fontsize=6.5 * FS, va="top", color="#333")

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
            # the canonical transcript is drawn once, under the last (bottom) track
            for k, (name, handle) in enumerate(tracks):
                fig, top, h = pager.next_row()
                alignment_panel(fig, full_left, full_w, top, h,
                                "",
                                handle, chrom, gstart, gend, fasta, tx, args, args.seed,
                                marks=var_marks, spans=ev_spans, track_label=name,
                                show_tx=(k == len(tracks) - 1))
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
                                    track_label="bulk")

            # 5. splicing events: one row (bulk); a haplotype event gets its own page with bulk,
            # hap1 and hap2 rows. One title, on the bulk row. The IGV view takes the first two
            # columns; the third holds a sashimi plot of this sample (bulk, and hap1 / hap2 for a
            # haplotype event) and a representative cohort sample (chosen once per gene).
            rep, rep_bam = None, None
            if events:
                rep = choose_representative(gene, str(sample), [str(x) for x in args.samples], sample_to_type,
                                            sample_bams, amalgam_paths, amalgam_cache)
                if rep:
                    rep_bam = pysam.AlignmentFile(sample_bams[rep["sample"]], "rb")
                    r_txt = "n/a" if np.isnan(rep["r"]) else f"{rep['r']:.3f}"
                    print(f"  {gene}: representative {rep['sample']} ({rep['sample_type']}, r = {r_txt}, "
                          f"{rep['n']} candidate(s))")
                else:
                    print(f"  {gene}: no representative sample (no AMALGAM data or no other samples)")
            cols3 = panel_columns(3, 3)
            igv_w = cols3[1][0] + cols3[1][1] - full_left
            sash_left, sash_w = cols3[2]
            for e in events:
                s = min(_jxn_coords(j)[1] for j in e["junctions"]) - 1
                t = max(_jxn_coords(j)[2] for j in e["junctions"])
                flank = max(int(round((t - s) * args.junction_flank)), 5)
                title = (f"{' / '.join(e['labels'])}  —  "
                         + "; ".join(f"{j} ({junction_annotation(tx, j)})" for j in e["junctions"]))
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
                    alignment_panel(fig, full_left, igv_w, top, h, title if k == 0 else "", handle, chrom,
                                    s - flank, t + flank, fasta, tx, args, args.seed,
                                    arcs_junctions=e["junctions"], track_label=name,
                                    show_tx=(k == len(ev_tracks) - 1), title_width=full_w)
                # sashimi in the third column, over the same rows (all on one page)
                fig0, top0, _ = slots[0]
                same = [sl for sl in slots if sl[0] is fig0]
                sash_top = top0 - 0.19 - 0.06
                sash_bottom = same[-1][1] - same[-1][2] + 0.10
                sash_tracks = [dict(label=f"{sample} {name}", bam=handle, color=type_color(sample_to_type.get(str(sample))))
                               for name, handle in handles[:len(same)]]
                if rep_bam is not None:
                    r_txt = "n/a" if np.isnan(rep["r"]) else f"{rep['r']:.3f}"
                    sash_tracks.append(dict(label=f"{rep['sample']} (representative {rep['sample_type']}, r = {r_txt})",
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

            # 6. metric distributions, up to 4 per row
            if pairs and manifests:
                values = load_metric_values(manifests, gene, pairs, args.coverage_threshold)
                for i in range(0, len(pairs), DIST_PANELS_PER_ROW):
                    chunk = pairs[i:i + DIST_PANELS_PER_ROW]
                    fig, top, h = pager.next_row()
                    for (left, w), (metric, jxn) in zip(panel_columns(len(chunk), DIST_PANELS_PER_ROW), chunk):
                        ax = _ax(fig, left + 0.25, top - h + 0.3, w - 0.3, h - 0.55)
                        draw_distribution(ax, metric, jxn, values.get((metric, jxn), {}), sample_to_type, type_colors,
                                          outliers_by_pair.get((metric, jxn), set()), sample)
            elif pairs:
                print(f"  {gene}: no cohort manifests given; distribution panels skipped.")
        pager.flush()
    bulk.close()
    fasta.close()
    print(f"Saved: {args.outfile}")


if __name__ == "__main__":
    main()
