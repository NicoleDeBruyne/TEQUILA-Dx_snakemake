"""
make_sample_hit_report.py

One PDF per sample (US Letter, landscape) with IGV-style views of every candidate gene
hit in merged_all_hits.tsv (default filter set, hit == TRUE, tier <= --max-tier, ordered by `ranking`).

Every gene starts on a new page and is laid out as rows, 3 rows per page:
  1. Bulk alignments over the whole gene body (the gene's region in the BED).
  2. hap1 and 3. hap2 alignments over the gene body, if the gene was phased and has ASE or a
     hap1 / hap2 outlier junction.
  4. Variant zooms (15 bp either side of each candidate variant), up to 3 panels per row, bulk.
  5. Outlier splicing events, one per row, bulk: the event's junctions plus 10% of their span on
     either side, with sashimi arcs for the event's junctions only. An event called in hap1 or
     hap2 gets its own page with bulk, hap1 and hap2 rows. Title: event label(s) and each
     junction with its annotation (canonical / annotated / unannotated).
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
variant zooms, and the gene's canonical transcript (thick = CDS, thin = UTR).

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
from sample_alias import add_alias_map_arg, parse_alias_map, resolve

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
    p = argparse.ArgumentParser(description="IGV-style PDF report of one sample's top candidate genes.")
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
    p.add_argument("--outfile", required=True)
    add_alias_map_arg(p)
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
        man = pd.read_csv(mpath, sep="\t", dtype=str)
        paths = man.loc[man["gene"] == gene, "result_path"].tolist()
        for rp in paths:
            if not rp or rp == "None" or not os.path.isfile(rp):
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


def draw_transcript(ax, tx, start, end, label=None, marks=None, spans=None):
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.axis("off")
    for a, b in spans or []:
        ax.axvspan(a, b, color=SPAN_COLOR, zorder=0)
    if tx is None:
        ax.text(0.5, 0.5, "no transcript in GTF", transform=ax.transAxes, ha="center", va="center", fontsize=5 * FS)
        return
    exons, cds = tx["exons"], tx["cds"]
    t0, t1 = exons[0][0], exons[-1][1]
    a, b = max(t0, start), min(t1, end)
    if a < b:
        ax.plot([a, b], [0.5, 0.5], color="#1f3a93", linewidth=0.6, zorder=1)
        n_arrows = 12
        marker = ">" if tx["strand"] == "+" else "<"
        for xa in np.linspace(a, b, n_arrows + 2)[1:-1]:
            if not any(s <= xa < e for s, e in exons):
                ax.plot(xa, 0.5, marker=marker, markersize=2.2, color="#1f3a93", zorder=2)
    for s, e in exons:
        if e <= start or s >= end:
            continue
        ax.add_patch(Rectangle((s, 0.32), e - s, 0.36, color="#1f3a93", linewidth=0, zorder=3))
    for s, e in cds:
        if e <= start or s >= end:
            continue
        ax.add_patch(Rectangle((s, 0.15), e - s, 0.70, color="#1f3a93", linewidth=0, zorder=3))
    for x in marks or []:
        ax.plot(x, 0.95, marker="v", markersize=3, color="#d62728", zorder=4, clip_on=False)
    if label:
        ax.text(-0.003, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=5.5 * FS)


def draw_distribution(ax, metric, junction, values, sample_to_type, type_colors, outliers, sample, alias_map):
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
                    ax.annotate(resolve(s, alias_map), (x, y), xytext=(4, 0), textcoords="offset points",
                                fontsize=5 * FS, va="center", fontweight="bold")
                elif s in outliers:
                    ax.scatter(x, y, s=9, color=col, edgecolors="black", linewidths=0.4, zorder=3)
                    ax.annotate(resolve(s, alias_map), (x, y), xytext=(3, 0), textcoords="offset points",
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
                    marks=None, spans=None, track_label=""):
    max_chars = int(width * 19 / FS)
    if len(title) > max_chars:
        title = title[:max_chars - 3] + "..."
    title_h, tick_h, cov_h, ann_h, ref_h, pad = 0.19, 0.13, 0.36, 0.20, (0.15 if ref_track else 0.0), 0.03
    aln_h = h - title_h - tick_h - cov_h - ann_h - ref_h - 4 * pad - 0.06
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
    ax = _ax(fig, left, y - ann_h, width, ann_h)
    draw_transcript(ax, tx, start, end, label=None,
                    marks=marks, spans=spans)


# ----------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------
def main():
    args = parse_args()
    alias_map = parse_alias_map(args.alias_map)
    sample = args.sample
    shown_name = resolve(sample, alias_map)
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
    fasta = pysam.FastaFile(args.genome)
    bulk = pysam.AlignmentFile(args.bam, "rb")

    with PdfPages(args.outfile) as pdf:
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
            for name, handle in tracks:
                fig, top, h = pager.next_row()
                alignment_panel(fig, full_left, full_w, top, h,
                                "",
                                handle, chrom, gstart, gend, fasta, tx, args, args.seed,
                                marks=var_marks, spans=ev_spans, track_label=name)
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
            # hap1 and hap2 rows. One title, on the bulk row.
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
                for k, (name, path) in enumerate(ev_tracks):
                    handle = bulk if path is None else pysam.AlignmentFile(path, "rb")
                    fig, top, h = pager.next_row()
                    alignment_panel(fig, full_left, full_w, top, h, title if k == 0 else "", handle, chrom,
                                    s - flank, t + flank, fasta, tx, args, args.seed,
                                    arcs_junctions=e["junctions"], track_label=name)
                    if path is not None:
                        handle.close()
                if len(ev_tracks) > 1:
                    pager.fresh_page()

            # 6. metric distributions, up to 4 per row
            if pairs and manifests:
                values = load_metric_values(manifests, gene, pairs, args.coverage_threshold)
                for i in range(0, len(pairs), DIST_PANELS_PER_ROW):
                    chunk = pairs[i:i + DIST_PANELS_PER_ROW]
                    fig, top, h = pager.next_row()
                    for (left, w), (metric, jxn) in zip(panel_columns(len(chunk), DIST_PANELS_PER_ROW), chunk):
                        ax = _ax(fig, left + 0.25, top - h + 0.3, w - 0.3, h - 0.55)
                        draw_distribution(ax, metric, jxn, values.get((metric, jxn), {}), sample_to_type, type_colors,
                                          outliers_by_pair.get((metric, jxn), set()), sample, alias_map)
            elif pairs:
                print(f"  {gene}: no cohort manifests given; distribution panels skipped.")
        pager.flush()
    bulk.close()
    fasta.close()
    print(f"Saved: {args.outfile}")


if __name__ == "__main__":
    main()