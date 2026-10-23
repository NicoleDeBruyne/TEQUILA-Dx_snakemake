import argparse
import gzip
import io
import itertools
import os
import random
import re
import textwrap
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.path import Path
from matplotlib.patches import Patch, PathPatch, Rectangle, FancyBboxPatch
from matplotlib import rcParams, font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.colors import LinearSegmentedColormap, to_rgb, to_hex

import tier_rules
from palette import LIGHT_GREY as PALETTE_LIGHT_GREY, SELECT_ORDER as PALETTE_ORDER
from scipy.cluster.hierarchy import leaves_list, linkage


rcParams["pdf.fonttype"] = 42
FONT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts")
if os.path.isdir(FONT_DIR):
    for _f in sorted(os.listdir(FONT_DIR)):
        if _f.lower().endswith((".ttf", ".otf")):
            font_manager.fontManager.addfont(os.path.join(FONT_DIR, _f))
rcParams["font.family"] = "sans-serif"
NARROW_FONT = ["Liberation Sans Narrow", "Liberation Sans", "DejaVu Sans"]
WIDE_FONT = ["Liberation Sans", "DejaVu Sans"]
rcParams["font.sans-serif"] = NARROW_FONT
FS = 1.7
rcParams["font.size"] = 6 * FS

PAGE_W, PAGE_H = 11.0, 8.5
MARGIN_L, MARGIN_R = 0.75, 0.25
MARGIN_T, MARGIN_B = 0.70, 0.25
ROW_FRAC = {"quarter": 0.25, "third": 1 / 3, "half": 0.5, "two_thirds": 2 / 3}
IGV_LABEL_W = 0.75
BOX_BOTTOM, BOX_TOP = 0.42, 0.33
VARIANT_PANELS_PER_ROW = 3
DIST_PANELS_PER_ROW = 4
PANEL_GAP = 0.18

BASE_COLORS = {"A": "#00a600", "C": "#0000ff", "G": "#d17105", "T": "#ff0000", "N": "#888888"}
READ_COLOR = "#bdbdbd"
DEL_COLOR = "#000000"
INS_COLOR = "#8a2be2"
GAP_COLOR = "#9e9e9e"
COV_COLOR = "#a6a6a6"
ARC_COLOR = "#b2182b"
SPAN_COLOR = "#fff2a8"

RASTER_DPI = 300
ALIGN_DPI = 600
AF_THRESHOLD = 0.2
AF_MIN_DEPTH = 10
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400
MIN_INDEL = 10
EVENT_IGV_FRAC = 2 / 3
INDEL_POS_TOL = 5

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
    p.add_argument("--max-reads", type=int, default=1000, help="Reads drawn per alignment track.")
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
                   help="QC page: _10A {bed_id}_on_target_rates.tsv (sample, mapped, on_target_rate, ...).")
    p.add_argument("--read-attributes-tsv", default=None,
                   help="QC page: _10B {bed_id}_read_attributes.tsv (sample, targetType, min, q1, median, q3, max).")
    p.add_argument("--flr-matrix", default=None,
                   help="QC page: _10E {bed_id}_full_length_ratio_matrix.tsv (gene x sample avgFLR).")
    p.add_argument("--sample-coverage-tsv", default=None,
                   help="The sample's own qc/{sample}_avg_coverage.tsv (gene, length, ...); its canonical-transcript "
                        "lengths label the per-gene coverage box plots.")
    p.add_argument("--coverage-matrix", default=None,
                   help="QC: _10D {bed_id}_breadth_coverage_matrix.tsv (gene x sample fraction of canonical-"
                        "transcript bases at depth >= merge_min_dp_snv).")
    p.add_argument("--coverage-depth", type=int, default=50,
                   help="Depth threshold behind --coverage-matrix (merge_min_dp_snv), for labels only.")
    p.add_argument("--expression-matrices", nargs="*", default=[], metavar="SAMPLE_TYPE:PATH",
                   help="Each sample type's by-assignment gene_assignment_matrix_motr.tsv (gene x sample), "
                        "for the expression box plot beside each gene's bulk view.")
    p.add_argument("--pca-coords", default=None,
                   help="QC page: _10C2 {bed_id}_PCA_coords.tsv (sample-type validation PCA coordinates).")
    p.add_argument("--hpo-terms-tsv", default=None,
                   help="The sample's resolved HPO terms from rule _6E ({sample}_hpo_terms.tsv), listed on the "
                        "overview page.")
    p.add_argument("--wgs-concordance-tsv", default=None,
                   help="QC page: _10P {bed_id}_wgs_concordance_summary.tsv (RNA vs WGS sample identity check).")
    p.add_argument("--outfile", required=True)
    p.add_argument("--html-out", default=None,
                   help="Also write an offline HTML review page (the same pages as images, with a review "
                        "status and notes per gene) to this path.")
    p.add_argument("--html-dpi", type=int, default=200, help="Resolution of the page images in --html-out.")
    p.add_argument("--previous-reviews", default=None,
                   help="Reviews TSV exported from an earlier review page for this sample; its statuses and "
                        "notes are carried into --html-out and marked as coming from an earlier run.")
    p.add_argument("--tx-usage-matrix", default=None,
                   help="Isoform pages: _9F <sample_type>_transcript_usage_matrix.tsv (isoform usage per sample).")
    p.add_argument("--tx-gene-status", default=None,
                   help="Isoform pages: _9F <sample_type>_transcript_usage_gene_status.tsv (why a gene was not "
                        "tested, and each sample's gene read count).")
    p.add_argument("--tx-isoform-outliers", default=None, help="Isoform pages: _9G individual isoform outliers.")
    p.add_argument("--tx-landscape-outliers", default=None, help="Isoform pages: _9G isoform mix outliers.")
    p.add_argument("--tx-gtf", default=None,
                   help="Isoform pages: the AMALGAM transcript GTF the isoforms were quantified against "
                        "(_9B annotation/filtered.gtf.gz, tabix-indexed); draws each isoform's structure.")
    p.add_argument("--tx-method", default=None, help="The transcript analysis method (beta_binomial or modified_zscore).")
    p.add_argument("--tx-n-threshold", type=int, default=10,
                   help="Samples with enough gene reads a gene needed to be tested (for the 'not available' text).")
    p.add_argument("--tx-min-gene-reads", type=float, default=50)
    p.add_argument("--tx-min-usage", type=float, default=0.05)
    p.add_argument("--tx-min-reads", type=float, default=10)
    p.add_argument("--tx-unavailable", default=None,
                   help="Why the transcript analysis is not available for this sample's cohort; shown on the "
                        "isoform pages instead of the plots.")
    return p.parse_args()


def _nonempty(v):
    return v is not None and not (isinstance(v, float) and np.isnan(v)) and str(v).strip().lower() not in _MISSING


def _items(cell):
    if not _nonempty(cell):
        return []
    return [x.strip() for x in str(cell).split(";")]


_TISSUE_RE = re.compile(r"\s*([^,()]+?)\s*\(([^)]+)\)")


def _tissue_values(cell):
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
    genes = set(genes)
    txs = defaultdict(lambda: defaultdict(lambda: {"exons": [], "cds": [], "canonical": False}))
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            i = line.find('gene_name "')
            if i < 0:
                continue
            i += 11
            if line[i:line.find('"', i)] not in genes:
                continue
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
        merged = []
        for a, b in sorted(e for tr in cands.values() for e in tr["exons"]):
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        for t_id, tr in cands.items():
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
    m = _EID_RE.search(str(text))
    return (str(text)[:m.start()], m.group(1)) if m else (str(text), None)


def junction_calls(row):
    calls = []
    for ph in ("bulk", "hap1", "hap2"):
        jxns = _items(row.get(f"{ph}_jxns"))
        events = _items(row.get(f"{ph}_jxn_event"))
        deltas = _items(row.get(f"{ph}_deltaPSI"))
        nsamp = _items(row.get(f"{ph}_jxn_nsamples"))
        for i, j in enumerate(jxns):
            if not _nonempty(j):
                continue
            dmap = dict((t, v) for v, t in _tissue_values(deltas[i] if i < len(deltas) else ""))
            nmap = dict((t, v) for v, t in _tissue_values(nsamp[i] if i < len(nsamp) else ""))
            for lab, tissue in _tissue_values(events[i] if i < len(events) else ""):
                lab, eid = _split_event_id(lab)
                for l in lab.split(","):
                    l = l.strip()
                    if l in GTEX_LABELS:
                        calls.append(dict(source="GTEx", tissue=tissue, phasing=ph, junction=j, label=l,
                                          metric="junction_PSI_approx", delta=dmap.get(tissue, ""), eid=eid,
                                          nsamples=nmap.get(tissue, "")))
        cj = _items(row.get(f"cohort_{ph}_jxns"))
        cev = _items(row.get(f"cohort_{ph}_jxn_event"))
        cdel = {m: _items(row.get(f"cohort_{ph}_{s}")) for m, s in COHORT_DELTA_SUFFIX.items()}
        cns = _items(row.get(f"cohort_{ph}_jxn_nsamples"))
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
                                  metric=metric, delta=dl[i] if i < len(dl) else "", eid=eid,
                                  nsamples=cns[i] if i < len(cns) else ""))
    return calls


def build_events(calls):
    events = {}

    def _add_event(g, labels, src, phasing, cs):
        key = frozenset(g)
        ev = events.setdefault(key, dict(junctions=sorted(g, key=lambda j: _jxn_coords(j)[1:]),
                                         labels=[], calls=[], deltas=defaultdict(list), nsamples=0))
        for label in labels:
            if label not in ev["labels"]:
                ev["labels"].append(label)
        ev["calls"].append(f"{'/'.join(dict.fromkeys(labels))} [{src}, {phasing}]")
        for c in cs:
            if c["junction"] in key and _nonempty(c["delta"]):
                ev["deltas"][c["junction"]].append(c["delta"])
            if c["junction"] in key:
                n = pd.to_numeric(c.get("nsamples", ""), errors="coerce")
                if pd.notna(n):
                    ev["nsamples"] = max(ev["nsamples"], int(n))

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


_MANIFEST_CACHE = {}
_RESULT_CACHE = {}


def _manifest_paths(mpath, gene):
    man = _MANIFEST_CACHE.get(mpath)
    if man is None:
        df = pd.read_csv(mpath, sep="\t", dtype=str, keep_default_na=False)
        man = defaultdict(list)
        for g, rp in zip(df["gene"], df["result_path"]):
            man[g].append(rp)
        _MANIFEST_CACHE[mpath] = man
    return man.get(gene, [])


def _result_table(rp):
    if rp not in _RESULT_CACHE:
        _RESULT_CACHE[rp] = pd.read_csv(rp, sep="\t", dtype=str)
    return _RESULT_CACHE[rp]


def load_metric_values(manifests, gene, pairs, coverage_threshold):
    vals = defaultdict(dict)
    needed = {j for _, j in pairs}
    for stype, mpath in manifests.items():
        if not os.path.isfile(mpath):
            continue
        for rp in _manifest_paths(mpath, gene):
            rp = str(rp).strip()
            if not rp or rp in ("None", "nan", "NA") or not os.path.isfile(rp):
                continue
            df = _result_table(rp)
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


def _keep(read):
    return not (read.flag & SKIP_FLAGS)


def _is_spliced(read):
    return any(op == 3 for op, _ in (read.cigartuples or ()))


def _keep_spliced(read):
    return _keep(read) and _is_spliced(read)


def sample_reads(bam, chrom, start, end, max_reads, seed, junctions=None):
    rng = random.Random(seed)
    res, n = [], 0
    want = {}
    for j in junctions or ():
        _, s1, s2 = _jxn_coords(j)
        want[(s1, s2)] = j
    counts = {j: 0 for j in junctions or ()}
    for r in bam.fetch(chrom, max(0, start), end):
        if r.flag & SKIP_FLAGS:
            continue
        n += 1
        if len(res) < max_reads:
            res.append(r)
        else:
            k = rng.randrange(n)
            if k < max_reads:
                res[k] = r
        if want:
            cig = r.cigartuples
            if cig is None:
                continue
            pos = r.reference_start
            for op, ln in cig:
                if op in (0, 2, 7, 8):
                    pos += ln
                elif op == 3:
                    j = want.get((pos + 1, pos + ln))
                    if j is not None:
                        counts[j] += 1
                    pos += ln
    res.sort(key=lambda r: (r.reference_start, r.reference_end or 0))
    return res, n, counts


def coverage_counts(bam, chrom, start, end, keep=_keep):
    acgt = bam.count_coverage(chrom, max(0, start), end, quality_threshold=0, read_callback=keep)
    return np.array(acgt, dtype=np.int64)


def _fmt_pos(x, _=None):
    return f"{int(x):,}"


_BASE_INDEX = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate("ACGT"):
    _BASE_INDEX[ord(_b)] = _i


def draw_coverage(ax, acgt, ref_seq, start, end, label, arcs=None, show_xticks=True):
    depth = acgt.sum(axis=0)
    L = end - start
    x = np.arange(start, end)
    if L > 4000:
        nb = 2000
        edges = np.linspace(0, L, nb + 1).astype(int)
        mx = np.maximum.reduceat(depth, edges[:-1]) if L else depth
        xs = np.repeat(start + edges, 2)[1:-1]
        ys = np.repeat(mx, 2)
        ax.fill_between(xs, ys, color=COV_COLOR, linewidth=0, rasterized=True)
    else:
        x0 = x.astype(float)
        x1 = x0 + 1.0
        d = depth.astype(float)
        zero = np.zeros_like(d)
        verts = np.stack([np.column_stack(p_) for p_ in ((x0, zero), (x1, zero), (x1, d), (x0, d))], axis=1)
        ax.add_collection(PolyCollection(verts, closed=True, facecolors=COV_COLOR, edgecolors=COV_COLOR,
                                         linewidths=0), autolim=False)
    top = max(int(depth.max()) if depth.size else 0, 1)

    if ref_seq:
        ref_idx = _BASE_INDEX[np.frombuffer(ref_seq.upper().encode("latin-1", "replace"), dtype=np.uint8)]
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
            xa, xb = a + 0.5, b + 0.5
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


def fit_size(text, width_in, size, weight="normal", min_size=3.0):
    prop = FontProperties(family=WIDE_FONT, size=size, weight=weight)
    longest = max(TextPath((0, 0), line, prop=prop).get_extents().width / 72.0 for line in str(text).split("\n"))
    return size if longest <= width_in else max(min_size, size * width_in / longest)


def _track_label(ax, label):
    if "\n" not in label:
        ax.set_ylabel(label, fontsize=5.5 * FS, rotation=0, ha="right", va="center", labelpad=4,
                      fontweight="bold", fontfamily=WIDE_FONT)
        return
    name, rest = label.split("\n", 1)
    tr = matplotlib.transforms.offset_copy(ax.transAxes, fig=ax.figure, x=-4, y=0, units="points")
    size = fit_size("alignments", ax.get_position().x0 * PAGE_W - 4 / 72 - 0.06, 5.5 * FS, "bold")
    ax.text(0, 0.5, name, transform=tr, fontsize=size, fontweight="bold", ha="right", va="bottom",
            fontfamily=WIDE_FONT)
    ax.text(0, 0.5, rest, transform=tr, fontsize=size, fontweight="bold", ha="right", va="top",
            fontfamily=WIDE_FONT)


def indel_specs(variants):
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
    return tuple(int(round(v * 255)) for v in matplotlib.colors.to_rgb(c)) + (255,)


def draw_alignments(ax, reads, n_total, ref_seq, start, end, label, keep_indels=()):
    L = max(end - start, 1)
    gap = max(1, int(L * 0.003))
    row_ends, placed = [], []
    for r in reads:
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
    img = np.full((H, W, 4), 255, dtype=np.uint8)
    row_px = H / nrows
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
    mm_w = max(1, int(round(px_per_base)), int(W / 1500))
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
                    else:
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
        for x, base in mism:
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
    down = f"downsampled\n{n_total:,} \u2192 {len(reads):,}" if n_total > len(reads) else ""
    if label and down:
        tr = matplotlib.transforms.offset_copy(ax.transAxes, fig=ax.figure, x=-4, y=0, units="points")
        room = ax.get_position().x0 * PAGE_W - 4 / 72 - 0.06
        ax.text(0, 0.5, label, transform=tr, fontsize=fit_size(label, room, 5.5 * FS, "bold"), ha="right",
                va="bottom", fontweight="bold", fontfamily=WIDE_FONT)
        ax.text(0, 0.5, down, transform=tr, fontsize=fit_size(down, room, 4.5 * FS), color="#555", ha="right",
                va="top", linespacing=1.1, fontfamily=WIDE_FONT)
    elif label:
        ax.set_ylabel(label, fontsize=5.5 * FS, rotation=0, ha="right", va="center", labelpad=4, fontweight="bold",
                      fontfamily=WIDE_FONT)
    elif down:
        ax.text(0.995, 0.98, down.replace("\n", " "), transform=ax.transAxes, ha="right", va="top",
                fontsize=4.5 * FS, color="#333", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#bbb", lw=0.4,
                                                           alpha=0.9))


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
        n_arrows = max(3, int(round(ax.get_position().width * PAGE_W * 1.3)))
        marker = ">" if tx["strand"] == "+" else "<"
        if xform is not None:
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


METRIC_COUNT_COLS = {
    "junction_PSI":           ("junction_usage", "junction_coverage"),
    "junction_PSI_approx":    ("junction_usage", "junction_coverage_approx"),
    "5ss_IR_ratio":           ("5ss_IR_count", "5ss_coverage"),
    "3ss_IR_ratio":           ("3ss_IR_count", "3ss_coverage"),
    "junction_full_IR_ratio": ("junction_full_IR_count", "junction_coverage"),
    "junction_IPA_ratio":     ("junction_IPA_count", "5ss_coverage"),
}


def load_sample_counts(manifests, gene, pairs, sample):
    out = {}
    needed = {j for _, j in pairs}
    for mpath in manifests.values():
        if not os.path.isfile(mpath):
            continue
        for rp in _manifest_paths(mpath, gene):
            rp = str(rp).strip()
            if not rp or rp in ("None", "nan", "NA") or not os.path.isfile(rp):
                continue
            df = _result_table(rp)
            if not {"sample", "phasing", "junction"} <= set(df.columns):
                continue
            df = df[(df["sample"].astype(str) == str(sample)) & (df["phasing"] == "bulk") & df["junction"].isin(needed)]
            for metric, jxn in pairs:
                num, den = METRIC_COUNT_COLS.get(metric, (None, None))
                cols = [num, num + "_unique", den, den + "_unique"] if num else []
                if (metric, jxn) in out or not cols or not set(cols) <= set(df.columns):
                    continue
                sub = df[df["junction"] == jxn]
                if sub.empty:
                    continue
                vals = pd.to_numeric(sub.iloc[0][cols], errors="coerce")
                if vals.notna().all():
                    out[(metric, jxn)] = tuple(int(v) for v in vals)
    return out


def load_sample_hap_values(manifests, gene, pairs, sample, coverage_threshold):
    out = defaultdict(dict)
    needed = {j for _, j in pairs}
    for mpath in manifests.values():
        if not os.path.isfile(mpath):
            continue
        for rp in _manifest_paths(mpath, gene):
            rp = str(rp).strip()
            if not rp or rp in ("None", "nan", "NA") or not os.path.isfile(rp):
                continue
            df = _result_table(rp)
            if not {"sample", "phasing", "junction"} <= set(df.columns):
                continue
            df = df[(df["sample"].astype(str) == str(sample)) & df["phasing"].isin(["hap1", "hap2"])
                    & df["junction"].isin(needed)]
            for metric, jxn in pairs:
                rc, dc = _METRIC_RESCALED_COL.get(metric), _METRIC_DENOM_COL.get(metric)
                if rc not in df.columns or dc not in df.columns:
                    continue
                for _, r in df[df["junction"] == jxn].iterrows():
                    v, d = pd.to_numeric(r[rc], errors="coerce"), pd.to_numeric(r[dc], errors="coerce")
                    if pd.notna(v) and pd.notna(d) and d >= coverage_threshold:
                        out[(metric, jxn)].setdefault(r["phasing"], float(v))
    return out


def _counts_line(counts):
    num, num_u, den, den_u = counts
    return f"{num:,} ({num_u:,} unique) out of {den:,} ({den_u:,} unique)"


def draw_distribution(ax, metric, junction, values, sample, shown, tgroups, counts=None, hap=None):
    title = f"{METRIC_SHORT.get(metric, metric)} | {junction}"
    if not values:
        _qc_unavailable(ax, title, "no samples with\nenough coverage")
    else:
        sample_box(ax, pd.Series(values, dtype=float), sample, shown, tgroups, title)
        for dx, name in ((-0.27, "hap1"), (0.27, "hap2")):
            v = (hap or {}).get(name)
            if v is None or np.isnan(v):
                continue
            ax.scatter([dx], [v], s=STAR_SIZE * 0.4, marker="*", color=QC_SAMPLE_COLOR, edgecolors="black",
                       linewidths=0.4, zorder=6)
            ax.annotate(name[-1], (dx, v), xytext=(4 if dx > 0 else -4, 0), textcoords="offset points",
                        ha="left" if dx > 0 else "right", va="center", fontsize=4 * FS, color="#222",
                        fontfamily=WIDE_FONT, zorder=6)
    if counts is None:
        ax.set_title(title, fontsize=5.5 * FS, fontweight="normal", pad=2, fontfamily=NARROW_FONT)
        return
    ax.set_title(title, fontsize=5.5 * FS, fontweight="normal", pad=4 + 6.5 * FS, fontfamily=NARROW_FONT)
    line = ax.text(0.5, 1.0, _counts_line(counts), transform=ax.transAxes, ha="center", va="bottom",
                   fontsize=4.8 * FS, color="#333")
    line.set_transform(matplotlib.transforms.offset_copy(ax.transAxes, fig=ax.figure, y=2, units="points"))
    renderer = ax.figure.canvas.get_renderer()
    allowed = (ax.get_position().width * PAGE_W + PANEL_GAP + 0.3) * ax.figure.dpi
    size = 4.8 * FS
    while line.get_window_extent(renderer=renderer).width > allowed and size > 3.4 * FS:
        size -= 0.2 * FS
        line.set_fontsize(size)


def _flush_pdf_images(pdf):
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


REVIEW_CTX = ("overview", None)


def link_text(fig, x_in, y_in, text, target, ha="right"):
    t = fig.text(x_in / PAGE_W, y_in / PAGE_H, text, fontsize=6 * FS, color="#0068a9", ha=ha, va="top")
    t._pdf_link = target
    return t


def _link_rects(fig):
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
    y, x = PAGE_H - 0.16, PAGE_W - MARGIN_R
    items = ([(f"next gene: {next_gene} \u00bb", ("gene", next_gene))] if next_gene else []) \
        + [("contents", ("contents",))] \
        + ([(f"\u00ab previous gene: {prev_gene}", ("gene", prev_gene))] if prev_gene else [])
    r = fig.canvas.get_renderer()
    for text, target in items:
        t = link_text(fig, x, y, text, target)
        x -= t.get_window_extent(renderer=r).width / fig.dpi + 0.25


def draw_contents(pdf, shown, genes):
    global REVIEW_CTX
    fp_ = pdf.first_page
    qc_sections = sorted(((n, i) for (k, n), i in fp_.items() if k == "qc"), key=lambda x: x[1])
    rows = [("Overview", None, None, None, ("overview", None))]
    rows += [(g["gene"], g["rank"], g["tier"], g["hits"], ("gene", g["gene"])) for g in genes
             if ("gene", g["gene"]) in fp_]
    for n, _ in qc_sections:
        rows.append((n, None, None, None, ("qc", n)))
        if n == OTHER_HITS:
            rows += [(g, None, None, "isoform outlier", ("ohit", g)) for g in pdf.ohit_pages]
    per_page = 30
    n_cont = max(1, -(-len(rows) // per_page))
    REVIEW_CTX = ("contents", None)
    for k in range(n_cont):
        fig = plt.figure(figsize=(PAGE_W, PAGE_H))
        fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H, f"{shown}  |  Contents" + (f"  (page {k + 1})" if k else ""),
                 fontsize=11 * FS, fontweight="bold", va="top", fontfamily=WIDE_FONT)
        fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H, "click a name to jump to its pages; every page links "
                 "back here (top right)", fontsize=6.5 * FS, va="top", color="#333", fontfamily=WIDE_FONT)
        y = PAGE_H - 0.95
        for x, h in ((MARGIN_L, "Rank"), (MARGIN_L + 0.5, "Gene / section"), (MARGIN_L + 2.3, "Tier"),
                     (MARGIN_L + 2.8, "Hits"), (PAGE_W - MARGIN_R - 0.4, "Page")):
            fig.text(x / PAGE_W, y / PAGE_H, h, fontsize=6.5 * FS, fontweight="bold", va="top")
        y -= 0.24
        prev_kind = "overview"
        for name, rank, tier, hits, target in rows[k * per_page:(k + 1) * per_page]:
            kind = "qc" if target[0] == "ohit" else target[0]
            if kind != prev_kind:
                y -= 0.06
                prev_kind = kind
            i = pdf.ohit_pages.get(target[1]) if target[0] == "ohit" else fp_.get(target)
            if rank:
                fig.text(MARGIN_L / PAGE_W, y / PAGE_H, rank, fontsize=6.5 * FS, va="top")
            link_text(fig, MARGIN_L + 0.5, y, name, target, ha="left").set_fontsize(6.5 * FS)
            if tier:
                fig.text((MARGIN_L + 2.3) / PAGE_W, y / PAGE_H, tier, fontsize=6.5 * FS, va="top")
            if hits:
                h = hits if len(hits) <= 95 else hits[:92] + "..."
                fig.text((MARGIN_L + 2.8) / PAGE_W, y / PAGE_H, h, fontsize=5.5 * FS, va="top", color="#333")
            if i is not None:
                fig.text((PAGE_W - MARGIN_R - 0.4) / PAGE_W, y / PAGE_H, str(i + 1 + (n_cont if i > 0 else 0)),
                         fontsize=6.5 * FS, va="top")
            y -= 0.205
        pdf.savefig(fig)
        plt.close(fig)
    return rows


_OBJ_AT = re.compile(rb"(\d+) 0 obj\s*(.*?)\s*endobj", re.S)


def _xref_offsets(data):
    offsets = {}
    try:
        start = int(re.search(rb"startxref\s+(\d+)", data[data.rindex(b"startxref"):]).group(1))
        if data[start:start + 4] != b"xref":
            return offsets
        end = data.index(b"trailer", start)
        tokens = data[start + 4:end].split()
        i = 0
        while i + 1 < len(tokens):
            first, count = int(tokens[i]), int(tokens[i + 1])
            i += 2
            for k in range(count):
                off, _gen, kind = tokens[i:i + 3]
                i += 3
                if kind == b"n":
                    offsets[first + k] = int(off)
    except (ValueError, AttributeError):
        return {}
    return offsets


def add_pdf_navigation(path, cap, rows):
    with open(path, "rb") as fh:
        data = fh.read()
    offsets = _xref_offsets(data)

    def obj(num):
        off = offsets.get(num)
        if off is not None:
            m = _OBJ_AT.match(data, off)
            if m and int(m.group(1)) == num:
                return m.group(2)
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
    except Exception as e:
        print(f"  [WARNING] PDF navigation not added ({e!r}).")
        return
    if len(kids) != cap.n:
        print(f"  [WARNING] PDF navigation not added (expected {cap.n} pages, found {len(kids)}).")
        return

    cont = [i for i in range(cap.n) if cap.page_ctx[i][0] == "contents"]
    order = [0] + cont + [i for i in range(1, cap.n) if i not in cont]
    pos = {orig: new for new, orig in enumerate(order)}

    letter_pages = getattr(cap, "letter_pages", {})
    ohit_pages = getattr(cap, "ohit_pages", {})

    def dest(target):
        if target[0] == "letter":
            i = letter_pages.get(target[1])
        elif target[0] == "ohit":
            i = ohit_pages.get(target[1])
        else:
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

    def esc(t):
        return t.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)").encode("latin-1", "replace")
    tops = [("Overview", ("overview", None), [])]
    if cont:
        tops.append(("Contents", ("contents",), []))
    genes_ = [(f"{r[1]}. {r[0]} (tier {r[2]})", r[4]) for r in rows if r[4][0] == "gene" and dest(r[4]) is not None]
    qcs = [(r[0], r[4]) for r in rows
           if r[4][0] == "qc" and r[4][1] not in (OTHER_GENES, OTHER_HITS) and dest(r[4]) is not None]
    ohits = [(g, ("ohit", g)) for g in ohit_pages if dest(("ohit", g)) is not None]
    letters = [(f"{c}", ("letter", c)) for c in sorted(letter_pages, key=lambda c: (c != "#", c))]
    if genes_:
        tops.append(("Candidate genes", genes_[0][1], genes_))
    if ohits and dest(("qc", OTHER_HITS)) is not None:
        tops.append((OTHER_HITS, ("qc", OTHER_HITS), ohits))
    if qcs:
        tops.append(("QC and expression", qcs[0][1], qcs))
    if letters and dest(("qc", OTHER_GENES)) is not None:
        tops.append((OTHER_GENES, ("qc", OTHER_GENES), letters))
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

    base = len(data) if data.endswith(b"\n") else len(data) + 1
    out = bytearray()
    new_offsets = {}
    for n in sorted(new_objs):
        new_offsets[n] = base + len(out)
        out += b"%d 0 obj\n" % n + new_objs[n] + b"\nendobj\n"
    xref = base + len(out)
    out += b"xref\n0 1\n0000000000 65535 f \n"
    for n in sorted(new_offsets):
        out += b"%d 1\n%010d 00000 n \n" % (n, new_offsets[n])
    out += b"trailer\n<< /Size %d /Root %d 0 R" % (nxt[0], root)
    if info:
        out += b" /Info %s 0 R" % info.group(1)
    out += b" /Prev %d >>\nstartxref\n%d\n%%%%EOF\n" % (prev, xref)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(data)
        if not data.endswith(b"\n"):
            fh.write(b"\n")
        fh.write(out)
    os.replace(tmp, path)
    print(f"  contents page, bookmarks and {sum(len(v) for v in per_page.values())} page link(s) added")


def _shrink_images(fig, dpi):
    saved = []
    for ax in fig.axes:
        for im in ax.images:
            a = im.get_array()
            if not isinstance(a, np.ndarray) or a.dtype != np.uint8 or a.ndim != 3:
                continue
            bb = ax.get_window_extent()
            target_w = bb.width / fig.dpi * dpi
            f = int(a.shape[1] // max(target_w, 1.0))
            if f < 2:
                continue
            h, w = a.shape[0] // f * f, a.shape[1] // f * f
            small = a[:h, :w].reshape(h // f, f, w // f, f, a.shape[2]).mean(axis=(1, 3))
            saved.append((im, a))
            im.set_data(np.round(small).astype(np.uint8))
    return saved


def _compact_png(png):
    try:
        from PIL import Image
    except ImportError:
        return png
    im = Image.open(io.BytesIO(png)).convert("RGB")
    q = im.quantize(colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE)
    pal = q.getpalette()[:768]
    q.putpalette([255 if min(pal[i:i + 3]) >= 246 else v for i in range(0, len(pal), 3) for v in pal[i:i + 3]])
    out = io.BytesIO()
    q.save(out, format="PNG", compress_level=6)
    return out.getvalue() if out.tell() < len(png) else png


class CapturePdf:

    def __init__(self, pdf, dpi):
        self._pdf, self._dpi, self.pages = pdf, dpi, []
        self.n = 0
        self.first_page = {}
        self.links = []
        self.page_ctx = []
        self.letter_pages = {}
        self.ohit_pages = {}

    def savefig(self, fig, *a, **k):
        ctx = getattr(fig, "_review_ctx", REVIEW_CTX)
        self.first_page.setdefault(ctx, self.n)
        self.page_ctx.append(ctx)
        self.links += [(self.n, r, t) for r, t in _link_rects(fig)]
        self._pdf.savefig(fig, *a, **k)
        self.n += 1
        if self._dpi and ctx[0] != "contents":
            nav = [t for t in fig.texts if getattr(t, "_pdf_link", None) is not None]
            for t in nav:
                t.set_visible(False)
            ims = [im for ax in fig.axes for im in ax.images]
            old = [im.get_interpolation() for im in ims]
            for im in ims:
                im.set_interpolation("antialiased")
            shrunk = _shrink_images(fig, self._dpi)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=self._dpi)
            for im, a in shrunk:
                im.set_data(a)
            for im, o in zip(ims, old):
                im.set_interpolation(o)
            for t in nav:
                t.set_visible(True)
            self.pages.append((ctx, _compact_png(buf.getvalue())))

    def __getattr__(self, name):
        return getattr(self._pdf, name)


def pdf_gene_header(shown, gene, row, n_variants, n_events):
    title = f"{shown} | {gene} | Rank {row.get('ranking', '.')}, Tier {row.get('tier', '.')}"
    parts = []
    if _nonempty(row.get("inheritance_patterns")):
        inh = str(row.get("inheritance_patterns"))
        if ("AD" in inh or "XLD" in inh) and str(row.get("haploinsufficient", "")).strip().upper() == "TRUE":
            inh += " (haploinsufficient)"
        parts.append(f"inheritance: {inh}")
    if str(row.get("ASE", "")).strip().upper() == "TRUE":
        parts.append("Allelic imbalance")
    if str(row.get("isoform_outlier", "")).strip().upper() == "TRUE":
        parts.append("Isoform outlier")
    if n_variants:
        parts.append(f"{n_variants} candidate variant(s)")
    if n_events:
        parts.append(f"{n_events} outlier splicing event(s)")
    return title, " | ".join(parts)


def amelie_summary(row):
    rank, n = str(row.get("AMELIE_rank", ".")), str(row.get("AMELIE_n_genes", "."))
    if rank == "error":
        return "AMELIE unavailable", []
    if not _nonempty(n):
        return "", []
    if not _nonempty(rank):
        return "AMELIE: no literature match", []
    text = f"AMELIE rank {rank}/{n}"
    if _nonempty(row.get("AMELIE_score")):
        text += f" (score {row.get('AMELIE_score')})"
    pmids = [p for p in str(row.get("AMELIE_pmids", ".")).split(",") if _nonempty(p)]
    return text, pmids


def load_hpo_labels(path):
    df = _read_tsv(path, dtype=str, keep_default_na=False)
    if df is None or "hpo_id" not in df.columns:
        return []
    labels = []
    for row in df.itertuples(index=False):
        if not _nonempty(row.hpo_id):
            continue
        label = row.hpo_name if _nonempty(row.hpo_name) else row.hpo_id
        if label not in labels:
            labels.append(label)
    return labels


def right_header_text(fig, left_artist, prefix, items, sep=", ", y_in=0.30, **kw):
    renderer = fig.canvas.get_renderer()
    left_limit = left_artist.get_window_extent(renderer=renderer).x1 / fig.dpi + 0.3
    for n in range(len(items), -1, -1):
        shown = sep.join(items[:n])
        extra = f" (+{len(items) - n} more)" if n < len(items) else ""
        text = prefix + shown + extra if (shown or not items) else prefix.rstrip(" :·") + extra
        t = fig.text((PAGE_W - MARGIN_R) / PAGE_W, 1 - y_in / PAGE_H, text, ha="right", va="top", **kw)
        if t.get_window_extent(renderer=renderer).x0 / fig.dpi >= left_limit or n == 0:
            return t
        t.remove()


class Pager:

    def __init__(self, pdf, header_fn):
        self.pdf, self.header_fn = pdf, header_fn
        self.fig, self.used = None, 1.0
        self.part = 0

    def _new_page(self):
        self.flush()
        self.fig = plt.figure(figsize=(PAGE_W, PAGE_H))
        self.fig._review_ctx = REVIEW_CTX
        self.part += 1
        self.header_fn(self.fig, self.part)
        self.used = 0.0

    def break_page(self):
        if self.fig is not None and self.used > 1e-9:
            self.used = 1.0

    def new_gene(self, header_fn, top_margin=MARGIN_T):
        self.flush()
        self.header_fn, self.part, self.used = header_fn, 0, 1.0
        self.top_margin = top_margin

    def next_row(self, kind="quarter"):
        frac = ROW_FRAC[kind] if isinstance(kind, str) else float(kind)
        if self.fig is None or self.used + frac > 1.0 + 1e-9:
            self._new_page()
        top_margin = getattr(self, "top_margin", MARGIN_T)
        avail = PAGE_H - top_margin - MARGIN_B
        top = PAGE_H - top_margin - self.used * avail
        self.used += frac
        return self.fig, top, frac * avail

    def flush(self):
        if self.fig is not None:
            self.pdf.savefig(self.fig, dpi=RASTER_DPI)
            _flush_pdf_images(self.pdf)
            plt.close(self.fig)
            self.fig = None


def fit_text(text, width_in, size, weight="normal", family=None):
    prop = FontProperties(size=size, weight=weight, family=family)

    def width(t):
        return TextPath((0, 0), t, size=size, prop=prop).get_extents().width / 72.0

    if not text or width(text) <= width_in:
        return text
    lo, hi = 0, len(text)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if width(text[:mid].rstrip() + "...") <= width_in:
            lo = mid
        else:
            hi = mid - 1
    return text[:lo].rstrip() + "..."


def _ax(fig, left, bottom, width, height):
    return fig.add_axes([left / PAGE_W, bottom / PAGE_H, width / PAGE_W, height / PAGE_H])


def panel_columns(n, per_row):
    total = PAGE_W - MARGIN_L - MARGIN_R
    w = (total - PANEL_GAP * (per_row - 1)) / per_row
    return [(MARGIN_L + i * (w + PANEL_GAP), w) for i in range(n)]


_REF_CACHE = {}


def _fetch_ref(fasta, chrom, start, end):
    key = (id(fasta), chrom, start, end)
    if key not in _REF_CACHE:
        if len(_REF_CACHE) > 64:
            _REF_CACHE.clear()
        _REF_CACHE[key] = fasta.fetch(chrom, start, end)
    return _REF_CACHE[key]


def alignment_panel(fig, left, width, top, h, title, bam, chrom, start, end, fasta, tx, args, seed,
                    label_tracks=True, arcs_junctions=None, ref_track=False, highlight=None,
                    marks=None, spans=None, track_label="", show_tx=True, title_width=None, keep_indels=()):
    title = fit_text(title, title_width or width, 6.5 * FS, "bold")
    title_h, tick_h, cov_h, ref_h, pad = 0.19, 0.13, 0.36, (0.15 if ref_track else 0.0), 0.03
    ann_h = 0.20 if show_tx else 0.0
    aln_h = h - title_h - tick_h - cov_h - ann_h - ref_h - (4 if show_tx else 3) * pad - 0.06
    y = top - title_h
    if title:
        fig.text(left / PAGE_W, (y + 0.02) / PAGE_H, title, fontsize=6.5 * FS, fontweight="bold", va="bottom", ha="left")
    y -= tick_h
    ref_seq = _fetch_ref(fasta, chrom, max(0, start), end) if fasta else ""
    if len(ref_seq) < end - max(0, start):
        ref_seq = ref_seq + "N" * (end - max(0, start) - len(ref_seq))
    start = max(0, start)
    acgt = coverage_counts(bam, chrom, start, end)
    reads, n, counts = sample_reads(bam, chrom, start, end, args.max_reads, seed, arcs_junctions)
    arcs = None
    if arcs_junctions:
        arcs = [(_jxn_coords(j)[1] - 2, _jxn_coords(j)[2], counts[j]) for j in arcs_junctions]
    ax = _ax(fig, left, y - cov_h, width, cov_h)
    draw_coverage(ax, acgt, ref_seq, start, end, f"{track_label}\ncoverage" if label_tracks else "", arcs=arcs)
    y -= cov_h + pad
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


REP_MIN_JUNCTIONS = 3
REP_DEPTH_MARGIN = 0.02


def load_gene_psi(manifests, gene, coverage_threshold):
    rows = []
    for stype, mpath in manifests.items():
        if not os.path.isfile(mpath):
            continue
        for rp in _manifest_paths(mpath, gene):
            rp = str(rp).strip()
            if not rp or rp in ("None", "nan", "NA") or not os.path.isfile(rp):
                continue
            df = _result_table(rp)
            need = {"sample", "phasing", "junction", "rescaled_junction_PSI", "junction_coverage"}
            if not need <= set(df.columns):
                continue
            df = df[df["phasing"] == "bulk"]
            rows.append(pd.DataFrame({"sample": df["sample"].astype(str), "junction": df["junction"],
                                      "psi": pd.to_numeric(df["rescaled_junction_PSI"], errors="coerce"),
                                      "cov": pd.to_numeric(df["junction_coverage"], errors="coerce")}))
    if not rows:
        return pd.DataFrame(), pd.DataFrame()
    long = pd.concat(rows, ignore_index=True).drop_duplicates(["sample", "junction"])
    cov = long.pivot(index="sample", columns="junction", values="cov")
    ok = long[(long["cov"] >= coverage_threshold) & long["psi"].notna()]
    psi = ok.pivot(index="sample", columns="junction", values="psi")
    return psi, cov


def junctions_in_window(junctions, start, end):
    out = []
    for j in junctions:
        try:
            _, s1, s2 = _jxn_coords(j)
        except (ValueError, IndexError):
            continue
        if s1 > start and s2 <= end:
            out.append(j)
    return out


def choose_representative(gene, sample, samples, sample_to_type, sample_bams, psi, cov, exclude, junctions=None):
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
    if junctions is not None:
        psi = psi[[j for j in psi.columns if j in set(junctions)]]
        if psi.shape[1] == 0:
            return None
    cands = [s for s in cands if s in psi.index and psi.loc[s].notna().any()]
    if not cands:
        return None
    sub = psi.loc[cands].dropna(axis=1, how="all")
    median = sub.median(axis=0, skipna=True)
    scored = []
    for min_shared in dict.fromkeys((min(REP_MIN_JUNCTIONS, sub.shape[1]), 1)):
        for c in cands:
            diff = (sub.loc[c] - median).abs().dropna()
            if len(diff) >= min_shared:
                scored.append((float(diff.mean()), c))
        if scored:
            break
    if not scored:
        return None
    if cov is None or cov.empty:
        depth = pd.Series(dtype=float)
    else:
        cols = list(sub.columns) if junctions is not None else list(cov.columns)
        depth = cov.reindex(columns=cols).sum(axis=1, min_count=1)
    best_d = min(d for d, _ in scored)
    close = [(d, c) for d, c in scored if d <= best_d + REP_DEPTH_MARGIN]
    d, best = max(close, key=lambda dc: (float(depth.get(dc[1], 0) or 0), -dc[0]))
    return dict(sample=best, sample_type=rtype, dist=d, n=len(cands), n_junctions=int(sub.shape[1]),
                scope="window" if junctions is not None else "gene-wide")


def sashimi_data(bam, chrom, start, end, keep, min_frac):
    L = end - start
    jx = defaultdict(int)
    los, his = [], []
    n = 0
    for r in bam.fetch(chrom, start, end):
        if not keep(r):
            continue
        cig = r.cigartuples
        if cig is None:
            continue
        n += 1
        pos = r.reference_start
        for op, ln in cig:
            if op in (0, 2, 7, 8):
                lo, hi = max(pos, start) - start, min(pos + ln, end) - start
                if lo < hi:
                    los.append(lo)
                    his.append(hi)
                pos += ln
            elif op == 3:
                ss1, ss2 = pos + 1, pos + ln
                if ss1 > start and ss2 <= end:
                    jx[f"{chrom}_{ss1}_{ss2}"] += 1
                pos += ln
    if los:
        cov = np.cumsum(np.bincount(los, minlength=L + 1) - np.bincount(his, minlength=L + 1))[:L].astype(np.int64)
    else:
        cov = np.zeros(L, dtype=np.int64)
    thr = max(1, int(np.floor(min_frac * cov.max()))) if L else 1
    return cov, {j: c for j, c in jx.items() if c >= thr}, n


def intron_transform(start, end, exons, scale):
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
    keep = _keep_spliced
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
        fig.text(ax_left / PAGE_W, (y - lab_h + 0.02) / PAGE_H, t["label"], fontsize=5.5 * FS, fontfamily=WIDE_FONT,
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
            i1, i2 = s1 - 2 - start, s2 - start
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
            peak = 0.125 * pts[0][1] + 0.375 * pts[1][1] + 0.375 * pts[2][1] + 0.125 * pts[3][1]
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


QC_SAMPLE_COLOR = "#333333"
_BLACK_LINE = dict(color="black", linewidth=0.8)
BOX_STYLE = dict(boxprops=_BLACK_LINE, whiskerprops=_BLACK_LINE, capprops=_BLACK_LINE, medianprops=_BLACK_LINE)


def darker(color, f=0.5):
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


def _qc_unavailable(ax, title, msg="not available", size=7):
    ax.set_title(title, fontsize=size * FS, fontweight="bold")
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, color="#777")
    ax.set_xticks([]); ax.set_yticks([])
    _qc_spines(ax)


def _xlabels(ax, labels):
    ax.set_xticklabels(labels, fontsize=5 * FS, rotation=30, ha="right", rotation_mode="anchor", fontfamily=WIDE_FONT)


def qc_read_counts(ax, df, sample, shown, tgroups):
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
            bars.append((f"{t} median", sub["mapped"].median(), sub["on_target_rate"].median() / 100, col))
    x = np.arange(len(bars))
    mapped = np.array([b[1] for b in bars], dtype=float)
    rate = np.array([b[2] for b in bars], dtype=float)
    scale, unit = (1e6, "M") if np.nanmax(mapped) >= 1e6 else (1e3, "k")
    on, off = mapped * rate / scale, mapped * (1 - rate) / scale
    cols = [b[3] for b in bars]
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
    ax.set_ylabel(f"Mapped reads ({unit})", fontfamily=WIDE_FONT)
    ax.set_ylim(0, np.nanmax(mapped) / scale * 1.4)
    ax.set_title(title, fontsize=7 * FS, fontweight="bold", fontfamily=WIDE_FONT)
    ax.legend(handles=legend_handles, prop=FontProperties(family=WIDE_FONT, size=5 * FS), frameon=False, loc="upper right", ncol=2,
              bbox_to_anchor=(1, 1.03))
    _qc_spines(ax)

def qc_read_lengths(ax, df, sample, shown, tgroups):
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
    b = ax.bxp(stats, positions=pos, widths=0.6, showfliers=False, patch_artist=True,
               boxprops=dict(edgecolor="black", linewidth=0.8), whiskerprops=_BLACK_LINE, capprops=_BLACK_LINE,
               medianprops=_BLACK_LINE)
    for patch, c in zip(b["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("black")
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], fontsize=(5 if n_box <= 3 else 4.2) * FS, rotation=30, ha="right",
                       rotation_mode="anchor", fontfamily=WIDE_FONT)
    ax.set_xlim(-0.6, 2 * n_box + 0.6)
    for xc, lab in groups:
        ax.text(xc, 1.0, lab, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=6 * FS, fontweight="bold")
    ax.set_ylabel("Read length (bp)", fontfamily=WIDE_FONT)
    ax.set_title(title, fontsize=7 * FS, fontweight="bold", pad=14 * FS, fontfamily=WIDE_FONT)
    _qc_spines(ax)

def qc_pca(ax, df, sample, shown):
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
    ax.set_xlabel(f"PC1 ({v1:.1f}%)", fontfamily=WIDE_FONT); ax.set_ylabel(f"PC2 ({v2:.1f}%)", fontfamily=WIDE_FONT)
    ax.set_title(title, fontsize=7 * FS, fontweight="bold", fontfamily=WIDE_FONT)
    _qc_spines(ax)
    refs = df[df["kind"] == "reference"]
    return [(str(g), gs["color"].iloc[0], len(gs)) for g, gs in refs.groupby("group", sort=False)]


POINT_ALPHA = 0.5
POINT_SIZE = 8
STAR_SIZE = 100
FRACTION_YLIM = (-0.03, 1.05)


def _qc_box(ax, entries, title, ylabel=None, ylim=FRACTION_YLIM):
    rng = np.random.default_rng(0)
    pos, data, shown_vals = [], [], []
    for i, (_lab, v, col, star) in enumerate(entries):
        v = np.asarray(v, dtype=float)
        v = v[~np.isnan(v)]
        if len(v):
            ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=POINT_SIZE if len(v) <= 50 else POINT_SIZE / 2,
                       color=col, alpha=POINT_ALPHA, linewidths=0, zorder=1)
            pos.append(i)
            data.append(v)
            shown_vals.append(v)
        if star is not None and not np.isnan(star):
            ax.scatter([i], [star], s=STAR_SIZE, marker="*", color=QC_SAMPLE_COLOR, edgecolors="black",
                       linewidths=0.5, zorder=5)
            shown_vals.append(np.array([star]))
    if data:
        ax.boxplot(data, positions=pos, widths=0.6, showfliers=False, zorder=2, **BOX_STYLE)
    ax.set_xticks(range(len(entries)))
    ax.set_xticklabels([e[0] for e in entries], fontsize=5 * FS, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.6, len(entries) - 0.4)
    if ylim is not None:
        ax.set_ylim(*ylim)
        if ylim == FRACTION_YLIM:
            ax.set_yticks([0, 0.5, 1.0])
            ax.set_yticklabels(["0.0", "0.5", "1.0"])
    elif shown_vals:
        vals = np.concatenate(shown_vals)
        lo, hi = vals.min(), vals.max()
        padv = max((hi - lo) * 0.15, 0.02)
        ax.set_ylim(max(0, lo - padv), min(1.0, hi + padv) if hi <= 1 else hi + padv)
    ax.tick_params(axis="y", labelsize=5 * FS)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=5 * FS)
    ax.set_title(title, fontsize=6 * FS, fontweight="bold")
    _qc_spines(ax)


def sample_box(ax, values, sample, shown, tgroups, title, ylabel=None, ylim=FRACTION_YLIM):
    values = pd.to_numeric(pd.Series(values), errors="coerce")
    values.index = values.index.astype(str)
    values = values[~values.index.duplicated()]
    own = float(values.get(sample, np.nan))
    entries = [(shown, [], QC_SAMPLE_COLOR, own)]
    entries += [(t, values.reindex([m_ for m_ in members if m_ in values.index]).to_numpy(dtype=float), col, None)
                for t, members, col in tgroups]
    _qc_box(ax, entries, title, ylabel=ylabel, ylim=ylim)


def draw_expression_box(ax, gene, sample, shown, tgroups, expr):
    title = "Expression"
    rows = [m.loc[gene] for m in expr.values() if m is not None and gene in m.index]
    rows = [r.iloc[0] if isinstance(r, pd.DataFrame) else r for r in rows]
    vals = pd.to_numeric(pd.concat(rows), errors="coerce") if rows else pd.Series(dtype=float)
    vals.index = vals.index.astype(str)
    vals = vals[~vals.index.duplicated()]
    if vals.empty or np.isnan(vals.get(sample, np.nan)):
        return _qc_unavailable(ax, title, "no expression data", size=6)
    sample_box(ax, np.log2(vals + 1), sample, shown, tgroups, title, ylabel="log2(MOTR + 1)", ylim=None)


def gene_metric_box(ax, matrix, gene, sample, shown, tgroups, title, ylabel):
    if matrix is None or gene not in matrix.index or sample not in matrix.columns:
        return _qc_unavailable(ax, title, "no data", size=6)
    row = matrix.loc[gene]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    sample_box(ax, row, sample, shown, tgroups, title, ylabel=ylabel)


def _ranked(values, grid):
    v = np.sort(np.asarray(values, dtype=float)[~np.isnan(np.asarray(values, dtype=float))])[::-1]
    if len(v) == 0:
        return np.full(len(grid), np.nan)
    if len(v) == 1:
        return np.full(len(grid), v[0])
    return np.interp(grid, np.linspace(0, 100, len(v)), v)


def qc_curve(ax, matrix, sample, shown, tgroups, ylim=(-0.02, 1.02), ylabel=None):
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
    ax.set_xlabel("Genes (ranked)", fontsize=5 * FS)
    ax.set_xticks([])
    ax.tick_params(labelsize=5 * FS)
    ax.set_title("ranked genes", fontsize=6 * FS, fontweight="bold")
    _qc_spines(ax)


def qc_mean_panel(ax, matrix, sample, shown, tgroups, title):
    if matrix is None or matrix.empty or sample not in matrix.columns:
        return _qc_unavailable(ax, title)
    sample_box(ax, matrix.mean(axis=0, skipna=True), sample, shown, tgroups, title)


def load_flr_matrix(path):
    df = _read_tsv(path, index_col=0)
    if df is None:
        return None
    df.columns = df.columns.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def load_testability_matrix(all_rows):
    if all_rows is None:
        return None
    cols = {"sample", "gene", "n_canonical_jxns", "n_canonical_jxns_tested"}
    df = all_rows[[c for c in all_rows.columns if c in cols]].copy()
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


_POSSIBLE_CELLS = frozenset(tier_rules.possible_cells())


def draw_hit_tier_grid(fig, rows, left=1.3, bottom=0.35, cw=0.62, rh=0.32, max_w=None):
    cats = ["variant", "ASE", "junction", "isoform"]
    order = ["variant", "junction", "isoform", "ASE"]
    combos = [frozenset(c) for k in range(1, len(order) + 1) for c in itertools.combinations(order, k)]
    if max_w:
        cw = min(cw, max_w / len(combos))
    tiers = [1, 2, 3, 4, 5, 6]
    grid = np.zeros((len(tiers), len(combos)), dtype=int)
    no_cat = 0
    for _, r in rows.iterrows():
        have = frozenset(c for c, ok in (
            ("variant", _hit_flag(r.get("variant"))),
            ("ASE", _hit_flag(r.get("ASE"))),
            ("junction", _hit_flag(r.get("outlier_junction")) or _hit_flag(r.get("cohort_outlier_junction"))),
            ("isoform", _hit_flag(r.get("isoform_outlier"))),
        ) if ok)
        tier = pd.to_numeric(pd.Series([r.get("tier")]), errors="coerce").iloc[0]
        if pd.isna(tier) or int(tier) not in tiers:
            continue
        if not have:
            no_cat += 1
            continue
        grid[tiers.index(int(tier)), combos.index(have)] += 1
    if no_cat:
        print(f"  [WARNING] {no_cat} candidate gene(s) have a tier but no variant / ASE / junction / isoform hit; "
              "left out of the hit-category grid.")
    n_genes = int(grid.sum())

    gw, gh = cw * len(combos), rh * len(tiers)
    grid_bottom = bottom + 1.0
    ax = _ax(fig, left, grid_bottom, gw, gh)
    vmax = max(int(grid.max()), 1)
    cmap = LinearSegmentedColormap.from_list("hits", ["#ffffff", QC_SAMPLE_COLOR])
    for i in range(len(tiers)):
        for j in range(len(combos)):
            v = grid[i, j]
            impossible = (combos[j], tiers[i]) not in _POSSIBLE_CELLS
            frac = np.sqrt(v / vmax)
            if impossible and not v:
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="#eeeeee", edgecolor="#dddddd", linewidth=0.5))
                ax.text(j + 0.5, i + 0.5, "\u2013", ha="center", va="center", fontsize=6.5 * FS, color="#aaaaaa")
                continue
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(frac), edgecolor="#dddddd", linewidth=0.5))
            if v:
                ax.text(j + 0.5, i + 0.5, str(v), ha="center", va="center", fontsize=6.5 * FS,
                        color="white" if frac > 0.6 else "black")
            if impossible and v:
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="none", edgecolor="#c23637", linewidth=1.2))
                print(f"  [WARNING] {v} gene(s) in tier {tiers[i]} with hit categories "
                      f"{'+'.join(sorted(combos[j]))}, which the tier rules should not produce.")
    ax.set_xlim(0, len(combos)); ax.set_ylim(len(tiers), 0)
    ax.set_xticks([]); ax.set_yticks(np.arange(len(tiers)) + 0.5)
    ax.set_yticklabels([f"Tier {t}" for t in tiers], fontsize=6 * FS, fontfamily=WIDE_FONT)
    ax.tick_params(axis="y", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    col_tot, row_tot = grid.sum(axis=0), grid.sum(axis=1)
    axt = _ax(fig, left, grid_bottom + gh + 0.08, gw, 0.7)
    axt.bar(np.arange(len(combos)) + 0.5, col_tot, width=0.6, color=QC_SAMPLE_COLOR)
    for j, v in enumerate(col_tot):
        axt.text(j + 0.5, v, str(v), ha="center", va="bottom", fontsize=5.5 * FS)
    axt.set_xlim(0, len(combos)); axt.set_ylim(0, max(col_tot.max(), 1) * 1.3)
    axt.axis("off")
    axt.set_title(f"Candidate genes by hit category and tier ({n_genes} genes)", fontsize=7 * FS,
                  fontweight="bold", pad=4, fontfamily=WIDE_FONT)
    axr = _ax(fig, left + gw + 0.08, grid_bottom, 0.9, gh)
    axr.barh(np.arange(len(tiers)) + 0.5, row_tot, height=0.6, color=QC_SAMPLE_COLOR)
    for i, v in enumerate(row_tot):
        axr.text(v, i + 0.5, f" {v}", ha="left", va="center", fontsize=5.5 * FS)
    axr.set_ylim(len(tiers), 0); axr.set_xlim(0, max(row_tot.max(), 1) * 1.35)
    axr.axis("off")
    axm = _ax(fig, left, bottom, gw, 0.85)
    cat_y = {c: len(cats) - k for k, c in enumerate(cats)}
    for j, c in enumerate(combos):
        ys = [cat_y[k] for k in c]
        if len(ys) > 1:
            axm.plot([j + 0.5, j + 0.5], [min(ys), max(ys)], color="black", linewidth=1.5, zorder=2)
        for k in cats:
            axm.scatter([j + 0.5], [cat_y[k]], s=40, zorder=3, color="black" if k in c else "#dddddd")
    axm.set_xlim(0, len(combos)); axm.set_ylim(0.5, len(cats) + 0.5)
    axm.set_yticks(list(cat_y.values())); axm.set_yticklabels(list(cat_y.keys()), fontsize=6 * FS, fontfamily=WIDE_FONT)
    axm.tick_params(axis="y", length=0)
    axm.set_xticks([])
    for sp in axm.spines.values():
        sp.set_visible(False)


def draw_page_legend(fig, left, top, shown, tgroups, refs=None):
    width, height = PAGE_W - MARGIN_R - left, 3.0
    ax = _ax(fig, left, top - height, width, height)
    ax.set_xlim(0, width); ax.set_ylim(top - height, top); ax.axis("off")
    mx, tx, step = 0.15, 0.3, 0.16
    y = top
    x0 = 0.08
    ax.text(x0, y - 0.06, "TEQUILA-Dx", fontsize=6 * FS, fontweight="bold", va="top", fontfamily=WIDE_FONT)
    y -= 0.06
    y -= step * 1.5
    entries = [("star", QC_SAMPLE_COLOR, shown)] + [("diamond", col, f"{t} (n={len(m_)})") for t, m_, col in tgroups]
    if refs:
        entries.append(("head", None, "GTEx"))
        entries += [("dot", col, textwrap.shorten(str(g), 24, placeholder="...")) for g, col, _ in refs]
    for kind, col, lab in entries:
        if kind == "head":
            y -= step * 0.3
            ax.text(0.08, y, lab, fontsize=6 * FS, fontweight="bold", va="center", fontfamily=WIDE_FONT)
        else:
            if kind == "star":
                ax.scatter([mx], [y], marker="*", s=80, color=col, edgecolors="black", linewidths=0.5, clip_on=False)
            elif kind == "diamond":
                ax.scatter([mx], [y], marker="D", s=30, color=col, edgecolors="white", linewidths=0.3, clip_on=False)
            else:
                ax.scatter([mx], [y], marker="o", s=14, color=col, alpha=0.5, linewidths=0, clip_on=False)
            ax.text(tx, y, lab, fontsize=5.5 * FS, va="center", fontfamily=WIDE_FONT)
        y -= step
    wide = FontProperties(family=WIDE_FONT, size=5.5 * FS)
    bold = FontProperties(family=WIDE_FONT, size=6 * FS, weight="bold")
    longest = max([TextPath((0, 0), e[2], prop=wide).get_extents().width / 72.0 + tx for e in entries if e[0] != "head"]
                  + [TextPath((0, 0), t, prop=bold).get_extents().width / 72.0 + x0
                     for t in ["TEQUILA-Dx"] + [e[2] for e in entries if e[0] == "head"]])
    box_w = min(width, longest + 0.12)
    ax.add_patch(FancyBboxPatch((0, y + step * 0.4), box_w, top - (y + step * 0.4),
                                boxstyle="round,pad=0,rounding_size=0.05", facecolor="none",
                                edgecolor="#999999", linewidth=0.6, clip_on=False))
    return y + step * 0.4


def _wrap_to_width(text, prop, width_in):
    lines, cur = [], ""
    for word in text.split():
        trial = f"{cur} {word}".strip()
        if cur and TextPath((0, 0), trial, prop=prop).get_extents().width / 72.0 > width_in:
            lines.append(cur)
            cur = word
        else:
            cur = trial
    if cur:
        lines.append(cur)
    return lines


def draw_hpo_block(fig, left, top, bottom, labels):
    width = PAGE_W - MARGIN_R - left
    body = FontProperties(family=WIDE_FONT, size=5.5 * FS)
    step = 0.132
    fig.text(left / PAGE_W, top / PAGE_H, "Phenotype (HPO)", fontsize=6 * FS, fontweight="bold", va="top",
             fontfamily=WIDE_FONT)
    y = top - step * 1.3
    for i, label in enumerate(labels):
        lines = _wrap_to_width(label, body, width - 0.12)
        if y - step * len(lines) < bottom:
            fig.text(left / PAGE_W, y / PAGE_H, f"(+{len(labels) - i} more)", fontsize=5.5 * FS, va="top",
                     color="#777", fontfamily=WIDE_FONT)
            return
        for k, line in enumerate(lines):
            fig.text((left + (0 if k == 0 else 0.12)) / PAGE_W, y / PAGE_H, line, fontsize=5.5 * FS, va="top",
                     color="#333", fontfamily=WIDE_FONT)
            y -= step


def wgs_concordance_text(path, sample):
    df = _read_tsv(path, dtype=str, keep_default_na=False)
    if df is None or "sample" not in df.columns:
        return None, None
    row = df[df["sample"] == str(sample)]
    if row.empty:
        return None, None
    r = row.iloc[0]
    pct = lambda v: f"{float(v) * 100:.0f}%"
    n = r["n_rna_snvs"]
    flag = r["flag"]
    if flag == "ok":
        return f"WGS concordance: {pct(r['own_wgs_concordance'])} of {n} RNA SNVs found in own WGS", "#333"
    if flag == "mismatch":
        return (f"WGS MISMATCH: best match is {r['best_wgs_match']}'s WGS ({pct(r['best_wgs_concordance'])}); "
                f"own WGS {pct(r['own_wgs_concordance'])} of {n} RNA SNVs"), "#c62828"
    if flag == "low":
        return f"WGS concordance LOW: {pct(r['own_wgs_concordance'])} of {n} RNA SNVs found in own WGS", "#c62828"
    if flag == "matches_other_wgs":
        return (f"No WGS, but RNA matches {r['best_wgs_match']}'s WGS "
                f"({pct(r['best_wgs_concordance'])} of {n} RNA SNVs)"), "#c62828"
    if flag == "too_few_sites" and r["has_wgs"] == "True":
        return f"WGS concordance: not assessed (only {n} high-confidence RNA SNVs)", "#777"
    return None, None


def draw_qc_page(pdf, args, sample, shown, sample_to_type, type_colors, hit_rows, all_rows=None):
    stype = sample_to_type.get(sample, "unknown")
    global QC_SAMPLE_COLOR
    QC_SAMPLE_COLOR = darker(type_colors.get(stype, "#4c72b0"))
    global ARC_COLOR
    ARC_COLOR = QC_SAMPLE_COLOR
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
             fontsize=11 * FS, fontweight="bold", va="top", fontfamily=WIDE_FONT)
    fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H,
             ", ".join(f"{t} (n = {len(m_)}{' others' if t == stype else ''})" for t, m_, _ in tgroups),
             fontsize=6.5 * FS, va="top", color="#333", fontfamily=WIDE_FONT)
    concord_text, concord_color = wgs_concordance_text(args.wgs_concordance_tsv, sample)
    if concord_text:
        fig.text((PAGE_W - MARGIN_R) / PAGE_W, 1 - 0.48 / PAGE_H, concord_text, ha="right",
                 fontsize=6.5 * FS, va="top", color=concord_color, fontfamily=WIDE_FONT,
                 fontweight="bold" if concord_color == "#c62828" else "normal")
    qc_read_counts(_ax(fig, 1.05, 5.05, 1.95, 2.45), _read_tsv(args.on_target_tsv), sample, shown, tgroups)
    qc_read_lengths(_ax(fig, 3.65, 5.05, 2.2, 2.25), _read_tsv(args.read_attributes_tsv), sample, shown, tgroups)
    refs = qc_pca(_ax(fig, 6.45, 5.05, 2.0, 2.45), _read_tsv(args.pca_coords), sample, shown)
    legend_bottom = draw_page_legend(fig, 8.85, 7.6, shown, tgroups, refs)
    hpo_labels = load_hpo_labels(args.hpo_terms_tsv)
    if hpo_labels:
        draw_hpo_block(fig, 8.85, legend_bottom - 0.2, 4.45, hpo_labels)
    flr = load_flr_matrix(args.flr_matrix)
    testability = load_testability_matrix(all_rows)
    coverage = load_flr_matrix(args.coverage_matrix)
    canonical_bases, canonical_jxns = {}, {}
    sample_cov = _read_tsv(args.sample_coverage_tsv, dtype=str, keep_default_na=False)
    if sample_cov is not None and {"gene", "length"} <= set(sample_cov.columns):
        canonical_bases = _int_lookup(sample_cov["gene"], sample_cov["length"])
    if all_rows is not None and {"gene", "sample", "n_canonical_jxns"} <= set(all_rows.columns):
        own = all_rows[all_rows["sample"].astype(str) == str(sample)]
        rest = all_rows[all_rows["sample"].astype(str) != str(sample)]
        canonical_jxns = {**_int_lookup(rest["gene"], rest["n_canonical_jxns"]),
                          **_int_lookup(own["gene"], own["n_canonical_jxns"])}

    draw_hit_tier_grid(fig, hit_rows, left=1.05, bottom=0.3, cw=0.5, rh=0.28, max_w=4.3)
    for k, (mat, name, ylab) in enumerate((
            (coverage, "Coverage", f"Fraction of canonical bases\nwith >= {args.coverage_depth}x coverage"),
            (testability, "Junction testability", "Fraction of canonical junctions tested"),
            (flr, "Full-length ratio", "Average proportion of canonical transcript\ncovered by a single read"))):
        left = 6.55 + k * 1.5
        bax = _ax(fig, left, 2.55, 0.92, 1.45)
        qc_mean_panel(bax, mat, sample, shown, tgroups, f"{name}\n(per-sample mean)")
        bax.title.set_fontfamily(WIDE_FONT)
        bax.set_ylabel(ylab, fontsize=4.5 * FS, labelpad=2, linespacing=1.1)
        fit_ylabel(bax, 4.5 * FS)
        cax = _ax(fig, left, 0.55, 0.92, 1.25)
        if mat is None or mat.empty or sample not in mat.columns:
            _qc_unavailable(cax, "(per gene)", size=6)
        else:
            qc_curve(cax, mat, sample, shown, tgroups, ylabel=ylab)
            cax.set_title("(per gene)", fontsize=6 * FS, fontweight="bold")
            cax.set_xticks([0, 50, 100])
            cax.yaxis.label.set_linespacing(1.1)
            fit_ylabel(cax, 4.5 * FS)
        cax.title.set_fontfamily(WIDE_FONT)
    pdf.savefig(fig, dpi=RASTER_DPI)
    plt.close(fig)

    return dict(stype=stype, tgroups=tgroups, flr=flr, testability=testability, coverage=coverage,
                canonical_bases=canonical_bases, canonical_jxns=canonical_jxns)


OTHER_GENES = "Other panel genes"
OTHER_HITS = "Other gene hits"
OTHER_ISOFORMS = "other"


def load_tx(args, sample):
    tx = dict(unavailable=args.tx_unavailable, usage={}, status=None, iso=pd.DataFrame(), land=pd.DataFrame(),
              args=args)
    if tx["unavailable"]:
        return tx
    usage = _read_tsv(args.tx_usage_matrix, dtype={"gene": str, "transcript_id": str})
    status = _read_tsv(args.tx_gene_status, dtype={"gene": str, "status": str})
    if usage is None or status is None:
        tx["unavailable"] = "the transcript analysis results were not found"
        return tx
    usage.columns = usage.columns.astype(str)
    status.columns = status.columns.astype(str)
    tx["usage"] = {g: d.drop(columns="gene").set_index("transcript_id")
                   for g, d in usage.groupby("gene", sort=False)}
    tx["status"] = status.drop_duplicates("gene").set_index("gene")
    for key, path in (("iso", args.tx_isoform_outliers), ("land", args.tx_landscape_outliers)):
        df = _read_tsv(path, dtype={"sample": str, "gene": str, "transcript_id": str})
        if df is not None and not df.empty:
            tx[key] = df[df["sample"].astype(str) == str(sample)].reset_index(drop=True)
    return tx


def tx_outlier_genes(tx):
    return {str(g) for key in ("iso", "land") if not tx[key].empty for g in tx[key]["gene"]}


def tx_reason(tx, gene, sample, stype):
    if tx["unavailable"]:
        return "not available: " + tx["unavailable"]
    a = tx["args"]
    if gene not in tx["status"].index:
        return "not available: AMALGAM quantified no transcripts for this gene"
    row = tx["status"].loc[gene]
    status = str(row["status"])
    if status == "low_n":
        return (f"not available: only {int(row['n_samples_used'])} {stype} sample(s) had >= {a.tx_min_gene_reads:g} "
                f"reads on this gene ({a.tx_n_threshold} needed)")
    if status == "no_isoform":
        return (f"not available: no isoform reached >= {a.tx_min_usage:.0%} usage with >= {a.tx_min_reads:g} "
                f"reads in any {stype} sample")
    if status == "single_isoform":
        return "not available: only one isoform of this gene is expressed"
    u = tx["usage"].get(gene)
    if u is None or sample not in u.columns or u[sample].isna().all():
        reads = pd.to_numeric(row.get(sample, np.nan), errors="coerce")
        if pd.isna(reads):
            return "not available: this sample was not in the transcript analysis"
        return f"not available: this sample has {reads:.0f} reads on this gene (fewer than {a.tx_min_gene_reads:g})"
    return None


def _tx_stat(rec):
    if _nonempty(rec.get("padj")):
        return f"padj {float(rec['padj']):.1e}"
    if _nonempty(rec.get("modz")):
        return f"z {float(rec['modz']):+.1f}"
    return ""


def tx_gene_outliers(tx, gene):
    iso = [r for r in tx["iso"].to_dict("records") if str(r["gene"]) == gene] if not tx["iso"].empty else []
    land = [r for r in tx["land"].to_dict("records") if str(r["gene"]) == gene] if not tx["land"].empty else []
    return iso, land


def _lighter_colors(n):
    tiers = len(PALETTE_ORDER) // 7
    return [PALETTE_ORDER[(i % 7) + 7 * (1 + (i // 7) % (tiers - 1))] for i in range(n)]


def _cluster_order(u, main):
    cols = list(u.columns)
    if len(cols) < 3:
        return sorted(cols, key=lambda c: (-u.loc[main, c], c))
    x = np.sqrt(np.clip(u.fillna(0).to_numpy(dtype=float).T, 0, None))
    order = [cols[i] for i in leaves_list(linkage(x, method="average", metric="euclidean", optimal_ordering=True))]
    if u.loc[main, order[0]] < u.loc[main, order[-1]]:
        order = order[::-1]
    return order


def draw_isoform_stack(ax, gene, sample, shown, stype, tx, out_tx):
    u = tx["usage"][gene]
    u = u.loc[:, [c for c in u.columns if u[c].notna().any()]]
    names = sorted((n for n in u.index if n != OTHER_ISOFORMS), key=lambda n: -u.loc[n].mean())
    names += [OTHER_ISOFORMS] if OTHER_ISOFORMS in u.index else []
    u = u.loc[names]
    order = _cluster_order(u, names[0])
    cols = _lighter_colors(len(names))
    colors = {n: (PALETTE_LIGHT_GREY if n == OTHER_ISOFORMS else c) for n, c in zip(names, cols)}
    x = np.arange(len(order))
    width = 0.85 if len(order) <= 60 else 1.0
    bottom = np.zeros(len(order))
    for n in names:
        v = u.loc[n, order].to_numpy(dtype=float)
        ax.bar(x, v, bottom=bottom, width=width, color=colors[n], edgecolor="white",
               linewidth=0.3 if len(order) <= 60 else 0, zorder=2)
        bottom += v
    i = order.index(sample)
    ax.bar([i], [1.0], width=width, facecolor="none", edgecolor="black", linewidth=1.0, zorder=3)
    ax.scatter([i], [1.07], marker="*", s=STAR_SIZE, color=QC_SAMPLE_COLOR, edgecolors="black", linewidths=0.5,
               zorder=5, clip_on=False)
    ax.set_xticks([i])
    ax.set_xticklabels([shown], fontsize=5 * FS, fontweight="bold", color=QC_SAMPLE_COLOR, fontfamily=WIDE_FONT)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0.0", "0.5", "1.0"])
    ax.tick_params(axis="y", labelsize=5 * FS)
    ax.set_ylabel("Isoform usage", fontsize=5 * FS)
    ax.set_xlabel(f"{stype} samples (n = {len(order)}), clustered by isoform mix", fontsize=5 * FS,
                  fontfamily=WIDE_FONT)
    _qc_spines(ax)
    return colors


def add_isoform_legend(ax, colors, out_tx):
    names = list(colors)
    handles = [Patch(facecolor=colors[n], edgecolor="none", label=n + ("  (outlier)" if n in out_tx else ""))
               for n in names]
    leg = ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                    handlelength=1.0, handleheight=1.0, borderaxespad=0, fontsize=5 * FS)
    for t, n in zip(leg.get_texts(), names):
        t.set_fontfamily(WIDE_FONT)
        if n in out_tx:
            t.set_fontweight("bold")
            t.set_color(QC_SAMPLE_COLOR)


_GTF_TID = re.compile(r'transcript_id "([^"]+)"')
_STRUCT_CACHE = {}


def load_tx_structures(gtf, region, tids):
    key = (gtf, region, tuple(sorted(tids)))
    if key in _STRUCT_CACHE:
        return _STRUCT_CACHE[key]
    out = {}
    if not gtf or not os.path.isfile(gtf) or region is None:
        return out
    chrom, start, end = region
    try:
        tbx = pysam.TabixFile(gtf)
    except (OSError, ValueError):
        return out
    contigs = set(tbx.contigs)
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    contig = chrom if chrom in contigs else (alt if alt in contigs else None)
    if contig is None:
        tbx.close()
        return out
    want = set(tids)

    def scan(s, e):
        recs = defaultdict(lambda: {"exons": [], "span": [None, None], "strand": "+"})
        for line in tbx.fetch(contig, max(0, s), e):
            p = line.split("\t")
            if len(p) < 9:
                continue
            m = _GTF_TID.search(p[8])
            if not m or m.group(1) not in want:
                continue
            r = recs[m.group(1)]
            a, b = int(p[3]) - 1, int(p[4])
            r["strand"] = p[6]
            r["span"] = [a if r["span"][0] is None else min(r["span"][0], a),
                         b if r["span"][1] is None else max(r["span"][1], b)]
            if p[2] == "exon":
                r["exons"].append((a, b))
        return recs

    first = scan(start, end)
    spans = [r["span"] for r in first.values() if r["span"][0] is not None]
    if spans:
        recs = scan(min(s for s, _ in spans), max(e for _, e in spans))
        for tid, r in recs.items():
            if r["exons"]:
                exons = sorted(set(r["exons"]))
                out[tid] = dict(exons=exons, strand=r["strand"])
    tbx.close()
    _STRUCT_CACHE[key] = out
    return out


def _compressed_axis(structures, intron_scale=0.25):
    blocks = []
    for a, b in sorted(e for s in structures.values() for e in s["exons"]):
        if blocks and a <= blocks[-1][1]:
            blocks[-1][1] = max(blocks[-1][1], b)
        else:
            blocks.append([a, b])
    xs, ys, pos = [], [], 0.0
    for i, (a, b) in enumerate(blocks):
        if i:
            pos += (a - blocks[i - 1][1]) * intron_scale
        xs += [a, b]
        ys += [pos, pos + (b - a)]
        pos += b - a
    return np.array(xs, dtype=float), np.array(ys, dtype=float), pos


def draw_isoform_structures(fig, left, width, top, h, names, colors, structures, out_tx, other_note=""):
    shown = [n for n in names if n != OTHER_ISOFORMS]
    found = {n: structures[n] for n in shown if n in structures}
    has_other = OTHER_ISOFORMS in names
    title_h = 0.22
    fig.text(left / PAGE_W, (top - 0.15) / PAGE_H,
             "Isoform structures (5\u2032 \u2192 3\u2032; introns shared by every isoform drawn at 25%)",
             fontsize=5.5 * FS, va="bottom", ha="left", color="#333", fontfamily=WIDE_FONT)
    ax = _ax(fig, left, top - h + 0.05, width, h - title_h - 0.05)
    ax.axis("off")
    n = len(shown) + (1 if has_other else 0)
    ax.set_ylim(n - 0.5, -0.5)
    if not found:
        ax.set_xlim(0, 1)
        ax.text(0.5, 0.5, "isoform structures not found in the AMALGAM GTF", transform=ax.transAxes,
                ha="center", va="center", fontsize=5.5 * FS, color="#777")
        return
    xs, ys, total = _compressed_axis(found)
    minus = Counter(s["strand"] for s in found.values()).most_common(1)[0][0] == "-"

    def X(g):
        v = np.interp(g, xs, ys)
        return total - v if minus else v

    ax.set_xlim(-0.01 * total, total * 1.01)
    for i, name in enumerate(shown):
        col = colors.get(name, "#777")
        s = found.get(name)
        if s is None:
            ax.text(0, i, "structure not in GTF", va="center", ha="left", fontsize=5 * FS, color="#777")
        else:
            a, b = X(s["exons"][0][0]), X(s["exons"][-1][1])
            ax.plot([min(a, b), max(a, b)], [i, i], color=col, linewidth=0.8, zorder=1, solid_capstyle="butt")
            for e0, e1 in s["exons"]:
                x0, x1 = sorted((X(e0), X(e1)))
                ax.add_patch(Rectangle((x0, i - 0.32), max(x1 - x0, total * 0.0015), 0.64, facecolor=col,
                                       edgecolor="none", zorder=2))
        ax.text(1.01, i, name + ("  (outlier)" if name in out_tx else ""), transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=5 * FS, color=QC_SAMPLE_COLOR if name in out_tx else "#333",
                fontweight="bold" if name in out_tx else "normal", fontfamily=WIDE_FONT)
    if has_other:
        i = len(shown)
        ax.add_patch(Rectangle((0, i - 0.32), total * 0.03, 0.64, facecolor=colors.get(OTHER_ISOFORMS, "#ccc"),
                               edgecolor="none", zorder=2))
        ax.text(total * 0.04, i, "other: minor isoforms pooled" + (f" ({other_note})" if other_note else ""),
                va="center", ha="left", fontsize=5 * FS, color="#555", fontfamily=WIDE_FONT)
        ax.text(1.01, i, OTHER_ISOFORMS, transform=ax.get_yaxis_transform(), va="center", ha="left",
                fontsize=5 * FS, color="#333", fontfamily=WIDE_FONT)


def draw_isoform_box(ax, values, rec, sample, shown, tgroups):
    sample_box(ax, values, sample, shown, tgroups, str(rec["transcript_id"]), ylabel="Isoform usage")
    stat = _tx_stat(rec)
    ax.text(0.5, 1.0, f"\u0394 = {float(rec['delta_usage']):+.2f}" + (f", {stat}" if stat else ""),
            transform=ax.transAxes, ha="center", va="bottom", fontsize=5 * FS, color=QC_SAMPLE_COLOR,
            fontweight="bold")
    ax.set_title(fit_text(str(rec["transcript_id"]), ax.get_position().width * PAGE_W, 6 * FS, "bold"),
                 fontsize=6 * FS, fontweight="bold", pad=10 * FS, fontfamily=WIDE_FONT)


def draw_isoform_rows(pager, gene, sample, shown, ctx, tx, kind="half"):
    stype = ctx["stype"]
    avail = PAGE_H - MARGIN_T - MARGIN_B
    frac = ROW_FRAC[kind]
    if not tx_reason(tx, gene, sample, stype):
        n_names = len(tx["usage"][gene].index)
        iso_rows = -(-len(tx_gene_outliers(tx, gene)[0]) // DIST_PANELS_PER_ROW)
        need = ((0.3 + 0.15 * n_names) / avail if tx["args"].tx_gtf else 0.0) + (ROW_FRAC["third"] if iso_rows else 0.0)
        left_on_page = 1.0 - pager.used if pager.fig is not None and pager.used < 1.0 else 1.0
        frac = max(0.22, min(frac, left_on_page - need))
    fig, top, h = pager.next_row(frac)
    left, width = MARGIN_L, PAGE_W - MARGIN_L - MARGIN_R
    title_artist = fig.text(left / PAGE_W, (top - 0.17) / PAGE_H, "Isoform usage", fontsize=6.5 * FS,
                            fontweight="bold", va="bottom", ha="left", fontfamily=WIDE_FONT)
    reason = tx_reason(tx, gene, sample, stype)
    if reason:
        fig.text(left / PAGE_W, (top - 0.42) / PAGE_H, reason, fontsize=6 * FS, va="bottom", ha="left",
                 color="#777", fontfamily=WIDE_FONT)
        return
    iso, land = tx_gene_outliers(tx, gene)
    def _recur(r):
        n_out, n_all = r.get("n_outlier_samples"), r.get("n")
        return f"; in {int(n_out)}/{int(n_all)} samples" if _nonempty(n_out) and _nonempty(n_all) else ""

    parts = [f"{r['transcript_id']} usage {float(r['usage']):.2f} (\u0394 {float(r['delta_usage']):+.2f}"
             + (f", {_tx_stat(r)}" if _tx_stat(r) else "") + _recur(r) + ")" for r in iso]
    parts += [f"isoform mix outlier: TVD {float(r['tvd']):.2f} from the typical mix, "
              + (f"{float(r['nearest_tvd']):.2f} from the nearest sample" if _nonempty(r.get("nearest_tvd")) else "")
              + " (" + "; ".join(x for x in (_tx_stat(r), _recur(r)[2:]) if x) + ")" for r in land]
    if parts:
        renderer = fig.canvas.get_renderer()
        x0 = title_artist.get_window_extent(renderer=renderer).x1 / fig.dpi + 0.2
        fig.text(x0 / PAGE_W, (top - 0.17) / PAGE_H,
                 fit_text("  |  ".join(parts), left + width - x0, 5.5 * FS, "bold", WIDE_FONT), fontsize=5.5 * FS,
                 va="bottom",
                 ha="left", color=QC_SAMPLE_COLOR, fontweight="bold", fontfamily=WIDE_FONT)
    out_tx = [str(r["transcript_id"]) for r in iso]
    bar_left, bar_w = left + 0.4, width - 0.4 - 2.2
    ax = _ax(fig, bar_left, top - h + 0.42, bar_w, h - 0.42 - 0.45)
    colors = draw_isoform_stack(ax, gene, sample, shown, stype, tx, set(out_tx))
    names = [n for n in colors if n != OTHER_ISOFORMS]
    structures = load_tx_structures(tx["args"].tx_gtf, tx.get("regions", {}).get(gene), names)
    if names and structures:
        rows = len(colors)
        fig, top, h = pager.next_row(min(1.0, (0.3 + 0.15 * rows) / avail))
        a = tx["args"]
        note = f"none reached {a.tx_min_usage:.0%} usage with {a.tx_min_reads:g} reads in any sample"
        draw_isoform_structures(fig, bar_left, bar_w, top, h, list(colors), colors, structures, set(out_tx), note)
    else:
        add_isoform_legend(ax, colors, set(out_tx))
    own = [(t, m_, c) for t, m_, c in ctx["tgroups"] if t == stype]
    usage = tx["usage"][gene]
    for i in range(0, len(iso), DIST_PANELS_PER_ROW):
        chunk = iso[i:i + DIST_PANELS_PER_ROW]
        fig, top, h = pager.next_row("third")
        axes = box_axes(fig, top, h, panel_columns(len(chunk), DIST_PANELS_PER_ROW))
        for bax, rec in zip(axes, chunk):
            draw_isoform_box(bax, usage.loc[str(rec["transcript_id"])], rec, sample, shown, own)
        fit_axes_to_row(fig, axes, top, h)


def draw_other_hit_pages(pdf, args, sample, shown, ctx, genes, hit_info, expr, tx):
    global REVIEW_CTX
    if not genes:
        return
    REVIEW_CTX = ("qc", OTHER_HITS)
    pager = Pager(pdf, None)
    for gene in genes:
        sub = "isoform outlier"
        if gene in hit_info:
            rank, tier = hit_info[gene]
            sub += f"; hit: rank {rank}, tier {tier} (below the report cut-off)"

        def header(fig, part, gene=gene, sub=sub):
            nav_links(fig)
            fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H,
                     f"{shown}  |  {OTHER_HITS}  |  {gene}" + (f"  (page {part})" if part > 1 else ""),
                     fontsize=11 * FS, fontweight="bold", va="top", fontfamily=WIDE_FONT)
            fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H, sub, fontsize=6.5 * FS, va="top", color="#333",
                     fontfamily=WIDE_FONT)

        pager.new_gene(header)
        fig, top, h = pager.next_row("third")
        pdf.ohit_pages.setdefault(gene, pdf.n)
        gene_box_row(fig, top, h, gene, sample, shown, ctx, expr, args)
        draw_isoform_rows(pager, gene, sample, shown, ctx, tx, "third")
        print(f"  {gene}: isoform outlier page")
    pager.flush()
    _flush_pdf_images(pdf)


def box_axes(fig, top, h, cols, label_h=0.0):
    return [_ax(fig, left + 0.4, top - h + BOX_BOTTOM, w - 0.45, h - BOX_BOTTOM - BOX_TOP - label_h)
            for left, w in cols]


def gene_box_row(fig, top, h, gene, sample, shown, ctx, expr, args):
    axes = box_axes(fig, top, h, panel_columns(4, 4))
    draw_gene_qc(axes[:3], gene, sample, shown, ctx, args)
    draw_expression_box(axes[3], gene, sample, shown, ctx["tgroups"], expr)
    axes[3].title.set_fontfamily(WIDE_FONT)
    fit_axes_to_row(fig, axes, top, h)
    qc_frame(fig, axes[:3])
    return axes


def draw_gene_qc(axes, gene, sample, shown, ctx, args, wrap=False):
    tgroups = ctx["tgroups"]
    n_bases = ctx.get("canonical_bases", {}).get(str(gene))
    n_jxns = ctx.get("canonical_jxns", {}).get(str(gene))
    for ax, (mat, title, ylab) in zip(axes, (
            (ctx["coverage"], "Coverage", _out_of(f"Fraction of canonical bases with >={args.coverage_depth}x coverage", n_bases, "bases", wrap)),
            (ctx["testability"], "Junction testability",
             _out_of("Fraction of canonical junctions tested", n_jxns, "junctions", wrap)),
            (ctx["flr"], "Full-length ratio",
             textwrap.fill(FLR_LABEL.replace("\n", " "), 24) if wrap else FLR_LABEL))):
        gene_metric_box(ax, mat, gene, sample, shown, tgroups, title, ylab)
        ax.title.set_fontfamily(WIDE_FONT)
        fit_ylabel(ax, 4.5 * FS)


def fit_ylabel(ax, size, min_size=3.0):
    text = ax.get_ylabel()
    if not text:
        return
    avail = ax.get_position().height * PAGE_H * 0.98
    prop = FontProperties(size=size)
    longest = max(TextPath((0, 0), line, size=size, prop=prop).get_extents().width / 72.0
                  for line in text.split("\n"))
    if longest > avail:
        size = max(min_size, size * avail / longest)
    ax.yaxis.label.set_fontsize(size)


FLR_LABEL = "Average proportion of canonical transcript\ncovered by a single read"


def _out_of(label, n, unit, wrap=False):
    if wrap:
        label = textwrap.fill(label, 24)
    return label + (f"\n(out of {n:,} {unit})" if n is not None else "")


def fit_axes_to_row(fig, axes, top, h, pad=0.12):
    renderer = fig.canvas.get_renderer()
    for _ in range(4):
        bb = matplotlib.transforms.Bbox.union([ax.get_tightbbox(renderer) for ax in axes])
        lo = (top - h + pad) - bb.y0 / fig.dpi
        hi = bb.y1 / fig.dpi - (top - pad)
        lo, hi = max(lo, 0.0), max(hi, 0.0)
        if lo < 0.005 and hi < 0.005:
            return
        for ax in axes:
            x, y, w, ht = ax.get_position().bounds
            ax.set_position([x, y + lo / PAGE_H, w, max(ht - (lo + hi) / PAGE_H, 0.05)])


def qc_frame(fig, axes, pad=0.07):
    renderer = fig.canvas.get_renderer()
    bb = matplotlib.transforms.Bbox.union([ax.get_tightbbox(renderer) for ax in axes])
    x0, y0 = bb.x0 / fig.dpi - pad, bb.y0 / fig.dpi - pad
    x1, y1 = bb.x1 / fig.dpi + pad, bb.y1 / fig.dpi + pad
    fig.add_artist(FancyBboxPatch((x0 / PAGE_W, y0 / PAGE_H), (x1 - x0) / PAGE_W, (y1 - y0) / PAGE_H,
                                  boxstyle="round,pad=0,rounding_size=0.006", transform=fig.transFigure,
                                  facecolor="none", edgecolor="#8a8a8a", linewidth=0.9, zorder=0))


def _int_lookup(keys, values):
    out = {}
    for k, v in zip(keys, pd.to_numeric(values, errors="coerce")):
        if pd.notna(v) and str(k) not in out:
            out[str(k)] = int(v)
    return out


def _natural_key(gene):
    parts = re.split(r"(\d+)", str(gene).upper())
    return [(0, int(p_), "") if p_.isdigit() else (1, 0, p_) for p_ in parts if p_], str(gene)


def _letter(gene):
    c = str(gene)[:1].upper()
    return c if c.isalpha() else "#"


def draw_other_gene_pages(pdf, args, sample, shown, ctx, genes, hit_info, expr):
    global REVIEW_CTX
    if not genes:
        return
    REVIEW_CTX = ("qc", OTHER_GENES)

    def header(fig, part):
        nav_links(fig)
        fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H, f"{shown}  |  {OTHER_GENES}" + (f"  (page {part})" if part > 1 else ""),
                 fontsize=11 * FS, fontweight="bold", va="top", fontfamily=WIDE_FONT)

    pager = Pager(pdf, header)
    for gene in genes:
        fig, top, h = pager.next_row("third")
        pdf.letter_pages.setdefault(_letter(gene), pdf.n)
        label = str(gene)
        if gene in hit_info:
            rank, tier = hit_info[gene]
            label += f"   (hit: rank {rank}, tier {tier}; below the report cut-off)"
        fig.text(MARGIN_L / PAGE_W, (top - 0.17) / PAGE_H, label, fontsize=7 * FS, fontweight="bold",
                 va="bottom", ha="left")
        gene_box_row(fig, top - 0.2, h - 0.2, gene, sample, shown, ctx, expr, args)
    pager.flush()
    _flush_pdf_images(pdf)


def main():
    args = parse_args()
    sample = args.sample
    shown_name = sample
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
    all_rows = pd.read_csv(args.hits_tsv, sep="\t", dtype=str, keep_default_na=False) \
        if os.path.isfile(args.hits_tsv) else None
    if all_rows is None:
        hits = pd.DataFrame()
    elif "hit" in all_rows.columns:
        hits = all_rows[all_rows["hit"].astype(str).str.strip().str.upper() == "TRUE"].reset_index(drop=True)
    else:
        hits = all_rows
    jxn_outlier_samples = defaultdict(set)
    for _, hrow in hits.iterrows():
        if junction_calls(hrow):
            jxn_outlier_samples[str(hrow["gene"])].add(str(hrow["sample"]))
    mine = hits[hits["sample"].astype(str) == str(sample)].copy() if not hits.empty else hits
    if not mine.empty:
        mine["_rank"] = pd.to_numeric(mine["ranking"], errors="coerce")
        mine = mine[pd.to_numeric(mine["tier"], errors="coerce") <= args.max_tier]
        mine = mine.sort_values("_rank")
        if args.top_n > 0:
            mine = mine.head(args.top_n)
    genes = list(mine["gene"]) if not mine.empty else []
    hit_info = {}
    if not hits.empty:
        for _, hrow in hits[hits["sample"].astype(str) == str(sample)].iterrows():
            hit_info.setdefault(hrow["gene"], (hrow.get("ranking", "."), hrow.get("tier", ".")))
    print(f"{sample}: {len(genes)} top hit gene(s): {', '.join(genes) if genes else 'none'}")
    tx_data = load_tx(args, str(sample))
    tx_data["regions"] = load_bed(args.bed)
    tx_genes = tx_outlier_genes(tx_data)
    print(f"{sample}: {len(tx_genes)} gene(s) with an isoform outlier"
          + (f" (isoform results not available: {tx_data['unavailable']})" if tx_data["unavailable"] else ""))

    regions = load_bed(args.bed)
    tx_by_gene = load_canonical_transcripts(args.gtf, genes) if genes else {}
    hap_bams = load_gene_bams(args.gene_bam_mapping)
    sample_bams = dict(zip(map(str, args.samples), args.sample_bams)) if args.sample_bams else {}
    expr = {}
    for item in args.expression_matrices:
        st, _, path = item.partition(":")
        if path and os.path.isfile(path):
            expr[st] = load_flr_matrix(path)
    type_color = lambda st: type_colors.get(st, "#7f7f7f")
    fasta = pysam.FastaFile(args.genome)
    bulk = pysam.AlignmentFile(args.bam, "rb")

    global REVIEW_CTX
    review_genes = []
    with PdfPages(args.outfile) as real_pdf:
        pdf = CapturePdf(real_pdf, args.html_dpi if args.html_out else 0)
        REVIEW_CTX = ("overview", None)
        qc_ctx = draw_qc_page(pdf, args, str(sample), shown_name, sample_to_type, type_colors,
                              hits[hits["sample"].astype(str) == str(sample)] if not hits.empty else hits,
                              all_rows)
        if not genes:
            fig = plt.figure(figsize=(PAGE_W, PAGE_H))
            fig.text(0.5, 0.5, f"{shown_name}: no candidate gene hits.", ha="center", va="center", fontsize=14 * FS)
            pdf.savefig(fig)
            plt.close(fig)
            full = set()
        else:
            full = _draw_gene_pages(pdf, args, sample, shown_name, mine, regions, tx_by_gene, hap_bams,
                                    sample_bams, expr, jxn_outlier_samples, type_color, fasta,
                                    bulk, sample_to_type, type_colors, manifests, qc_ctx, review_genes, tx_data)
        other_hits = sorted((g for g in regions if g not in full and g in tx_genes), key=_natural_key)
        draw_other_hit_pages(pdf, args, str(sample), shown_name, qc_ctx, other_hits, hit_info, expr, tx_data)
        others = sorted((g for g in regions if g not in full and g not in tx_genes), key=_natural_key)
        draw_other_gene_pages(pdf, args, str(sample), shown_name, qc_ctx, others, hit_info, expr)
        pages = pdf.pages
        contents_rows = draw_contents(pdf, shown_name, review_genes)
    bulk.close()
    fasta.close()
    add_pdf_navigation(args.outfile, pdf, contents_rows)
    print(f"Saved: {args.outfile}")
    if args.html_out:
        import _10Q_make_review_html as review_html
        previous = review_html.load_previous_reviews(args.previous_reviews, shown_name)
        if previous:
            print(f"  {len(previous)} earlier review(s) loaded from {args.previous_reviews}")
        review_html.write_review_page(args.html_out, shown_name, qc_ctx["stype"], review_genes, pages, previous)
        print(f"Saved: {args.html_out}")


def _draw_gene_pages(pdf, args, sample, shown_name, mine, regions, tx_by_gene, hap_bams, sample_bams, expr,
                     jxn_outlier_samples, type_color, fasta, bulk, sample_to_type, type_colors,
                     manifests, qc_ctx, review_genes, tx_data):
    global REVIEW_CTX
    gene_order = [g for g in mine["gene"] if g in regions]
    if True:

        pager = Pager(pdf, None)
        for _, row in mine.iterrows():
            gene = row["gene"]
            _RESULT_CACHE.clear()
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

            gene_junctions = sorted({c["junction"] for c in calls}, key=lambda j: (_jxn_coords(j)[1:], j))
            title, sub = pdf_gene_header(shown_name, gene, row, len(variants), len(events))
            amelie_text, amelie_pmids = amelie_summary(row)

            def header(fig, part, title=title, sub=sub, prev_g=prev_g, next_g=next_g,
                       amelie_text=amelie_text, amelie_pmids=amelie_pmids):
                nav_links(fig, prev_g, next_g)
                title_artist = fig.text(MARGIN_L / PAGE_W, 1 - 0.22 / PAGE_H,
                                        title + (f"  (page {part})" if part > 1 else ""),
                                        fontsize=11 * FS, fontweight="bold", va="top", fontfamily=WIDE_FONT)
                if amelie_text:
                    right_header_text(fig, title_artist,
                                      amelie_text + (" · PMID " if amelie_pmids else ""), amelie_pmids,
                                      fontsize=6.5 * FS, fontfamily=WIDE_FONT,
                                      color="#777" if amelie_text == "AMELIE unavailable" else "#333")
                if sub:
                    fig.text(MARGIN_L / PAGE_W, 1 - 0.48 / PAGE_H, sub, fontsize=6.5 * FS, va="top", color="#333", fontfamily=WIDE_FONT)

            REVIEW_CTX = ("gene", gene)
            review_genes.append(dict(
                gene=gene, rank=str(row.get("ranking", ".")), tier=str(row.get("tier", ".")),
                hits=", ".join(f"{c}={row.get(c)}" for c in cats), cats=cats,
                inheritance=str(row.get("inheritance_patterns", "")) if _nonempty(row.get("inheritance_patterns")) else "",
                haploinsufficient=str(row.get("haploinsufficient", "")).strip().upper() == "TRUE",
                ase=str(row.get("ASE", "")).strip().upper() == "TRUE",
                variants=[dict(id=v["id"], clnsig=v["clnsig"]) for v in variants],
                junctions=gene_junctions))
            pager.new_gene(header)
            print(f"  {gene}: {len(variants)} variant(s), {len(events)} event(s), {len(pairs)} distribution(s)")
            full_left, full_w = MARGIN_L, PAGE_W - MARGIN_L - MARGIN_R
            var_marks = [v["pos"] - 1 for v in variants if v["chrom"] == chrom]
            ev_spans = [(min(_jxn_coords(j)[1] for j in e["junctions"]) - 1,
                         max(_jxn_coords(j)[2] for j in e["junctions"])) for e in events]

            tracks = [("bulk", bulk)]
            hap_handles = []
            has_ase = str(row.get("ASE", "")).strip().upper() == "TRUE"
            has_hap_jxn = any(c["phasing"] in ("hap1", "hap2") for c in calls)
            has_trans = len(variants) >= 2 and str(row.get("variants_in_trans", "")).strip().upper() == "TRUE"
            if gene in hap_bams and (has_ase or has_hap_jxn or has_trans):
                for name, path in zip(("hap1", "hap2"), hap_bams[gene]):
                    hh = pysam.AlignmentFile(path, "rb")
                    hap_handles.append(hh)
                    tracks.append((name, hh))
            fig, top, h = pager.next_row("third")
            gene_box_row(fig, top, h, gene, str(sample), shown_name, qc_ctx, expr, args)
            for k, (name, handle) in enumerate(tracks):
                if k == 1:
                    pager.break_page()
                fig, top, h = pager.next_row("two_thirds" if k == 0 else "half")
                left, width = full_left, full_w
                alignment_panel(fig, left, width, top, h, "", handle, chrom, gstart, gend, fasta, tx, args, args.seed,
                                marks=var_marks, spans=ev_spans, track_label=name,
                                show_tx=(k == 0 or k == len(tracks) - 1),
                                keep_indels=indel_specs(variants))
            for hh in hap_handles:
                hh.close()

            if variants:
                pager.break_page()
            for i in range(0, len(variants), VARIANT_PANELS_PER_ROW):
                chunk = variants[i:i + VARIANT_PANELS_PER_ROW]
                fig, top, h = pager.next_row("half")
                for k, ((left, w), v) in enumerate(zip(panel_columns(len(chunk), VARIANT_PANELS_PER_ROW), chunk)):
                    vs = v["pos"] - 1
                    ve = vs + max(len(v["ref"]), 1)
                    title = v["id"] + (f"  {v['clnsig']}" if v["clnsig"] else "")
                    alignment_panel(fig, left, w, top, h, title, bulk, v["chrom"],
                                    vs - args.variant_flank, ve + args.variant_flank, fasta, tx, args,
                                    args.seed, label_tracks=(k == 0), ref_track=True, highlight=(vs, ve),
                                    track_label="bulk", keep_indels=indel_specs([v]))

            values = load_metric_values(manifests, gene, pairs, args.coverage_threshold) if (pairs and manifests) else {}
            own_counts = load_sample_counts(manifests, gene, pairs, str(sample)) if (pairs and manifests) else {}
            hap_sig = {(c["metric"], c["junction"]) for c in calls if c["phasing"] in ("hap1", "hap2")}
            own_hap = (load_sample_hap_values(manifests, gene, sorted(hap_sig), str(sample), args.coverage_threshold)
                       if (hap_sig and manifests) else {})
            if pairs and not manifests:
                print(f"  {gene}: no cohort manifests given; distribution panels skipped.")

            def dist_rows(chunk_pairs, title=None):
                fig, top, h = pager.next_row("third")
                t_off = 0.0
                if title:
                    fig.text(full_left / PAGE_W, (top - 0.17) / PAGE_H, title, fontsize=6.5 * FS,
                             fontweight="bold", va="bottom", ha="left")
                    t_off = 0.2
                n_sub = -(-len(chunk_pairs) // DIST_PANELS_PER_ROW)
                sub_h = (h - t_off) / n_sub
                for r in range(n_sub):
                    chunk = chunk_pairs[r * DIST_PANELS_PER_ROW:(r + 1) * DIST_PANELS_PER_ROW]
                    sub_top = top - t_off - r * sub_h
                    for ax, (metric, jxn) in zip(box_axes(fig, sub_top, sub_h, panel_columns(len(chunk), DIST_PANELS_PER_ROW)),
                                                 chunk):
                        draw_distribution(ax, metric, jxn, values.get((metric, jxn), {}), str(sample), shown_name,
                                          qc_ctx["tgroups"], own_counts.get((metric, jxn)), own_hap.get((metric, jxn)))
                return fig

            psi, jcov = load_gene_psi(manifests, gene, args.coverage_threshold) if events else (None, None)
            rep_args = (gene, str(sample), [str(x) for x in args.samples], sample_to_type, sample_bams, psi, jcov,
                        jxn_outlier_samples.get(gene, set()))
            gene_rep = []
            rep_handles = {}
            igv_w = full_w * EVENT_IGV_FRAC - PANEL_GAP
            sash_left, sash_w = full_left + full_w * EVENT_IGV_FRAC, full_w * (1 - EVENT_IGV_FRAC)
            for ev_i, e in enumerate(events, 1):
                s = min(_jxn_coords(j)[1] for j in e["junctions"]) - 1
                t = max(_jxn_coords(j)[2] for j in e["junctions"])
                flank = max(int(round((t - s) * args.junction_flank)), 5)
                win = junctions_in_window(psi.columns if psi is not None else [], max(0, s - flank), t + flank)
                rep = choose_representative(*rep_args, junctions=win) if win else None
                if rep is None:
                    if not gene_rep:
                        gene_rep.append(choose_representative(*rep_args))
                    rep = gene_rep[0]
                if rep:
                    print(f"  {gene} event {ev_i}: representative {rep['sample']} ({rep['sample_type']}, {rep['scope']}: "
                          f"mean |dPSI| from median = {rep['dist']:.3f} over {rep['n_junctions']} junction(s), "
                          f"{rep['n']} candidate(s))")
                    if rep["sample"] not in rep_handles:
                        rep_handles[rep["sample"]] = pysam.AlignmentFile(sample_bams[rep["sample"]], "rb")
                    rep_bam = rep_handles[rep["sample"]]
                else:
                    print(f"  {gene} event {ev_i}: no representative sample (no covered junction PSI in other samples)")
                    rep_bam = None
                n_group = 1 + sum(len(m_) for t_, m_, _c in qc_ctx["tgroups"] if t_ == qc_ctx["stype"])
                recur = (f"  (outlier in {e['nsamples']}/{n_group} {qc_ctx['stype']} samples)"
                         if e.get("nsamples") else "")
                title = (f"{' / '.join(e['labels'])}{recur}  —  "
                         + "; ".join(f"{j} ({junction_annotation(tx, j)})" for j in e["junctions"]))
                title = fit_text(title, full_w, 6.5 * FS, "bold")
                ev_pairs = [p_ for p_ in pairs if p_[1] in e["junctions"]] if manifests else []
                is_hap_event = any(f", hap1]" in c or f", hap2]" in c for c in e["calls"])
                ev_tracks = [("bulk", None)]
                if is_hap_event and gene in hap_bams:
                    ev_tracks += list(zip(("hap1", "hap2"), hap_bams[gene]))
                pager.break_page()
                dist_fig = dist_rows(ev_pairs, title) if ev_pairs else None
                own_color = darker(type_color(sample_to_type.get(str(sample))))
                handles = [(name, bulk if path is None else pysam.AlignmentFile(path, "rb")) for name, path in ev_tracks]
                fig, top, h = pager.next_row("two_thirds")
                alignment_panel(fig, full_left, igv_w, top, h, title if fig is not dist_fig else "", bulk, chrom,
                                s - flank, t + flank, fasta, tx, args, args.seed,
                                arcs_junctions=e["junctions"], track_label="bulk", show_tx=True, title_width=full_w)
                sash_tracks = [dict(label=f"{sample} {name}", bam=handle, color=own_color) for name, handle in handles]
                if rep_bam is not None:
                    sash_tracks.append(dict(label=f"{rep['sample']} (representative {rep['sample_type']}"
                                                  + (", gene-wide)" if rep["scope"] == "gene-wide" else ")"),
                                            bam=rep_bam, color=type_color(rep["sample_type"])))
                else:
                    sash_tracks[0]["label"] += "  (no representative sample)"
                draw_sashimi(fig, sash_left, sash_w, top - 0.19 - 0.06, top - h + 0.10, sash_tracks, tx, chrom,
                             s - flank, t + flank, args)
                if len(handles) > 1:
                    pager.break_page()
                    for k, (name, handle) in enumerate(handles[1:]):
                        fig, top, h = pager.next_row("half")
                        alignment_panel(fig, full_left, igv_w, top, h, title if k == 0 else "", handle, chrom,
                                        s - flank, t + flank, fasta, tx, args, args.seed,
                                        arcs_junctions=e["junctions"], track_label=name,
                                        show_tx=(k == len(handles) - 2), title_width=full_w)
                for name, handle in handles:
                    if handle is not bulk:
                        handle.close()
            for hh in rep_handles.values():
                hh.close()
            pager.break_page()
            draw_isoform_rows(pager, gene, str(sample), shown_name, qc_ctx, tx_data, "half")

        pager.flush()
        return set(gene_order)


if __name__ == "__main__":
    main()