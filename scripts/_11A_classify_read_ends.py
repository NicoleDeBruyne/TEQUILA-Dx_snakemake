import argparse
import gzip
import importlib.util
import os
import re
from bisect import bisect_right
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pysam

CATEGORIES = ["oligo-dT artifact", "TSO HH artifact", "TSO LH artifact", "TSO NTA artifact",
              "oligo-dT", "alu-polyA", "TSO", "unknown"]
CLUSTER_COLS = ["sample", "gene", "chrom", "strand", "side", "end", "pos", "n_clipped", "depth", "clip_fraction",
                "median_clip_len", "max_clip_len", "consensus", "consensus_identity", "null_identity",
                "hap1_clipped", "hap2_clipped", "unphased_clipped", "location"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Label both ends of every read on a panel gene as a real end (TSO / oligo-dT), a priming "
                    "artifact, or unknown, and gather unknown soft clips of at least --min-clip bases into clusters "
                    "of reads clipped at exactly the same base on the same side. Depth is the number of aligned "
                    "bases at that base among reads assigned to the gene. The clipped sequences (the "
                    "--consensus-len bases next to the breakpoint) are aligned with partial-order alignment "
                    "(abPOA); consensus_identity is the median identity of each clip to the consensus, and "
                    "null_identity is the same measure for as many clips taken one each from different, randomly chosen "
                    "breakpoints elsewhere in the sample, as a reference for unrelated sequence.")
    p.add_argument("--sample", required=True)
    p.add_argument("--bam", required=True, help="Untrimmed long-read BAM (TSO and poly(A) still soft-clipped).")
    p.add_argument("--bed", required=True, help="Panel BED; column 4 is the gene symbol.")
    p.add_argument("--gtf", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--alu", default=None, help="Alu annotation (BED or GTF), optional; enables alu-polyA.")
    p.add_argument("--gene-bam-mapping", default=None,
                   help="Rule _3A's gene_bam_mapping_file.tsv; used to count clipped reads per haplotype.")
    p.add_argument("--genes", nargs="*", default=[],
                   help="Only analyse these panel genes (gene symbols). Default: every gene in --bed.")
    p.add_argument("--min-clip", type=int, default=50)
    p.add_argument("--min-cluster-reads", type=int, default=2,
                   help="Clusters with fewer clipped reads are not written.")
    p.add_argument("--consensus-len", type=int, default=50)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--nta-nuc", default="G")
    p.add_argument("--fwd-tso-seq", default="AAGCAGTGGTATCAACGCAGAGTACAT")
    p.add_argument("--threads", type=int, default=1, help="Genes are processed in parallel on this many cores.")
    p.add_argument("--clusters-outfile", required=True)
    p.add_argument("--counts-outfile", required=True)
    return p.parse_args()


def _complement(s):
    return s.translate(str.maketrans("ATCGatcg", "TAGCtagc"))


def _variants_1mm(m):
    b = "ATCG"
    n = len(m)
    return ([m[:i] + x + m[i + 1:] for i in range(n) for x in b if x != m[i]] +
            [m[:i] + m[i + 1:] for i in range(1, n - 1)] +
            [m[:i] + x + m[i:] for i in range(1, n) for x in b] + [m])


def build_tso_patterns(nta, tso):
    cnta = _complement(nta)
    cmp_ = _complement(tso)
    rev = tso[::-1]
    rcmp = cmp_[::-1]
    P = {"nta_dimer": nta * 2, "cnta_dimer": cnta * 2,
         "exact_fwd": tso + nta * 3, "exact_cmp": cmp_ + cnta * 3,
         "exact_rev": nta * 3 + rev, "exact_rcmp": cnta * 3 + rcmp}
    lead = ({s + nta * k for s in _variants_1mm(tso[-6:]) for k in range(2, 6)} |
            {s + cnta * k for s in _variants_1mm(cmp_[-6:]) for k in range(2, 6)})
    trail = ({nta * k + s for s in _variants_1mm(rev[:6]) for k in range(2, 6)} |
             {cnta * k + s for s in _variants_1mm(rcmp[:6]) for k in range(2, 6)})
    P["lead"], P["trail"] = defaultdict(set), defaultdict(set)
    for p in lead:
        P["lead"][len(p)].add(p)
    for p in trail:
        P["trail"][len(p)].add(p)
    return P


def _oligo_dt(sc, g, leading):
    frag = sc[-10:] if leading else sc[:10]
    if len(frag) < 10 or not (frag.count("A") >= 8 or frag.count("T") >= 8):
        return None
    f5 = g[:5] if leading else g[-5:]
    f10 = g[:10] if leading else g[-10:]
    if f5 in ("AAAAA", "TTTTT") or f10.count("A") >= 8 or f10.count("T") >= 8:
        return "oligo-dT artifact"
    return "oligo-dT"


def _edit(g, p, m):
    if g == p:
        return 0
    if abs(len(g) - len(p)) > m:
        return -1
    prev = list(range(len(p) + 1))
    for i, cg in enumerate(g, 1):
        cur = [i] + [0] * len(p)
        for j, cp in enumerate(p, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (cg != cp))
        if min(cur) > m:
            return -1
        prev = cur
    return prev[-1] if prev[-1] <= m else -1


def _tso(sc, g, leading, P):
    if not sc:
        return None
    pats = P["lead"] if leading else P["trail"]
    if any(len(sc) >= n and (sc[-n:] if leading else sc[:n]) in s for n, s in pats.items()):
        return "TSO"
    if len(sc) < 6:
        return None
    if leading:
        anchor, keys, nd = sc[-6:], ("exact_fwd", "exact_cmp"), P["nta_dimer"]
    else:
        anchor, keys, nd = sc[:6], ("exact_rev", "exact_rcmp"), P["cnta_dimer"]
    for k in keys:
        pat = P[k]
        idx = pat.find(anchor)
        if idx == -1:
            continue
        n = len(pat) - idx - 6 if leading else idx
        if n < 2:
            continue
        gp = g[:n] if leading else g[-n:]
        pp = pat[-n:] if leading else pat[:n]
        hh, lh = max(0, (n - 3) // 3), max(0, (n - 2) // 2)
        if n == 2 and gp == nd:
            return "TSO NTA artifact"
        ed = _edit(gp, pp, lh)
        if 0 <= ed <= hh:
            return "TSO HH artifact"
        if ed != -1:
            return "TSO LH artifact"
    return None


def _in_alu(chrom, pos0, alu):
    if alu is None or chrom not in alu:
        return False
    starts, ends = alu[chrom]
    i = bisect_right(starts, pos0) - 1
    return i >= 0 and pos0 < ends[i]


def classify_end(sc, g, leading, P, chrom, pos0, alu):
    lab = _oligo_dt(sc, g, leading)
    if lab:
        if lab == "oligo-dT artifact" and _in_alu(chrom, pos0, alu):
            return "alu-polyA"
        return lab
    return _tso(sc, g, leading, P) or "unknown"


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_alu(path):
    ivs = defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                s, e = (int(f[3]) - 1, int(f[4])) if len(f) >= 9 and f[3].isdigit() and f[4].isdigit() \
                    else (int(f[1]), int(f[2]))
            except (ValueError, IndexError):
                continue
            ivs[f[0]].append((s, e))
    out = {}
    for c, lst in ivs.items():
        lst.sort()
        out[c] = ([s for s, _ in lst], [e for _, e in lst])
    return out


def load_panel_genes(bed):
    genes = []
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            g = line.rstrip("\n").split("\t")[3]
            if g not in genes:
                genes.append(g)
    return genes


_ATTR = {k: re.compile(k + r'\s+"([^"]+)"') for k in ("gene_name", "gene_id", "transcript_id")}


def load_gene_models(gtf, wanted):
    wanted = set(wanted)
    tx_exons = defaultdict(lambda: defaultdict(list))
    info = {}
    with _open(gtf) as fh:
        for line in fh:
            if "\texon\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = _ATTR["gene_name"].search(f[8]) or _ATTR["gene_id"].search(f[8])
            if not m or m.group(1) not in wanted:
                continue
            t = _ATTR["transcript_id"].search(f[8])
            if not t:
                continue
            g = m.group(1)
            info.setdefault(g, (f[0], f[6]))
            tx_exons[g][t.group(1)].append((int(f[3]) - 1, int(f[4])))
    genes = {}
    for g, txs in tx_exons.items():
        chrom, strand = info[g]
        introns = set()
        for ex in txs.values():
            ex = sorted(ex)
            for (_s1, e1), (s2, _e2) in zip(ex, ex[1:]):
                introns.add((e1, s2))
        merged = []
        for s, e in sorted(x for ex in txs.values() for x in ex):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        bounds = sorted({x for ex in txs.values() for iv in ex for x in iv})
        genes[g] = dict(chrom=chrom, strand=strand, start=merged[0][0], end=merged[-1][1],
                        exons=[tuple(x) for x in merged], introns=frozenset(introns), bounds=bounds)
    return genes


def read_blocks_introns(r):
    rp = r.reference_start
    blocks, introns = [], set()
    for op, ln in r.cigartuples:
        if op in (0, 7, 8):
            blocks.append((rp, rp + ln))
            rp += ln
        elif op == 2:
            rp += ln
        elif op == 3:
            introns.add((rp, rp + ln))
            rp += ln
    return blocks, introns


def exon_overlap(blocks, exons):
    tot = 0
    for s, e in blocks:
        for a, b in exons:
            if a >= e:
                break
            if b > s:
                tot += min(e, b) - max(s, a)
    return tot


def assign_gene(r, cands, genes):
    blocks, introns = read_blocks_introns(r)
    scored = []
    for g in cands:
        ov = exon_overlap(blocks, genes[g]["exons"])
        if ov > 0:
            scored.append(((len(introns & genes[g]["introns"]), ov), g))
    if not scored:
        return None
    best = max(sc for sc, _ in scored)
    top = [g for sc, g in scored if sc == best]
    return top[0] if len(top) == 1 else None


def location(pos0, gene):
    if any(a <= pos0 < b for a, b in gene["exons"]):
        return "exonic"
    i = bisect_right(gene["bounds"], pos0)
    near = [gene["bounds"][j] for j in (i - 1, i) if 0 <= j < len(gene["bounds"])]
    if any(abs(b - pos0) <= 5 for b in near):
        return "exon boundary"
    if gene["start"] <= pos0 < gene["end"]:
        return "intronic"
    return "outside gene"


def levenshtein(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return prev[-1]


def poa_consensus(seqs):
    import pyabpoa
    res = pyabpoa.msa_aligner().msa(list(seqs), out_cons=True, out_msa=False)
    cons = res.cons_seq[0] if res.cons_seq else ""
    ident = [1 - levenshtein(x, cons) / max(len(x), len(cons), 1) for x in seqs]
    return cons, round(float(np.median(ident)), 4)


def haplotype_reads(mapping, gene):
    if not mapping or not os.path.isfile(mapping):
        return {}
    m = pd.read_csv(mapping, sep="\t", dtype=str)
    row = m[m["gene"] == gene]
    out = {}
    if row.empty:
        return out
    for hap in ("hap1", "hap2"):
        path = row.iloc[0].get(f"{hap}_bam")
        if path and isinstance(path, str) and os.path.isfile(path):
            with pysam.AlignmentFile(path, "rb") as b:
                for r in b.fetch(until_eof=True):
                    out[r.query_name] = hap
    return out


_W = {}


def _init_worker(genes, args):
    _W.update(genes=genes, args=args, P=build_tso_patterns(args.nta_nuc, args.fwd_tso_seq),
              alu=load_alu(args.alu) if args.alu else None)


def process_gene(g):
    genes, a, P, alu = _W["genes"], _W["args"], _W["P"], _W["alu"]
    d = genes[g]
    chrom = d["chrom"]
    counts, events = Counter(), defaultdict(list)
    assigned = set()
    rows = []
    with pysam.AlignmentFile(a.bam, "rb") as bam, pysam.FastaFile(a.genome) as fa:
        if chrom not in bam.references:
            return g, counts, rows, 0, 0, []
        pad = 100_000
        r0 = max(0, d["start"] - pad)
        ref = fa.fetch(chrom, r0, d["end"] + pad).upper()
        cands_all = [x for x in genes if genes[x]["chrom"] == chrom and genes[x]["start"] < d["end"] + pad
                     and genes[x]["end"] > d["start"] - pad]

        def seq(s, e):
            if s >= r0 and e - r0 <= len(ref):
                return ref[s - r0:e - r0]
            return fa.fetch(chrom, max(0, s), e).upper()

        for r in bam.fetch(chrom, d["start"], d["end"]):
            if r.flag & 0x904 or not r.cigartuples:
                continue
            rs, re_ = r.reference_start, r.reference_end
            cands = [x for x in cands_all if genes[x]["start"] < re_ and genes[x]["end"] > rs]
            if len(cands) > 1 and assign_gene(r, cands, genes) != g:
                continue
            if len(cands) == 1 and exon_overlap(read_blocks_introns(r)[0], d["exons"]) == 0:
                continue
            q = r.query_sequence
            if q is None:
                continue
            assigned.add(r.query_name)
            sc_l, sc_r = q[:r.query_alignment_start], q[r.query_alignment_end:]
            left = classify_end(sc_l, seq(rs, rs + 20), True, P, chrom, rs, alu)
            right = classify_end(sc_r, seq(max(0, re_ - 20), re_), False, P, chrom, re_ - 1, alu)
            five, three = (right, left) if d["strand"] == "-" else (left, right)
            counts[(g, "5'", five)] += 1
            counts[(g, "3'", three)] += 1
            for side, lab, sc, pos in (("left", left, sc_l, rs), ("right", right, sc_r, re_ - 1)):
                if lab == "unknown" and len(sc) >= a.min_clip:
                    end = ("5'" if side == "left" else "3'") if d["strand"] == "+" else \
                        ("3'" if side == "left" else "5'")
                    clip = sc[::-1][:a.consensus_len] if side == "left" else sc[:a.consensus_len]
                    events[(side, pos)].append(dict(clip_len=len(sc), read=r.query_name, end=end, seq=clip))
        n_clips = sum(len(v) for v in events.values())
        clusters = [(k, v) for k, v in events.items() if len(v) >= a.min_cluster_reads]
        haps = haplotype_reads(a.gene_bam_mapping, g) if clusters else {}
        for (side, base), group in sorted(clusters):
            depth = int(np.sum(bam.count_coverage(chrom, base, base + 1, quality_threshold=0,
                                                  read_callback=lambda x: x.query_name in assigned
                                                  and not (x.flag & 0x904))))
            depth = max(depth, len(group))
            seqs = [e["seq"] for e in group]
            cons, ident = poa_consensus(seqs)
            hc = Counter(haps.get(e["read"], "unphased") for e in group)
            lens = [e["clip_len"] for e in group]
            rows.append(dict(sample=a.sample, gene=g, chrom=chrom, strand=d["strand"], side=side,
                             end=Counter(e["end"] for e in group).most_common(1)[0][0],
                             pos=base + 1, n_clipped=len(group), depth=depth,
                             clip_fraction=round(len(group) / depth, 4),
                             median_clip_len=int(np.median(lens)), max_clip_len=int(max(lens)),
                             consensus=cons[::-1] if side == "left" else cons, consensus_identity=ident,
                             null_identity="", _site=(g, side, base), hap1_clipped=hc.get("hap1", 0), hap2_clipped=hc.get("hap2", 0),
                             unphased_clipped=hc.get("unphased", 0), location=location(base, d)))
    pool = [((g, side, b), [e["seq"] for e in v]) for (side, b), v in events.items()]
    return g, counts, rows, len(assigned), n_clips, pool


def main():
    a = parse_args()
    if importlib.util.find_spec("pyabpoa") is None:
        raise SystemExit("ERROR: pyabpoa is needed for the clip consensus; install it in the pipeline's conda "
                         "environment (pip install pyabpoa, or conda install -c bioconda pyabpoa)")
    panel = load_panel_genes(a.bed)
    if a.genes:
        absent = [g for g in a.genes if g not in panel]
        if absent:
            print(f"  Not on this panel: {', '.join(absent)}")
        panel = [g for g in panel if g in set(a.genes)]
    genes = load_gene_models(a.gtf, panel)
    missing = [g for g in panel if g not in genes]
    print(f"{a.sample}: {len(genes)} of {len(panel)} panel genes found in {a.gtf}")
    if missing:
        print(f"  Not in the GTF: {', '.join(missing[:20])}{' ...' if len(missing) > 20 else ''}")
    order = sorted(genes, key=lambda g: (genes[g]["chrom"], genes[g]["start"]))
    counts, rows = Counter(), []
    n_reads = n_clips = 0
    if a.threads > 1:
        with ProcessPoolExecutor(max_workers=a.threads, initializer=_init_worker, initargs=(genes, a)) as ex:
            results = list(ex.map(process_gene, order, chunksize=4))
    else:
        _init_worker(genes, a)
        results = [process_gene(g) for g in order]
    pool = []
    for _g, c, r, n, k, pl in results:
        counts.update(c)
        rows.extend(r)
        n_reads += n
        n_clips += k
        pool.extend(pl)
    rng = np.random.default_rng(a.seed)
    for row in rows:
        site = row.pop("_site")
        idx = [i for i, (key, _s) in enumerate(pool) if key != site]
        if len(idx) >= row["n_clipped"]:
            pick = rng.choice(idx, row["n_clipped"], replace=False)
            row["null_identity"] = poa_consensus([pool[i][1][rng.integers(len(pool[i][1]))] for i in pick])[1]
    print(f"  {n_reads:,} reads assigned to panel genes; {n_clips:,} unknown soft clips >= {a.min_clip} bp")
    out = pd.DataFrame(rows, columns=CLUSTER_COLS).sort_values(["n_clipped", "gene"], ascending=[False, True])
    os.makedirs(os.path.dirname(a.clusters_outfile) or ".", exist_ok=True)
    out.to_csv(a.clusters_outfile, sep="\t", index=False)
    print(f"  {len(out):,} cluster(s) with >= {a.min_cluster_reads} clipped reads -> {a.clusters_outfile}")

    crow = [dict(sample=a.sample, gene=g, end=e, category=c, reads=n) for (g, e, c), n in counts.items()]
    cdf = pd.DataFrame(crow, columns=["sample", "gene", "end", "category", "reads"])
    if not cdf.empty:
        cdf = (cdf.pivot_table(index=["sample", "gene", "end"], columns="category", values="reads",
                               aggfunc="sum", fill_value=0)
               .reindex(columns=CATEGORIES, fill_value=0).reset_index())
        cdf["total"] = cdf[CATEGORIES].sum(axis=1)
    else:
        cdf = pd.DataFrame(columns=["sample", "gene", "end"] + CATEGORIES + ["total"])
    cdf.to_csv(a.counts_outfile, sep="\t", index=False)
    print(f"  Read-end label counts -> {a.counts_outfile}")


if __name__ == "__main__":
    main()