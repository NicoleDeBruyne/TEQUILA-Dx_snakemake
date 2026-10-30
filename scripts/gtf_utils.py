# GTF / alignment helpers shared between scripts

import gzip
import re
from collections import defaultdict
from typing import Dict, List, Tuple

# M, D, =, X. N is handled separately since it splits a read into blocks
_REF_CONSUMING_NON_N = {0, 2, 7, 8}
# regex version of the attribute parser, compiled once per key since it gets called a lot
_ATTR_RE_CACHE = {}


def get_attr(attr_str, key):
    rx = _ATTR_RE_CACHE.get(key)
    if rx is None:
        rx = re.compile(key + r'\s+"([^"]+)"')
        _ATTR_RE_CACHE[key] = rx
    m = rx.search(attr_str)
    return m.group(1) if m else None


def strip_version(s):
    return re.sub(r'\.\d+$', '', s) if s else s


def _open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _gtf_attr(attr_str, key):
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


def load_gene_list(bed):
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


# one transcript per gene: Ensembl_canonical if tagged, else the longest one.
# exons get merged so overlaps aren't double counted
def load_canonical_transcripts(gtf_path, genes):
    genes = set(genes)

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
            else:
                t = transcripts.setdefault(transcript_id, {
                    "gene": gene_name, "chrom": chrom, "strand": strand,
                    "canonical": False, "exons": [],
                })
                t["exons"].append((int(start) - 1, int(end)))

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


# aligned reference blocks of a read, broken at N (introns)
def ref_blocks(read):
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
            if op == 3:
                pos += length
    if block_start is not None:
        blocks.append((block_start, pos))
    return blocks


# annotated introns per panel gene, built from consecutive exons of each transcript
def parse_gtf_junctions(
    gtf_path: str,
    gene_names: List[str],
) -> Dict[str, Dict]:
    gene_set = set(gene_names)

    tx_canonical: set = set()
    gene_exons: Dict[str, Dict[str, List[Tuple[str, int, int]]]] = defaultdict(lambda: defaultdict(list))

    open_fn = gzip.open if gtf_path.endswith(".gz") else open

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

            gname = get_attr(attrs, "gene_name") or get_attr(attrs, "gene_symbol")
            gid   = strip_version(get_attr(attrs, "gene_id"))
            gene_key = None
            if gname and gname in gene_set:
                gene_key = gname
            elif gid and gid in gene_set:
                gene_key = gid
            if gene_key is None:
                continue

            tx_id = get_attr(attrs, "transcript_id")
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
