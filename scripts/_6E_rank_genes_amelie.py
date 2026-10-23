import argparse
import json
import os
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid

import pandas as pd

AMELIE_URL = "https://amelie.stanford.edu/api/gene_list_api/"
MAX_GENES = 1000
HPO_RE = re.compile(r"^HP:\d{7}$")


def parse_args():
    p = argparse.ArgumentParser(
        description="Rank the panel genes for one sample by how well the published literature links them to the "
                    "sample's HPO terms, using the AMELIE gene-list API. Only HPO IDs and gene symbols are sent; "
                    "the submission is labelled with a random ID, not the sample name. Panels larger than AMELIE's "
                    "limit are sent in batches and the genes are ranked together by their AMELIE score. If AMELIE "
                    "cannot be reached, every gene is written with status 'error' and the script still exits "
                    "successfully.")
    p.add_argument("--sample", required=True)
    p.add_argument("--hpo", required=True,
        help="Semicolon-separated HPO terms: IDs (HP:0000000), term names, or a mix. Names are matched "
             "case-insensitively against HPO term names and exact synonyms in --hpo-obo.")
    p.add_argument("--hpo-obo", default=None, help="HPO ontology (hp.obo); needed to resolve term names.")
    p.add_argument("--terms-outfile", required=True,
        help="Writes each input term with its resolved HPO ID and name, and how it was matched.")
    p.add_argument("--bed", required=True, help="Panel BED; column 4 is the gene symbol.")
    p.add_argument("--outfile", required=True)
    p.add_argument("--url", default=AMELIE_URL)
    p.add_argument("--timeout", type=int, default=600, help="Seconds to wait for each AMELIE request.")
    p.add_argument("--retries", type=int, default=3)
    p.add_argument("--n-pmids", type=int, default=3, help="Supporting PMIDs to keep per gene.")
    p.add_argument("--batch-size", type=int, default=MAX_GENES,
        help=f"Genes per AMELIE request (at most {MAX_GENES}).")
    return p.parse_args()


def load_panel_genes(bed):
    genes = []
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] not in genes:
                genes.append(fields[3])
    return genes


def normalize(text):
    return re.sub(r"\s+", " ", text.strip().lower()).rstrip(".")


def load_obo(path):
    names, by_name, by_synonym, alt_ids, replaced = {}, {}, {}, {}, {}
    term = None

    def finish(t):
        if not t or "id" not in t:
            return
        tid = t["id"]
        if t.get("obsolete"):
            if t.get("replaced_by"):
                replaced[tid] = t["replaced_by"]
            return
        names[tid] = t.get("name", tid)
        by_name.setdefault(normalize(t.get("name", "")), tid)
        for syn in t.get("synonyms", []):
            by_synonym.setdefault(normalize(syn), tid)
        for alt in t.get("alt_ids", []):
            alt_ids[alt] = tid

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                finish(term)
                term = {} if line == "[Term]" else None
                continue
            if term is None or ": " not in line:
                continue
            key, value = line.split(": ", 1)
            if key == "id":
                term["id"] = value.strip()
            elif key == "name":
                term["name"] = value.strip()
            elif key == "synonym":
                m = re.match(r'"(.*)"\s+EXACT', value)
                if m:
                    term.setdefault("synonyms", []).append(m.group(1))
            elif key == "alt_id":
                term.setdefault("alt_ids", []).append(value.strip())
            elif key == "is_obsolete" and value.strip() == "true":
                term["obsolete"] = True
            elif key == "replaced_by":
                term["replaced_by"] = value.strip()
    finish(term)
    return dict(names=names, by_name=by_name, by_synonym=by_synonym, alt_ids=alt_ids, replaced=replaced)


def split_terms(text):
    pieces = []
    for piece in text.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        parts = [p.strip() for p in piece.split(",") if p.strip()]
        if len(parts) > 1 and all(HPO_RE.match(p) for p in parts):
            pieces.extend(parts)
        else:
            pieces.append(piece)
    return pieces


def resolve_terms(text, obo):
    rows = []
    for raw in split_terms(text):
        tid, how = None, None
        if HPO_RE.match(raw):
            if obo is None:
                tid, how = raw, "id"
            elif raw in obo["names"]:
                tid, how = raw, "id"
            elif raw in obo["alt_ids"]:
                tid, how = obo["alt_ids"][raw], "alt_id"
            elif raw in obo["replaced"]:
                tid, how = obo["replaced"][raw], "replaced_by"
        elif obo is not None:
            key = normalize(raw)
            if key in obo["by_name"]:
                tid, how = obo["by_name"][key], "name"
            elif key in obo["by_synonym"]:
                tid, how = obo["by_synonym"][key], "synonym"
        name = obo["names"].get(tid, ".") if (obo is not None and tid) else "."
        rows.append(dict(input_term=raw, hpo_id=tid or ".", hpo_name=name, matched_by=how or "unmatched"))
    return rows


UNVERIFIED = {"used": False}


def unverified_context():
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def query_amelie(url, hpo, genes, timeout, retries):
    body = urllib.parse.urlencode({
        "patientName": str(uuid.uuid4()),
        "phenotypes": ",".join(hpo),
        "genes": ",".join(genes),
    }).encode()
    last_error = None
    attempt = 0
    while attempt < retries:
        attempt += 1
        try:
            req = urllib.request.Request(url, data=body, method="POST")
            ctx = unverified_context() if UNVERIFIED["used"] else None
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
                return json.loads(resp.read().decode()), None
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
            reason = getattr(e, "reason", e)
            if not UNVERIFIED["used"] and isinstance(reason, ssl.SSLCertVerificationError):
                UNVERIFIED["used"] = True
                attempt -= 1
                print(f"  WARNING: AMELIE's SSL certificate could not be verified ({reason.verify_message}); "
                      f"retrying without certificate verification.")
                continue
            last_error = f"{type(e).__name__}: {e}"
            print(f"  attempt {attempt}/{retries} failed: {last_error}")
            if attempt < retries:
                time.sleep(min(60, 10 * attempt))
    return None, last_error


def batches(genes, size):
    n = -(-len(genes) // size)
    step = -(-len(genes) // n)
    return [genes[i:i + step] for i in range(0, len(genes), step)]


def query_in_batches(url, hpo, genes, size, timeout, retries):
    chunks = batches(genes, size)
    data = []
    for i, chunk in enumerate(chunks, start=1):
        print(f"Batch {i}/{len(chunks)}: {len(chunk)} gene(s)")
        t0 = time.time()
        part, error = query_amelie(url, hpo, chunk, timeout, retries)
        if part is None:
            return None, f"batch {i}/{len(chunks)}: {error}"
        print(f"  answered in {time.time() - t0:.0f}s; {len(part)} gene(s) ranked")
        data.extend(part)
    return data, None


def to_table(sample, genes, data, n_pmids):
    ranked = []
    for entry in data:
        if not isinstance(entry, (list, tuple)) or len(entry) != 2:
            continue
        gene, papers = entry
        pairs = []
        for paper in papers or []:
            try:
                pairs.append((str(paper[0]), float(paper[1])))
            except (TypeError, ValueError, IndexError):
                continue
        pairs.sort(key=lambda x: -x[1])
        ranked.append((str(gene), pairs))
    ranked.sort(key=lambda x: (-(x[1][0][1] if x[1] else float("-inf")), x[0]))
    by_gene = {}
    for rank, (gene, pairs) in enumerate(ranked, start=1):
        by_gene[gene] = (rank, round(pairs[0][1], 2) if pairs else ".",
                         ",".join(p for p, _ in pairs[:n_pmids]) or ".")
    rows = []
    for gene in genes:
        rank, score, pmids = by_gene.get(gene, (".", ".", "."))
        rows.append(dict(sample=sample, gene=gene, AMELIE_rank=rank, AMELIE_score=score, AMELIE_pmids=pmids,
                         AMELIE_n_genes=len(genes), AMELIE_status="ok"))
    return pd.DataFrame(rows)


def main():
    print(f"\n\n\n******************************************************************************************")
    print(f"Ranking panel genes with AMELIE...")
    print(f"******************************************************************************************\n")
    args = parse_args()
    genes = load_panel_genes(args.bed)
    obo = None
    if args.hpo_obo and os.path.isfile(args.hpo_obo):
        obo = load_obo(args.hpo_obo)
        print(f"Loaded {len(obo['names'])} HPO terms from {args.hpo_obo}")
    else:
        print(f"WARNING: HPO ontology not found ({args.hpo_obo}); only HPO IDs can be used and term names "
              f"will not be shown.")
    terms = resolve_terms(args.hpo, obo)
    for t in terms:
        if t["matched_by"] == "unmatched":
            print(f"  WARNING: '{t['input_term']}' does not match any HPO term; it will not be sent to AMELIE.")
        else:
            print(f"  '{t['input_term']}' -> {t['hpo_id']} {t['hpo_name']} ({t['matched_by']})")
    hpo = list(dict.fromkeys(t["hpo_id"] for t in terms if t["hpo_id"] != "."))
    print(f"{args.sample}: {len(hpo)} HPO term(s) to query: {', '.join(hpo) or 'none'}")
    print(f"{len(genes)} panel gene(s) from {args.bed}")

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.terms_outfile)), exist_ok=True)
    pd.DataFrame(terms, columns=["input_term", "hpo_id", "hpo_name", "matched_by"]).to_csv(
        args.terms_outfile, sep="\t", index=False)
    error = None
    if not hpo:
        error = "no valid HPO terms"
    elif not genes:
        error = "no panel genes"
    data = None
    if error is None:
        size = max(1, min(args.batch_size, MAX_GENES))
        print(f"Querying {args.url} ...")
        data, error = query_in_batches(args.url, hpo, genes, size, args.timeout, args.retries)
        if data is not None:
            print(f"AMELIE ranked {len(data)} gene(s) in total.")

    if data is None:
        print(f"WARNING: AMELIE ranking unavailable for {args.sample} ({error}). "
              f"Writing status 'error'; rerun this rule later with --forcerun _6E_rank_genes_amelie.")
        df = pd.DataFrame([dict(sample=args.sample, gene=g, AMELIE_rank=".", AMELIE_score=".", AMELIE_pmids=".",
                                AMELIE_n_genes=len(genes), AMELIE_status="error") for g in genes])
    else:
        df = to_table(args.sample, genes, data, args.n_pmids)
        top = df[df["AMELIE_rank"] != "."].sort_values("AMELIE_rank").head(10)
        if not top.empty:
            print("\nTop-ranked genes:")
            print(top[["gene", "AMELIE_rank", "AMELIE_score", "AMELIE_pmids"]].to_string(index=False))

    df.to_csv(args.outfile, sep="\t", index=False)
    print(f"\nWrote {args.outfile}")


if __name__ == "__main__":
    main()
