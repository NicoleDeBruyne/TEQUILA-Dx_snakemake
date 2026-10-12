"""
review_html.py

Writes the offline HTML review page that goes next to each sample's PDF report: the same pages
(as embedded PNG images), one candidate gene at a time, with a review status and notes per gene.

Everything is in the one file (images as data URIs, inline CSS and JavaScript); nothing is loaded
from or sent to the network. Reviews are saved in the browser's local storage for this file (keyed
by sample and the time the page was generated, so a regenerated page starts fresh), and can be
exported / imported as TSV. "Export PDF" lays out every page with each gene's status and notes and
opens the browser's print dialog ("Save as PDF"); the pipeline's PDF is never modified.

Per gene: one of STATUSES (or none yet; clicking the selected status again clears it), a
"Reviewed externally" checkbox, and free-text notes.
"""
import base64
import datetime
import html
import json

STATUSES = [
    ("good", "Good candidate"),
    ("possible", "Possible candidate"),
    ("unlikely", "Unlikely candidate"),
    ("unclear", "Unclear - needs further review"),
]
NO_STATUS = ("none", "No status yet")       # every gene starts here; not a button


def write_review_page(path, sample, sample_type, genes, pages):
    """genes: list of dict(gene, rank, tier, hits, inheritance, haploinsufficient), in report order.
    pages: list of ((kind, name), png_bytes) with kind 'overview', 'gene' (name = gene) or 'qc'
    (name = section title), in PDF order."""
    page_list = [{"kind": k, "name": n, "src": "data:image/png;base64," + base64.b64encode(png).decode()}
                 for (k, n), png in pages]
    data = {
        "sample": str(sample),
        "sampleType": str(sample_type),
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "statuses": [{"key": k, "label": l} for k, l in STATUSES],
        "noStatus": {"key": NO_STATUS[0], "label": NO_STATUS[1]},
        "genes": genes,
        "pages": page_list,
    }
    blob = json.dumps(data).replace("</", "<\\/")
    out = (_TEMPLATE
           .replace("__TITLE__", html.escape(f"{sample} review"))
           .replace("__DATA__", blob))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(out)


_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--bg:#f7f7f5;--panel:#fff;--line:#dedcd5;--text:#222;--muted:#6b6a65;--accent:#0068a9;
 --good:#3d892e;--possible:#0068a9;--unlikely:#6e7f8f;--unclear:#c69528;--none:#9a9890}
*{box-sizing:border-box}
body{margin:0;font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
 color:var(--text);background:var(--bg)}
#app{display:grid;grid-template-columns:300px 1fr;grid-template-rows:auto 1fr;height:100vh}
header{grid-column:1/3;display:flex;align-items:center;gap:16px;padding:10px 16px;background:var(--panel);
 border-bottom:1px solid var(--line)}
header h1{font-size:17px;margin:0;font-weight:600}
header .meta{color:var(--muted)}
header .spacer{flex:1}
button,select{font:inherit;padding:5px 10px;border:1px solid var(--line);border-radius:6px;background:#fff;
 cursor:pointer;color:var(--text)}
button:hover{border-color:#aaa}
button.primary{background:var(--accent);border-color:var(--accent);color:#fff}
#progress{color:var(--muted);white-space:nowrap}
aside{overflow:auto;background:var(--panel);border-right:1px solid var(--line);padding:10px 0}
aside h2{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:14px 14px 6px}
aside .filter{margin:0 14px 6px;width:calc(100% - 28px)}
.item{display:flex;align-items:center;gap:8px;padding:6px 14px;cursor:pointer;border-left:3px solid transparent}
.item:hover{background:#f1f0ec}
.item.active{background:#e9f1f8;border-left-color:var(--accent)}
.item .rank{color:var(--muted);width:24px;text-align:right;font-variant-numeric:tabular-nums}
.item .name{flex:1;font-weight:500}
.tier{font-size:11px;padding:1px 6px;border-radius:9px;background:#eee;color:#444;white-space:nowrap}
.dot{width:10px;height:10px;border-radius:50%;flex:none;border:1px solid rgba(0,0,0,.15)}
.dot.none{background:#fff;border:2px solid var(--none)}
.dot.good{background:var(--good)}.dot.possible{background:var(--possible)}.dot.unlikely{background:var(--unlikely)}
.dot.unclear{background:var(--unclear)}
.ext{font-size:10px;padding:0 5px;border-radius:8px;border:1px solid #8a6db0;color:#6e2769;white-space:nowrap}
.extbox{display:flex;align-items:center;gap:6px;margin:0 0 10px;cursor:pointer;user-select:none}
.extbox input{width:16px;height:16px;cursor:pointer}
main{overflow:auto;padding:16px 20px 40px}
.gene-head{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px;
 margin-bottom:14px}
.gene-head h2{margin:0 0 2px;font-size:20px}
.gene-head .sub{color:var(--muted);margin-bottom:10px}
.statuses{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:10px}
.statuses button{display:flex;align-items:center;gap:6px}
.statuses button .key{font-size:11px;color:var(--muted);border:1px solid var(--line);border-radius:4px;padding:0 4px}
.statuses button.on{border-width:2px;font-weight:600}
.statuses button.on.good{border-color:var(--good);color:var(--good)}
.statuses button.on.possible{border-color:var(--possible);color:var(--possible)}
.statuses button.on.unlikely{border-color:var(--unlikely);color:var(--unlikely)}
.statuses button.on.unclear{border-color:var(--unclear);color:#8a6510}
textarea{width:100%;min-height:54px;font:inherit;padding:8px;border:1px solid var(--line);border-radius:6px;resize:vertical}
.nav{display:flex;gap:8px;margin-top:8px;align-items:center;color:var(--muted);font-size:12px}
.pages img{display:block;width:100%;max-width:1400px;margin:0 auto 14px;background:#fff;border:1px solid var(--line);
 border-radius:4px;cursor:zoom-in}
.pages img.zoom{width:200%;max-width:none;cursor:zoom-out}
.pages{overflow-x:auto}
.empty{color:var(--muted);padding:40px;text-align:center}
#print{display:none}
@media print{
 @page{size:letter landscape;margin:0.25in}
 body{background:#fff}
 #app{display:none}
 #print{display:block}
 .pp{break-before:page;page-break-before:always}
 .pp:first-child{break-before:auto;page-break-before:auto}
 .pp img{display:block;width:100%;max-height:7.9in;object-fit:contain}
 .pp.first img{max-height:7.2in}
 .rbar{font-size:11pt;border:1.5px solid #333;border-radius:6px;padding:6px 10px;margin-bottom:6px}
 .rbar b{margin-right:8px}
 table.sum{border-collapse:collapse;width:100%;font-size:10pt}
 table.sum th,table.sum td{border:1px solid #bbb;padding:4px 6px;text-align:left;vertical-align:top}
 table.sum th{background:#f0f0f0}
}
</style>
</head>
<body>
<div id="app">
 <header>
  <h1 id="title"></h1><span class="meta" id="meta"></span>
  <span class="spacer"></span>
  <span id="progress"></span>
  <button id="btnExportTsv" title="Download the statuses and notes as a TSV file">Export reviews (TSV)</button>
  <button id="btnImport" title="Load statuses and notes from an exported TSV file">Import reviews</button>
  <input type="file" id="fileImport" accept=".tsv,.txt" style="display:none">
  <button id="btnPdf" class="primary" title="Opens the print dialog; choose Save as PDF">Export PDF</button>
 </header>
 <aside>
  <h2>Summary</h2>
  <div id="navSummary"></div>
  <h2>Candidate genes</h2>
  <select id="filter" class="filter"></select>
  <div id="navGenes"></div>
  <h2>QC and expression</h2>
  <div id="navQc"></div>
 </aside>
 <main id="main"></main>
</div>
<div id="print"></div>
<script type="application/json" id="data">__DATA__</script>
<script>
(function(){
"use strict";
const D = JSON.parse(document.getElementById("data").textContent);
const ST = D.statuses, SKEYS = ST.map(s => s.key), SLABEL = Object.fromEntries(ST.map(s => [s.key, s.label]));
const STORE = "rnadx-review|" + D.sample + "|" + D.generated;
let reviews = {};
try { reviews = JSON.parse(localStorage.getItem(STORE) || "{}") || {}; } catch (e) { reviews = {}; }
function save(){ try { localStorage.setItem(STORE, JSON.stringify(reviews)); } catch (e) {} }
const DEF = D.noStatus.key;
SLABEL[DEF] = D.noStatus.label;
function rev(g){
  const r = reviews[g] || {};
  return {status: SKEYS.includes(r.status) ? r.status : DEF, external: !!r.external,
          notes: r.notes || "", updated: r.updated || ""};
}
function setRev(g, patch){
  const r = Object.assign({}, rev(g), patch, {updated: new Date().toISOString().slice(0, 19)});
  reviews[g] = r; save(); refreshNav(); refreshProgress();
}
const pagesBy = (kind, name) => D.pages.filter(p => p.kind === kind && (name === undefined || p.name === name));
const qcNames = [...new Set(D.pages.filter(p => p.kind === "qc").map(p => p.name))];
const el = (tag, attrs, ...kids) => {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === "class") e.className = v; else if (k === "text") e.textContent = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v); else e.setAttribute(k, v);
  }
  kids.forEach(c => c != null && e.append(c));
  return e;
};

document.getElementById("title").textContent = D.sample;
document.getElementById("meta").textContent = "sample type: " + D.sampleType + "  ·  generated " + D.generated.replace("T", " ");
document.title = D.sample + " review";

const filter = document.getElementById("filter");
filter.append(el("option", {value: "all", text: "All genes"}));
filter.append(el("option", {value: DEF, text: D.noStatus.label}));
ST.forEach(s => filter.append(el("option", {value: s.key, text: s.label})));
filter.append(el("option", {value: "external", text: "Reviewed externally"}));
filter.addEventListener("change", refreshNav);

let current = {kind: "overview"};

function navItem(label, sel, extra){
  const a = el("div", {class: "item" + (extra && extra.cls ? " " + extra.cls : ""), onclick: () => show(sel)});
  if (extra && extra.rank != null) a.append(el("span", {class: "rank", text: extra.rank}));
  if (extra && extra.status) a.append(el("span", {class: "dot " + extra.status, title: SLABEL[extra.status]}));
  a.append(el("span", {class: "name", text: label}));
  if (extra && extra.external) a.append(el("span", {class: "ext", text: "ext", title: "Reviewed externally"}));
  if (extra && extra.tier) a.append(el("span", {class: "tier", text: "tier " + extra.tier}));
  if (sel.kind === current.kind && sel.name === current.name) a.classList.add("active");
  return a;
}
function refreshNav(){
  const ns = document.getElementById("navSummary"); ns.replaceChildren(navItem("Overview", {kind: "overview"}));
  const ng = document.getElementById("navGenes"); ng.replaceChildren();
  const f = filter.value;
  const shown = D.genes.filter(g => matches(g, f));
  shown.forEach(g => {
    const r = rev(g.gene);
    ng.append(navItem(g.gene, {kind: "gene", name: g.gene},
      {rank: g.rank, tier: g.tier, status: r.status, external: r.external, cls: "st-" + r.status}));
  });
  if (!D.genes.length) ng.append(el("div", {class: "empty", text: "No candidate genes"}));
  else if (!shown.length) ng.append(el("div", {class: "empty", text: "No genes with this status"}));
  const nq = document.getElementById("navQc"); nq.replaceChildren();
  qcNames.forEach(n => nq.append(navItem(n, {kind: "qc", name: n})));
}
function matches(g, f){
  const r = rev(g.gene);
  return f === "all" || (f === "external" ? r.external : r.status === f);
}
function refreshProgress(){
  const n = D.genes.length, done = D.genes.filter(g => rev(g.gene).status !== DEF).length;
  document.getElementById("progress").textContent = n ? (done + " of " + n + " genes with a status") : "";
}
function pageImgs(list){
  const box = el("div", {class: "pages"});
  list.forEach(p => {
    const img = el("img", {src: p.src, alt: (p.name || p.kind) + " page", loading: "lazy"});
    img.addEventListener("click", () => img.classList.toggle("zoom"));
    box.append(img);
  });
  return box;
}
function geneHead(g){
  const r = rev(g.gene);
  const inh = g.inheritance ? ("inheritance: " + g.inheritance + (g.haploinsufficient && /AD|XLD/.test(g.inheritance) ? " (haploinsufficient)" : "")) : "";
  const head = el("div", {class: "gene-head"},
    el("h2", {text: g.gene}),
    el("div", {class: "sub", text: ["rank " + g.rank, "tier " + g.tier, g.hits, inh].filter(Boolean).join("   ·   ")}));
  const sb = el("div", {class: "statuses"});
  ST.forEach((s, i) => {
    const b = el("button", {class: s.key + (r.status === s.key ? " on" : ""), title: "Shortcut: " + (i + 1),
      onclick: () => { setRev(g.gene, {status: rev(g.gene).status === s.key ? DEF : s.key}); show(current, true); }},
      el("span", {class: "dot " + s.key}), el("span", {text: s.label}), el("span", {class: "key", text: String(i + 1)}));
    sb.append(b);
  });
  const cb = el("input", {type: "checkbox"});
  cb.checked = r.external;
  cb.addEventListener("change", () => setRev(g.gene, {external: cb.checked}));
  const ext = el("label", {class: "extbox", title: "Shortcut: e"}, cb, el("span", {text: "Reviewed externally"}));
  const ta = el("textarea", {placeholder: "Notes"});
  ta.value = r.notes || "";
  ta.addEventListener("input", () => setRev(g.gene, {notes: ta.value}));
  const idx = D.genes.findIndex(x => x.gene === g.gene);
  const nav = el("div", {class: "nav"},
    el("button", {text: "← Previous gene", onclick: () => step(-1)}),
    el("button", {text: "Next gene →", onclick: () => step(1)}),
    el("span", {text: "Gene " + (idx + 1) + " of " + D.genes.length + "   ·   keys: 1–4 set status (again to clear), e reviewed externally, ↑/↓ or j/k move, click a page to zoom"}));
  head.append(sb, ext, ta, nav);
  return head;
}
function show(sel, keepScroll){
  current = sel;
  const main = document.getElementById("main");
  const y = main.scrollTop;
  main.replaceChildren();
  if (sel.kind === "gene") {
    const g = D.genes.find(x => x.gene === sel.name);
    main.append(geneHead(g), pageImgs(pagesBy("gene", g.gene)));
  } else if (sel.kind === "qc") {
    main.append(pageImgs(pagesBy("qc", sel.name)));
  } else {
    main.append(pageImgs(pagesBy("overview")));
  }
  main.scrollTop = keepScroll ? y : 0;
  refreshNav();
}
function visibleGenes(){
  const f = filter.value;
  return D.genes.filter(g => matches(g, f) || (current.kind === "gene" && g.gene === current.name));
}
function step(d){
  const list = visibleGenes();
  if (!list.length) return;
  let i = current.kind === "gene" ? list.findIndex(g => g.gene === current.name) : -1;
  i = Math.min(Math.max(i + d, 0), list.length - 1);
  show({kind: "gene", name: list[i].gene});
}
document.addEventListener("keydown", e => {
  if (e.target.tagName === "TEXTAREA" || e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;
  if (e.metaKey || e.ctrlKey || e.altKey) return;
  if (e.key === "ArrowDown" || e.key === "j") { step(1); e.preventDefault(); }
  else if (e.key === "ArrowUp" || e.key === "k") { step(-1); e.preventDefault(); }
  else if (/^[1-9]$/.test(e.key) && +e.key <= SKEYS.length && current.kind === "gene") {
    const k = SKEYS[+e.key - 1];
    setRev(current.name, {status: rev(current.name).status === k ? DEF : k}); show(current, true);
  }
  else if (e.key === "e" && current.kind === "gene") {
    setRev(current.name, {external: !rev(current.name).external}); show(current, true);
  }
});

// ---- export / import reviews (TSV)
function clean(s){ return String(s == null ? "" : s).replace(/\t/g, " ").replace(/\r?\n/g, "\\n"); }
document.getElementById("btnExportTsv").addEventListener("click", () => {
  const rows = [["sample", "gene", "rank", "tier", "status", "status_label", "reviewed_externally", "notes", "updated"]];
  D.genes.forEach(g => { const r = rev(g.gene);
    rows.push([D.sample, g.gene, g.rank, g.tier, r.status, SLABEL[r.status], r.external ? "TRUE" : "FALSE",
               r.notes, r.updated].map(clean)); });
  const blob = new Blob([rows.map(r => r.join("\t")).join("\n") + "\n"], {type: "text/tab-separated-values"});
  const a = el("a", {href: URL.createObjectURL(blob), download: D.sample + "_reviews.tsv"});
  document.body.append(a); a.click(); a.remove();
});
const fileIn = document.getElementById("fileImport");
document.getElementById("btnImport").addEventListener("click", () => fileIn.click());
fileIn.addEventListener("change", () => {
  const f = fileIn.files[0]; if (!f) return;
  f.text().then(txt => {
    const lines = txt.split(/\r?\n/).filter(l => l.trim());
    const hdr = (lines.shift() || "").split("\t");
    const gi = hdr.indexOf("gene"), si = hdr.indexOf("status"), ni = hdr.indexOf("notes"), ui = hdr.indexOf("updated"),
          xi = hdr.indexOf("reviewed_externally");
    if (gi < 0 || si < 0) { alert("Not a reviews file: it needs 'gene' and 'status' columns."); return; }
    const known = new Set(D.genes.map(g => g.gene)); let n = 0;
    lines.forEach(l => { const c = l.split("\t"); const g = c[gi];
      if (!known.has(g)) return;
      reviews[g] = {status: SKEYS.includes(c[si]) ? c[si] : DEF,
                    external: xi >= 0 && /^(true|1|yes)$/i.test(c[xi] || ""),
                    notes: ni >= 0 ? (c[ni] || "").replace(/\\n/g, "\n") : "", updated: ui >= 0 ? (c[ui] || "") : ""};
      n++; });
    save(); show(current, true); refreshProgress();
    alert("Imported reviews for " + n + " gene(s).");
    fileIn.value = "";
  });
});

// ---- export PDF: every page in report order, each gene with its status and notes, via the print dialog
document.getElementById("btnPdf").addEventListener("click", () => {
  const pr = document.getElementById("print"); pr.replaceChildren();
  const sum = el("div", {class: "pp"}, el("h2", {text: D.sample + " — review (" + new Date().toISOString().slice(0, 10) + ")"}));
  const t = el("table", {class: "sum"});
  t.append(el("tr", {}, ...["Rank", "Gene", "Tier", "Status", "Reviewed externally", "Notes"].map(h => el("th", {text: h}))));
  D.genes.forEach(g => { const r = rev(g.gene);
    t.append(el("tr", {}, el("td", {text: g.rank}), el("td", {text: g.gene}), el("td", {text: g.tier}),
      el("td", {text: SLABEL[r.status]}), el("td", {text: r.external ? "yes" : ""}), el("td", {text: r.notes || ""}))); });
  if (!D.genes.length) t.append(el("tr", {}, el("td", {colspan: "6", text: "No candidate genes"})));
  sum.append(t); pr.append(sum);
  let lastGene = null;
  D.pages.forEach(p => {
    const box = el("div", {class: "pp"});
    if (p.kind === "gene" && p.name !== lastGene) {
      const r = rev(p.name);
      box.classList.add("first");
      box.append(el("div", {class: "rbar"}, el("b", {text: p.name + ": " + SLABEL[r.status]
                                                          + (r.external ? "  ·  reviewed externally" : "")}),
        el("span", {text: r.notes ? "Notes: " + r.notes : ""})));
    }
    if (p.kind === "gene") lastGene = p.name;
    box.append(el("img", {src: p.src}));
    pr.append(box);
  });
  const imgs = [...pr.querySelectorAll("img")];
  Promise.all(imgs.map(i => i.decode ? i.decode().catch(() => {}) : Promise.resolve()))
    .then(() => { window.print(); });
});
window.addEventListener("afterprint", () => document.getElementById("print").replaceChildren());

refreshProgress();
show({kind: "overview"});
})();
</script>
</body>
</html>
"""