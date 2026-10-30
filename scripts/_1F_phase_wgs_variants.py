# 1F - phase the WGS het sites. WGS phasing is treated as truth, the NanoTS RNA blocks are only used
# to join WGS blocks together or phase sites WGS left unphased

import argparse
import os
import subprocess
import sys
from collections import Counter, defaultdict

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add phase information to WGS heterozygous variants. WGS genotypes and WGS phasing are "
                    "treated as ground truth; NanoTS phase blocks are only used to phase WGS heterozygous "
                    "variants that WGS left unphased and to join WGS phase blocks. A NanoTS phase block that "
                    "disagrees with WGS phasing anywhere is not used at all.")
    parser.add_argument("--wgs-vcf", required=True,
        help="Normalized, panel-restricted WGS VCF from _1E_prepare_wgs_vcf.py.")
    parser.add_argument("--nanoTS-vcf", required=True,
        help="Normalized NanoTS VCF for the same sample.")
    parser.add_argument("--outfile", required=True,
        help="Output phased WGS VCF (.vcf.gz).")
    return parser.parse_args()


def fail(message):
    print(f"\nERROR: {message}", file=sys.stderr)
    sys.exit(1)


def het_alt_hap(sample):
    gt = sample.get("GT")
    if gt is None or len(gt) != 2 or None in gt or sorted(gt) != [0, 1]:
        return None
    return gt.index(1)


def ps_value(sample):
    try:
        ps = sample.get("PS")
    except (KeyError, ValueError):
        return None
    if ps in (None, ".", "", 0, "0"):
        return None
    return ps


def variant_key(rec):
    return (rec.chrom, rec.pos, rec.ref, rec.alts[0])


# union-find that also tracks parity, i.e. are two nodes on the same or opposite haplotype
class PhaseGraph:

    def __init__(self):
        self.parent = {}
        self.parity = {}

    def add(self, node):
        if node not in self.parent:
            self.parent[node] = node
            self.parity[node] = 0

    def find(self, node):
        path = []
        while self.parent[node] != node:
            path.append(node)
            node = self.parent[node]
        root, total = node, 0
        for n in reversed(path):
            total ^= self.parity[n]
            self.parity[n] = total
            self.parent[n] = root
        return root

    def relation(self, node):
        root = self.find(node)
        return root, self.parity[node] if node != root else 0

    def union(self, a, b, flip):
        ra, pa = self.relation(a)
        rb, pb = self.relation(b)
        if ra == rb:
            return (pa ^ pb) == flip
        self.parent[rb] = ra
        self.parity[rb] = pa ^ pb ^ flip
        return True

    def copy(self):
        other = PhaseGraph()
        other.parent = dict(self.parent)
        other.parity = dict(self.parity)
        return other


def load_wgs(path):
    sites = {}
    with pysam.VariantFile(path) as vcf:
        if len(vcf.header.samples) != 1:
            fail(f"{path} must contain exactly one sample.")
        for rec in vcf.fetch():
            if not rec.alts:
                continue
            s = rec.samples[0]
            hap = het_alt_hap(s)
            if hap is None:
                continue
            ps = ps_value(s) if s.phased else None
            if s.phased and ps is None:
                ps = f"chromosome:{rec.chrom}"
            sites[variant_key(rec)] = dict(hap=hap, wgs_block=(rec.chrom, ps) if s.phased else None)
    return sites


def load_nanoTS_blocks(path, wgs_sites, stats):
    blocks = defaultdict(dict)
    with pysam.VariantFile(path) as vcf:
        for rec in vcf.fetch():
            if not rec.alts:
                continue
            key = variant_key(rec)
            if key not in wgs_sites:
                continue
            s = rec.samples[0]
            hap = het_alt_hap(s)
            if hap is None:
                stats["WGS het sites not het in NanoTS (ignored)"] += 1
                continue
            ps = ps_value(s)
            if not s.phased or ps is None:
                continue
            blocks[(rec.chrom, ps)][key] = hap
    return blocks


# phased sites are represented by their WGS block, unphased ones stand alone
def site_node(key, site):
    return ("wgs_block",) + site["wgs_block"] if site["wgs_block"] else ("site",) + key


def build_wgs_graph(wgs_sites):
    graph = PhaseGraph()
    for key, site in wgs_sites.items():
        graph.add(site_node(key, site))
    return graph


def site_hap_in_node(site):
    return site["hap"] if site["wgs_block"] else 0


# returns None if the NanoTS block disagrees with WGS anywhere
def block_flips(block, wgs_sites, graph):
    flips = {}
    for key, nanoTS_hap in block.items():
        site = wgs_sites[key]
        root, parity = graph.relation(site_node(key, site))
        flip = nanoTS_hap ^ site_hap_in_node(site) ^ parity
        if flips.setdefault(root, flip) != flip:
            return None
    return flips


def apply_block(block, wgs_sites, graph):
    flips = block_flips(block, wgs_sites, graph)
    if flips is None:
        return False
    roots = list(flips)
    for other in roots[1:]:
        graph.union(roots[0], other, flips[roots[0]] ^ flips[other])
    return True


# add NanoTS blocks one by one, skipping any that contradict WGS or a block we already added
def harmonize(wgs_sites, blocks, stats):
    graph = build_wgs_graph(wgs_sites)
    wgs_only = graph.copy()
    ordered = sorted(blocks.items(), key=lambda kv: min(k[1] for k in kv[1]))
    for (chrom, ps), block in ordered:
        if len(block) < 2:
            stats["NanoTS blocks with <2 WGS het sites (not informative)"] += 1
            continue
        if block_flips(block, wgs_sites, wgs_only) is None:
            stats["NanoTS blocks disagreeing with WGS phasing (not used)"] += 1
            print(f"  NanoTS block {chrom}:{ps} disagrees with WGS phasing; not used.")
            continue
        trial = graph.copy()
        if not apply_block(block, wgs_sites, trial):
            stats["NanoTS blocks disagreeing with other NanoTS blocks (not used)"] += 1
            print(f"  NanoTS block {chrom}:{ps} disagrees with phasing from other NanoTS blocks; not used.")
            continue
        roots_before = {graph.find(site_node(k, wgs_sites[k])) for k in block}
        graph = trial
        stats["NanoTS blocks used"] += 1
        wgs_roots = {r for r in roots_before if r[0] == "wgs_block"}
        if len(wgs_roots) > 1:
            stats["WGS phase blocks joined by NanoTS"] += len(wgs_roots) - 1
    return graph


def assign_phase(wgs_sites, graph):
    members = defaultdict(list)
    for key, site in wgs_sites.items():
        root, parity = graph.relation(site_node(key, site))
        members[root].append((key, site_hap_in_node(site) ^ parity))
    phase = {}
    for root, sites in members.items():
        if len(sites) < 2 and root[0] == "site":
            continue
        by_chrom = defaultdict(list)
        for key, hap in sites:
            by_chrom[key[0]].append((key, hap))
        for chrom_sites in by_chrom.values():
            ps = min(k[1] for k, _ in chrom_sites)
            for key, hap in chrom_sites:
                phase[key] = (hap, ps)
    return phase


def write_output(wgs_path, outfile, wgs_sites, phase, stats):
    with pysam.VariantFile(wgs_path) as vcf_in:
        header = vcf_in.header.copy()
        if "PS" not in header.formats:
            header.formats.add("PS", 1, "Integer", "Phase set")
        header.formats.add("PHSRC", 1, "String",
                           "Source of the phase for heterozygous variants: WGS (phased in the WGS VCF) or "
                           "NanoTS (unphased in the WGS VCF, phased from a NanoTS phase block that agrees with WGS)")
        ps_is_int = header.formats["PS"].type == "Integer"
        tmp_out = outfile + ".tmp.vcf.gz"
        with pysam.VariantFile(tmp_out, "wz", header=header) as vcf_out:
            for rec in vcf_in.fetch():
                rec.translate(header)
                dst = rec.samples[0]
                key = variant_key(rec) if rec.alts else None
                if key in phase:
                    hap, ps = phase[key]
                    dst["GT"] = (1, 0) if hap == 0 else (0, 1)
                    dst.phased = True
                    dst["PS"] = ps if ps_is_int else str(ps)
                    source = "WGS" if wgs_sites[key]["wgs_block"] else "NanoTS"
                    dst["PHSRC"] = source
                    stats[f"WGS het sites phased ({source})"] += 1
                elif key in wgs_sites:
                    stats["WGS het sites left unphased"] += 1
                vcf_out.write(rec)
    os.replace(tmp_out, outfile)
    subprocess.run(["tabix", "-f", "-p", "vcf", outfile], check=True)


def main():
    print("\n\n\n******************************************************************************************")
    print("Phasing WGS variants...")
    print("******************************************************************************************\n")

    args = parse_args()
    for path in (args.wgs_vcf, args.nanoTS_vcf):
        if not os.path.isfile(path):
            fail(f"{path} not found.")

    # load, merge, write, then print a summary
    stats = Counter()
    wgs_sites = load_wgs(args.wgs_vcf)
    stats["WGS het sites"] = len(wgs_sites)
    stats["WGS het sites phased in WGS VCF"] = sum(1 for s in wgs_sites.values() if s["wgs_block"])
    blocks = load_nanoTS_blocks(args.nanoTS_vcf, wgs_sites, stats)
    stats["NanoTS phase blocks overlapping WGS het sites"] = len(blocks)

    graph = harmonize(wgs_sites, blocks, stats)
    phase = assign_phase(wgs_sites, graph)

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    write_output(args.wgs_vcf, args.outfile, wgs_sites, phase, stats)

    print("\nSummary:")
    for label in ("WGS het sites", "WGS het sites phased in WGS VCF",
                  "NanoTS phase blocks overlapping WGS het sites", "NanoTS blocks used",
                  "NanoTS blocks disagreeing with WGS phasing (not used)",
                  "NanoTS blocks disagreeing with other NanoTS blocks (not used)",
                  "NanoTS blocks with <2 WGS het sites (not informative)",
                  "WGS phase blocks joined by NanoTS", "WGS het sites not het in NanoTS (ignored)",
                  "WGS het sites phased (WGS)", "WGS het sites phased (NanoTS)", "WGS het sites left unphased"):
        print(f"  {label}: {stats[label]}")
    print(f"\nPhased WGS VCF: {args.outfile}")


if __name__ == "__main__":
    main()
