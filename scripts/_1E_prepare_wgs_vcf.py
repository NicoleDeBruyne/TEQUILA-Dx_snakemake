import argparse
import os
import re
import shlex
import subprocess
import sys
import tempfile
from collections import Counter

MAIN_CONTIG = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")
HEADER_CONTIG = re.compile(r"^##contig=<(.*)>$")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate and normalize a single-sample WGS VCF against the pipeline genome, then "
                    "restrict it to the panel genes plus flanking sequence.")
    parser.add_argument("--vcf", required=True,
        help="Input WGS VCF (vcf, vcf.gz or bcf). Must contain exactly one sample.")
    parser.add_argument("--genome", required=True,
        help="Genome FASTA used by the rest of the pipeline. A .fai index must exist next to it.")
    parser.add_argument("--bed", required=True,
        help="Panel BED with one gene region per row (chrom, start, end, gene).")
    parser.add_argument("--padding", type=int, default=10000,
        help="Bases added to each side of every panel gene region.")
    parser.add_argument("--norm-vcf", required=True,
        help="Output path for the whole-genome normalized VCF (.vcf.gz).")
    parser.add_argument("--panel-vcf", required=True,
        help="Output path for the panel-restricted normalized VCF (.vcf.gz).")
    return parser.parse_args()


def fail(message):
    print(f"\nERROR: {message}", file=sys.stderr)
    sys.exit(1)


def run(cmd):
    if subprocess.run(["bash", "-o", "pipefail", "-c", cmd]).returncode != 0:
        fail(f"command failed (see bcftools messages above; a reference allele mismatch means the VCF "
             f"is not from the pipeline's genome build):\n  {cmd}")


def capture(args):
    return subprocess.run(args, stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout


def read_fai(genome):
    fai = genome + ".fai"
    if not os.path.isfile(fai):
        fail(f"FASTA index {fai} not found. Create it with: samtools faidx {genome}")
    lengths = {}
    with open(fai) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def header_contig_lengths(vcf):
    lengths = {}
    for line in capture(["bcftools", "view", "-h", vcf]).splitlines():
        m = HEADER_CONTIG.match(line)
        if not m:
            continue
        attrs = dict(kv.split("=", 1) for kv in re.findall(r"[A-Za-z_]+=[^,>]+", m.group(1)))
        if "ID" in attrs and "length" in attrs:
            lengths[attrs["ID"]] = int(attrs["length"])
    return lengths


def scan_records(vcf):
    counts = Counter()
    filters = set()
    proc = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%FILTER\n", vcf],
                            stdout=subprocess.PIPE, universal_newlines=True)
    for line in proc.stdout:
        chrom, _, filt = line.rstrip("\n").partition("\t")
        counts[chrom] += 1
        filters.update(f for f in filt.split(";") if f not in ("", ".", "PASS"))
    if proc.wait() != 0:
        fail(f"bcftools could not read {vcf}. If it is bgzipped, check it is complete with: bgzip -t {vcf}")
    return counts, filters


def check_single_sample(vcf):
    samples = [s for s in capture(["bcftools", "query", "-l", vcf]).splitlines() if s]
    if len(samples) != 1:
        fail(f"{vcf} contains {len(samples)} samples ({', '.join(samples) or 'none'}). "
             f"WGS VCFs must contain exactly one sample.")
    print(f"WGS sample in VCF: {samples[0]}")


def contig_key(name):
    key = re.sub(r"^chr", "", name, flags=re.IGNORECASE).upper()
    return "M" if key == "MT" else key


def resolve_contigs(vcf, fai_lengths, record_counts, header_lengths):
    fasta_main = {contig_key(c): c for c in fai_lengths if MAIN_CONTIG.match(c)}
    mapping, missing = {}, []
    for c in record_counts:
        if c in fai_lengths:
            mapping[c] = c
        elif MAIN_CONTIG.match(c):
            if contig_key(c) in fasta_main:
                mapping[c] = fasta_main[contig_key(c)]
            else:
                missing.append(c)
    errors = []
    targets = Counter(mapping.values())
    collisions = sorted(c for c, f in mapping.items() if targets[f] > 1)
    if collisions:
        errors.append(f"contigs {', '.join(collisions)} in the VCF all map to the same genome FASTA contig")
    if missing:
        errors.append(f"main chromosome(s) {', '.join(sorted(missing))} in the VCF have no counterpart in the genome FASTA")
    mismatched = [f"{c} (VCF {header_lengths[c]}, FASTA {mapping[c]} {fai_lengths[mapping[c]]})"
                  for c in sorted(mapping)
                  if c in header_lengths and header_lengths[c] != fai_lengths[mapping[c]]]
    if mismatched:
        errors.append(f"chromosome length(s) differ from the genome FASTA, so the VCF is likely from "
                      f"a different genome build: {'; '.join(mismatched)}")
    if not any(MAIN_CONTIG.match(c) for c in mapping):
        errors.append("no main chromosome in the VCF matches the genome FASTA")
    if errors:
        fail(f"{vcf} does not match the genome:\n  - " + "\n  - ".join(errors))
    renamed = {c: f for c, f in mapping.items() if c != f}
    if renamed:
        print(f"Renaming {len(renamed)} contig(s) to match the genome FASTA: "
              f"{', '.join(f'{c}->{f}' for c, f in sorted(renamed.items()))}")
    if not any(c in header_lengths for c in mapping):
        print("WARNING: the VCF header has no contig lengths, so the genome build could not be checked by length. "
              "Normalization will still fail on reference allele mismatches.")
    print(f"Main chromosomes in VCF: {sum(1 for c in mapping if MAIN_CONTIG.match(c))}, all matched to the genome FASTA.")
    return mapping


def write_header(vcf, path, contigs, other_contigs, record_filters):
    lines = [l for l in capture(["bcftools", "view", "-h", vcf]).splitlines() if not l.startswith("##contig=")]
    defined_filters = {m.group(1) for l in lines for m in [re.match(r"^##FILTER=<ID=([^,>]+)", l)] if m}
    filter_lines = [f'##FILTER=<ID={f},Description="{f}">' for f in sorted(record_filters - defined_filters)]
    contig_lines = [f"##contig=<ID={c},length={n}>" for c, n in contigs] + [f"##contig=<ID={c}>" for c in other_contigs]
    with open(path, "w") as fh:
        fh.write("\n".join(lines[:-1] + filter_lines + contig_lines + lines[-1:]) + "\n")


def write_bed(path, intervals):
    with open(path, "w") as fh:
        for chrom, start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\n")


def padded_panel_intervals(bed, padding, fai_lengths):
    by_chrom = {}
    skipped = set()
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in fai_lengths:
                skipped.add(chrom)
                continue
            by_chrom.setdefault(chrom, []).append(
                (max(0, start - padding), min(fai_lengths[chrom], end + padding)))
    if skipped:
        print(f"WARNING: panel BED contig(s) not in the genome FASTA were skipped: {', '.join(sorted(skipped))}")
    order = {c: i for i, c in enumerate(fai_lengths)}
    merged = []
    for chrom in sorted(by_chrom, key=order.get):
        cur_start, cur_end = None, None
        for start, end in sorted(by_chrom[chrom]):
            if cur_end is not None and start <= cur_end:
                cur_end = max(cur_end, end)
                continue
            if cur_end is not None:
                merged.append((chrom, cur_start, cur_end))
            cur_start, cur_end = start, end
        if cur_end is not None:
            merged.append((chrom, cur_start, cur_end))
    return merged


def count_records(vcf):
    return int(capture(["bcftools", "index", "-n", vcf]).strip())


def main():
    print(f"\n\n\n******************************************************************************************")
    print(f"Preparing WGS VCF...")
    print(f"******************************************************************************************\n")

    args = parse_args()
    if not os.path.isfile(args.vcf):
        fail(f"WGS VCF {args.vcf} not found.")

    if args.vcf.endswith((".gz", ".bgz")):
        print("Checking the compressed VCF is complete...")
        if subprocess.run(["bgzip", "-t", args.vcf], stderr=subprocess.PIPE, universal_newlines=True).returncode != 0:
            fail(f"{args.vcf} is truncated or corrupted (bgzip -t failed). Re-copy or re-download it.")
    check_single_sample(args.vcf)
    fai_lengths = read_fai(args.genome)
    record_counts, filters = scan_records(args.vcf)
    header_lengths = header_contig_lengths(args.vcf)
    mapping = resolve_contigs(args.vcf, fai_lengths, record_counts, header_lengths)

    dropped = {c: n for c, n in record_counts.items() if c not in mapping}
    if dropped:
        print(f"Dropping {sum(dropped.values())} record(s) on {len(dropped)} contig(s) not in the genome FASTA: "
              f"{', '.join(sorted(dropped))}")
    else:
        print("No records on contigs outside the genome FASTA.")

    os.makedirs(os.path.dirname(os.path.abspath(args.norm_vcf)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.panel_vcf)), exist_ok=True)
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(args.norm_vcf))) as tmp:
        order = {c: i for i, c in enumerate(fai_lengths)}
        kept = sorted(mapping, key=lambda c: order[mapping[c]])
        keep_bed = os.path.join(tmp, "genome_contigs.bed")
        write_bed(keep_bed, [(c, 0, fai_lengths[mapping[c]]) for c in kept])
        header = os.path.join(tmp, "header.txt")
        write_header(args.vcf, header, [(c, fai_lengths[mapping[c]]) for c in kept],
                     sorted(c for c in record_counts if c not in mapping), filters)
        rename_file = os.path.join(tmp, "rename_chrs.txt")
        with open(rename_file, "w") as fh:
            for c in kept:
                fh.write(f"{c}\t{mapping[c]}\n")

        print(f"\nNormalizing against {args.genome} and keeping only variant genotypes...")
        run(f"bcftools reheader -h {shlex.quote(header)} {shlex.quote(args.vcf)}"
            f" | bcftools view -T {shlex.quote(keep_bed)} -i 'GT=\"alt\"' -Ou"
            f" | bcftools annotate --rename-chrs {shlex.quote(rename_file)} -Ou"
            f" | bcftools sort -T {shlex.quote(os.path.join(tmp, 'sort'))} -Ou"
            f" | bcftools norm -f {shlex.quote(args.genome)} -m -both -Ou"
            f" | bcftools view -i 'GT=\"alt\" && ALT!=\"<NON_REF>\" && ALT!=\"<*>\"' -Oz -o {shlex.quote(args.norm_vcf)}")
        run(f"tabix -f -p vcf {shlex.quote(args.norm_vcf)}")
        print(f"Normalized VCF: {args.norm_vcf} ({count_records(args.norm_vcf)} records)")

        panel_bed = os.path.join(tmp, "panel_padded.bed")
        intervals = padded_panel_intervals(args.bed, args.padding, fai_lengths)
        if not intervals:
            fail(f"no panel BED regions in {args.bed} are on contigs in the genome FASTA.")
        write_bed(panel_bed, intervals)
        print(f"\nRestricting to {len(intervals)} panel region(s) padded by {args.padding} bp...")
        run(f"bcftools view -R {shlex.quote(panel_bed)} -Oz -o {shlex.quote(args.panel_vcf)} {shlex.quote(args.norm_vcf)}")
        run(f"tabix -f -p vcf {shlex.quote(args.panel_vcf)}")
        print(f"Panel VCF: {args.panel_vcf} ({count_records(args.panel_vcf)} records)")

    print("\nDone.")


if __name__ == "__main__":
    main()
