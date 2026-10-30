# 1E: check the WGS VCF matches our genome build, normalize it, and pull out the panel regions

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


def fail(msg):
    print(f"\nERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def run(cmd):
    # pipefail, otherwise a bcftools error mid-pipe still exits 0
    if subprocess.run(["bash", "-o", "pipefail", "-c", cmd]).returncode != 0:
        fail(f"command failed (see bcftools messages above; a reference allele mismatch means the VCF "
             f"is not from the pipeline's genome build):\n  {cmd}")


def capture(args):
    return subprocess.run(args, stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout


def n_records(vcf):
    return int(capture(["bcftools", "index", "-n", vcf]).strip())


def write_bed(path, intervals):
    with open(path, "w") as fh:
        fh.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in intervals)


# header lines, contig lengths from the header, and the single sample name
def read_header(vcf):
    lines = capture(["bcftools", "view", "-h", vcf]).splitlines()
    lengths = {}
    for line in lines:
        m = HEADER_CONTIG.match(line)
        if not m:
            continue
        # contig lines look like ##contig=<ID=chr1,length=248956422,assembly=...>, just grab the key=value bits
        attrs = dict(kv.split("=", 1) for kv in re.findall(r"[A-Za-z_]+=[^,>]+", m.group(1)))
        if "ID" in attrs and "length" in attrs:
            lengths[attrs["ID"]] = int(attrs["length"])

    # last header line is #CHROM..., anything after FORMAT is a sample name
    samples = lines[-1].split("\t")[9:]
    if len(samples) != 1:
        fail(f"{vcf} contains {len(samples)} samples ({', '.join(samples) or 'none'}). "
             f"WGS VCFs must contain exactly one sample.")
    print(f"WGS sample in VCF: {samples[0]}")
    return lines, lengths


def scan_records(vcf):
    counts, filters = Counter(), set()
    # streaming this on purpose, WGS VCFs can be millions of records
    proc = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%FILTER\n", vcf],
                            stdout=subprocess.PIPE, universal_newlines=True)
    for line in proc.stdout:
        chrom, _, filt = line.rstrip("\n").partition("\t")
        counts[chrom] += 1
        # non-PASS filters used in records need a ##FILTER line or bcftools complains later
        filters.update(f for f in filt.split(";") if f not in ("", ".", "PASS"))
    if proc.wait() != 0:
        fail(f"bcftools could not read {vcf}. If it is bgzipped, check it is complete with: bgzip -t {vcf}")
    return counts, filters


# match up VCF contigs with FASTA contigs, and bail if it looks like the wrong build
def resolve_contigs(vcf, fai, counts, header_lengths):
    def key(name):  # chr1 == 1, chrM == MT
        k = re.sub(r"^chr", "", name, flags=re.IGNORECASE).upper()
        return "M" if k == "MT" else k

    fasta_main = {key(c): c for c in fai if MAIN_CONTIG.match(c)}
    mapping, missing = {}, []
    for c in counts:
        if c in fai:
            mapping[c] = c
        elif MAIN_CONTIG.match(c):
            if key(c) in fasta_main:
                mapping[c] = fasta_main[key(c)]
            else:
                missing.append(c)
        # decoys/alts/unplaced not in the FASTA are dropped, not an error

    errors = []
    hits = Counter(mapping.values())
    # e.g. a VCF with both chr1 and 1 in it, no way to know which one we're supposed to use
    collisions = sorted(c for c, f in mapping.items() if hits[f] > 1)
    if collisions:
        errors.append(f"contigs {', '.join(collisions)} in the VCF all map to the same genome FASTA contig")
    if missing:
        errors.append(f"main chromosome(s) {', '.join(sorted(missing))} in the VCF have no counterpart in the genome FASTA")
    # a length mismatch is how a hg19 VCF shows up against hg38
    mismatched = [f"{c} (VCF {header_lengths[c]}, FASTA {f} {fai[f]})"
                  for c, f in sorted(mapping.items())
                  if c in header_lengths and header_lengths[c] != fai[f]]
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


# panel BED -> padded, merged regions to pull out of the normalized VCF
def panel_intervals(bed, padding, fai):
    by_chrom, skipped = {}, set()
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            if chrom not in fai:
                skipped.add(chrom)
                continue
            # clip at 0 and the contig end, padding can run off either end for genes near telomeres
            by_chrom.setdefault(chrom, []).append((max(0, int(start) - padding), min(fai[chrom], int(end) + padding)))
    if skipped:
        print(f"WARNING: panel BED contig(s) not in the genome FASTA were skipped: {', '.join(sorted(skipped))}")

    order = {c: i for i, c in enumerate(fai)}
    merged = []
    for chrom in sorted(by_chrom, key=order.get):
        for s, e in sorted(by_chrom[chrom]):
            if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], e)
            else:
                merged.append([chrom, s, e])
    return merged


def main():
    print("\n\n\n******************************************************************************************")
    print("Preparing WGS VCF...")
    print("******************************************************************************************\n")

    args = parse_args()
    # basic checks on the input before anything slow
    if not os.path.isfile(args.vcf):
        fail(f"WGS VCF {args.vcf} not found.")
    # big WGS files often arrive partially copied; catch that before bcftools gives a cryptic error
    if args.vcf.endswith((".gz", ".bgz")):
        print("Checking the compressed VCF is complete...")
        if subprocess.run(["bgzip", "-t", args.vcf], stderr=subprocess.PIPE).returncode != 0:
            fail(f"{args.vcf} is truncated or corrupted (bgzip -t failed). Re-copy or re-download it.")

    # load everything we need and line the contigs up (this is where a wrong-build VCF dies)
    header_lines, header_lengths = read_header(args.vcf)
    if not os.path.isfile(args.genome + ".fai"):
        fail(f"FASTA index {args.genome}.fai not found. Create it with: samtools faidx {args.genome}")
    with open(args.genome + ".fai") as fh:
        fai = {f[0]: int(f[1]) for f in (line.split("\t") for line in fh)}
    counts, filters = scan_records(args.vcf)
    mapping = resolve_contigs(args.vcf, fai, counts, header_lengths)

    dropped = sorted(c for c in counts if c not in mapping)
    if dropped:
        print(f"Dropping {sum(counts[c] for c in dropped)} record(s) on {len(dropped)} contig(s) "
              f"not in the genome FASTA: {', '.join(dropped)}")
    else:
        print("No records on contigs outside the genome FASTA.")

    for path in (args.norm_vcf, args.panel_vcf):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    q = shlex.quote
    # intermediate files go in a temp dir next to the output and clean themselves up
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(args.norm_vcf))) as tmp:
        order = {c: i for i, c in enumerate(fai)}
        # keep contigs in FASTA order, bcftools sort follows the header order so this keeps everything consistent
        kept = sorted(mapping, key=lambda c: order[mapping[c]])
        keep_bed = os.path.join(tmp, "genome_contigs.bed")
        # whole-contig intervals, so -T below basically means 'only these contigs'
        write_bed(keep_bed, [(c, 0, fai[mapping[c]]) for c in kept])
        with open(os.path.join(tmp, "rename_chrs.txt"), "w") as fh:
            fh.writelines(f"{c}\t{mapping[c]}\n" for c in kept)

        # genome contig lines replace the originals; dropped contigs stay declared so bcftools can read
        # their records before -T filters them out
        header = [l for l in header_lines if not l.startswith("##contig=")]
        defined = {m.group(1) for m in (re.match(r"^##FILTER=<ID=([^,>]+)", l) for l in header) if m}
        extra = [f'##FILTER=<ID={f},Description="{f}">' for f in sorted(filters - defined)]
        extra += [f"##contig=<ID={c},length={fai[mapping[c]]}>" for c in kept]
        extra += [f"##contig=<ID={c}>" for c in dropped]
        with open(os.path.join(tmp, "header.txt"), "w") as fh:
            fh.write("\n".join(header[:-1] + extra + header[-1:]) + "\n")

        # step 1: whole-genome normalized VCF
        print(f"\nNormalizing against {args.genome} and keeping only variant genotypes...")
        # GT="alt" is applied twice: early to shrink the work, and again because splitting
        # multiallelics leaves records where the sample only has the other alt
        run(f"bcftools reheader -h {q(os.path.join(tmp, 'header.txt'))} {q(args.vcf)}"
            f" | bcftools view -T {q(keep_bed)} -i 'GT=\"alt\"' -Ou"
            f" | bcftools annotate --rename-chrs {q(os.path.join(tmp, 'rename_chrs.txt'))} -Ou"
            f" | bcftools sort -T {q(os.path.join(tmp, 'sort'))} -Ou"
            f" | bcftools norm -f {q(args.genome)} -m -both -Ou"
            f" | bcftools view -i 'GT=\"alt\" && ALT!=\"<NON_REF>\" && ALT!=\"<*>\"' -Oz -o {q(args.norm_vcf)}")
        run(f"tabix -f -p vcf {q(args.norm_vcf)}")
        print(f"Normalized VCF: {args.norm_vcf} ({n_records(args.norm_vcf)} records)")

        # step 2: cut that down to the panel
        intervals = panel_intervals(args.bed, args.padding, fai)
        if not intervals:
            fail(f"no panel BED regions in {args.bed} are on contigs in the genome FASTA.")
        panel_bed = os.path.join(tmp, "panel_padded.bed")
        write_bed(panel_bed, intervals)
        print(f"\nRestricting to {len(intervals)} panel region(s) padded by {args.padding} bp...")
        run(f"bcftools view -R {q(panel_bed)} -Oz -o {q(args.panel_vcf)} {q(args.norm_vcf)}")
        run(f"tabix -f -p vcf {q(args.panel_vcf)}")
        print(f"Panel VCF: {args.panel_vcf} ({n_records(args.panel_vcf)} records)")

    print("\nDone.")


if __name__ == "__main__":
    main()
