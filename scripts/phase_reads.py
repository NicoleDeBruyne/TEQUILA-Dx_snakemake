#!/usr/bin/env python3

import argparse
import os
import subprocess
import pysam
import numpy as np
import shutil
import concurrent.futures
import traceback


def parse_args():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Phases reads from an RNA-seq BAM file that cover gene regions of interest.")
    parser.add_argument("--bam", type=str, required=True)
    parser.add_argument("--nanoTS-vcf", type=str, required=True,
        help="NanoTS VCF (already internally phased -- has its own GT/PS). Per gene, the trusted "
             "variant set is selected directly from here: every PASS/DP-filtered variant overlapping "
             "the gene region, PLUS every PASS/DP-filtered variant anywhere else in this VCF that "
             "shares a phase set (PS) with one of those -- see select_nanoTS_variants_for_gene(). No AF "
             "filter is applied to NanoTS -- its genotypes are trusted directly, not used as a "
             "confidence-scored candidate the way Clair3/DeepVariant indels still are.")
    parser.add_argument("--clair3-vcf", type=str, required=True,
        help="Clair3 VCF, used only for indels shared with --deepvariant-vcf that NanoTS didn't already "
             "call -- see select_candidate_indels_for_gene().")
    parser.add_argument("--deepvariant-vcf", type=str, required=True,
        help="DeepVariant VCF, used the same way as --clair3-vcf (indel agreement partner).")
    parser.add_argument("--region", type=str)
    parser.add_argument("--bed", type=str)
    parser.add_argument("--genome", type=str, required=True)
    parser.add_argument("--outdir", type=str, required=True)
    parser.add_argument("--name", type=str, default="SAMPLE")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--snvs-only', action='store_true',
        help="Exclude NanoTS indels from the trusted set (SNVs only). Candidate indels "
             "(Clair3/DeepVariant-shared) are unaffected by this flag.")
    parser.add_argument('--min-dp', type=int, default=20,
        help="Minimum DP for a NanoTS or candidate-indel variant to be selected at all. Default: 20")
    parser.add_argument('--min-af', type=float, default=0.1,
        help="Minimum AF for a candidate indel (Clair3/DeepVariant-shared) to be selected. NOT applied "
             "to NanoTS variants. Default: 0.1")
    parser.add_argument('--terminal-variant-proportion', type=float, default=0.5)
    parser.add_argument('--min-distance-from-read-end', type=int, default=20)
    parser.add_argument('--phasing-threshold', type=float, default=0.5)
    parser.add_argument('--remove-monoexonic', action='store_true')
    parser.add_argument('--ignore-variants-list', type=str, nargs='+')
    parser.add_argument('--ignore-variants-bed', type=str)
    parser.add_argument('--samtools-exec', type=str, default='samtools')
    parser.add_argument('--bcftools-exec', type=str, default='bcftools')
    parser.add_argument('--tabix-exec', type=str, default='tabix')
    parser.add_argument('--bgzip-exec', type=str, default='bgzip')
    parser.add_argument('--minimap2-exec', type=str, default='minimap2')
    parser.add_argument('--whatshap-exec', type=str, default='whatshap')
    return parser.parse_args()


########################################################################################################################
# Helper functions for extracting gene regions
########################################################################################################################

def extract_gene_regions(bed):
    """Extract gene regions from the BED file in the format: {gene: (chrom, start, end)}."""
    gene_regions = {}

    with open(bed, 'r') as b:
        for line in b:
            fields = line.strip().split("\t")
            gene = fields[3]
            gene_regions[gene] = (fields[0], int(fields[1]), int(fields[2]))

    return gene_regions


########################################################################################################################
# Helper functions for phasing reads
########################################################################################################################

def format_variant(key):
    """Format a variant tuple for human-readable reports."""
    chrom, pos, ref, alt = key
    return f"{chrom}:{pos} {ref}>{alt}"


def write_variant_list(report_file, variants):
    """Write one variant per line."""
    for key in sorted(variants):
        report_file.write(f"    {format_variant(key)}\n")


def filter_bam_by_region(inbam, outbam, region, threads):
    """Filter a BAM file for primary and supplementary alignments over a region."""
    subprocess.run(
        [
            samtools_exec,
            'view',
            '-F', '256',
            '-@', str(threads),
            '-b',
            '-o', outbam,
            inbam,
            region
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )

    subprocess.run(
        [samtools_exec, 'index', '-@', str(threads), outbam],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )


def remove_monoexonic_reads(inbam, outbam, threads=1):
    """Remove reads lacking an 'N' CIGAR operation and write spliced reads to a new BAM."""
    with pysam.AlignmentFile(inbam, "rb") as infile, \
         pysam.AlignmentFile(outbam, "wb", header=infile.header, threads=threads) as outfile:

        for read in infile.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            if read.cigartuples and any(op == 3 for op, _ in read.cigartuples):
                outfile.write(read)

    sorted_bam = outbam.replace(".bam", ".sorted.bam")

    subprocess.run(
        [samtools_exec, "sort", "-@", str(threads), "-o", sorted_bam, outbam],
        check=True
    )

    os.replace(sorted_bam, outbam)

    subprocess.run(
        [samtools_exec, "index", "-@", str(threads), outbam],
        check=True
    )

    return outbam


def _load_ignore_positions(ignore_variants_bed):
    """Returns a set of (chrom, pos) 1-based positions to exclude from variant selection."""
    if not ignore_variants_bed:
        return None

    positions = set()

    with open(ignore_variants_bed) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])

            for p in range(start + 1, end + 1):
                positions.add((chrom, p))

    return positions


def write_variants_vcf(outfile, variants, sample_name, contigs):
    """Writes a minimal VCF for candidate indels."""
    if outfile.endswith(".gz"):
        outfile = outfile[:-3]

    sorted_variants = sorted(variants.items())

    with open(outfile, "w") as f:
        f.write("##fileformat=VCFv4.2\n")

        for c in contigs:
            f.write(f"##contig=<ID={c}>\n")

        f.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        f.write(
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
        )
        f.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )

        for (chrom, pos, ref, alt), (GT, DP) in sorted_variants:
            gt_str = "/".join(map(str, GT)) if GT else "./."

            f.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP\t"
                f"{gt_str}:{DP}\n"
            )

    subprocess.run(
        [bcftools_exec, "sort", "-Oz", "-o", f"{outfile}.gz", outfile],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )

    os.remove(outfile)

    subprocess.run(
        [tabix_exec, "-f", "-p", "vcf", f"{outfile}.gz"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )


def write_nanoTS_vcf(outfile, variants, sample_name, contigs):
    """Writes a minimal VCF for the trusted NanoTS set."""
    if outfile.endswith(".gz"):
        outfile = outfile[:-3]

    sorted_variants = sorted(variants.items())

    with open(outfile, "w") as f:
        f.write("##fileformat=VCFv4.2\n")

        for c in contigs:
            f.write(f"##contig=<ID={c}>\n")

        f.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        f.write(
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
        )
        f.write(
            '##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set">\n'
        )
        f.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )

        for (chrom, pos, ref, alt), (GT, phased, PS, DP) in sorted_variants:

            if not GT:
                gt_str = "./."
            else:
                sep = "|" if phased else "/"
                gt_str = sep.join(map(str, GT))

            ps_str = str(PS) if PS not in (None, "") else "."

            f.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP:PS\t"
                f"{gt_str}:{DP}:{ps_str}\n"
            )

    subprocess.run(
        [bcftools_exec, "sort", "-Oz", "-o", f"{outfile}.gz", outfile],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )

    os.remove(outfile)

    subprocess.run(
        [tabix_exec, "-f", "-p", "vcf", f"{outfile}.gz"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )


def read_nanoTS_vcf(vcf_path):
    """Read a NanoTS VCF into the internal variant dictionary."""
    variants = {}

    with pysam.VariantFile(vcf_path) as vcf:

        for rec in vcf.fetch():

            if not rec.alts:
                continue

            s = rec.samples[0]
            PS = s.get("PS")

            if PS in (".", "", 0, "0"):
                PS = None

            variants[
                (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            ] = (
                s.get("GT"),
                s.phased,
                PS,
                s.get("DP") or 0
            )

    return variants


def select_nanoTS_variants_for_gene(
    nanoTS_vcf_path,
    chrom,
    start,
    end,
    min_dp,
    snvs_only=False,
    ignore_positions=None
):
    """Select trusted heterozygous NanoTS variants."""

    def _passes(rec):

        if rec.filter.keys() != ["PASS"]:
            return None

        if len(rec.alts) != 1:
            return None

        ref, alt = rec.ref, rec.alts[0]

        if snvs_only and len(ref) == 1 and len(alt) == 1:
            pass
        elif snvs_only:
            return None

        if ignore_positions and (rec.chrom, rec.pos) in ignore_positions:
            return None

        s = rec.samples[0]
        GT = s.get("GT")

        if not GT or None in GT or len(set(GT)) < 2:
            return None

        DP = s.get("DP") or 0

        if DP < min_dp:
            return None

        return s

    variants = {}
    phase_sets = set()

    with pysam.VariantFile(nanoTS_vcf_path) as vcf:

        for rec in vcf.fetch(chrom, start, end):

            s = _passes(rec)

            if s is None:
                continue

            key = (
                rec.chrom,
                rec.pos,
                rec.ref,
                rec.alts[0]
            )

            PS = s.get("PS")

            variants[key] = (
                s.get("GT"),
                s.phased,
                PS,
                s.get("DP") or 0
            )

            if PS not in (None, ".", 0, "0"):
                phase_sets.add(PS)

        if phase_sets:

            for rec in vcf.fetch():

                if not rec.alts:
                    continue

                key = (
                    rec.chrom,
                    rec.pos,
                    rec.ref,
                    rec.alts[0]
                )

                if key in variants:
                    continue

                s = _passes(rec)

                if s is None:
                    continue

                if s.get("PS") in phase_sets:

                    variants[key] = (
                        s.get("GT"),
                        s.phased,
                        s.get("PS"),
                        s.get("DP") or 0
                    )

    return variants


def select_candidate_indels_for_gene(
    clair3_vcf_path,
    deepvariant_vcf_path,
    chrom,
    start,
    end,
    min_dp,
    min_af,
    exclude_keys,
    ignore_positions=None
):
    """Select candidate indels shared by Clair3 and DeepVariant."""

    def _get(vcf_path):

        out = {}

        with pysam.VariantFile(vcf_path) as vcf:

            for rec in vcf.fetch(chrom, start, end):

                if rec.filter.keys() != ["PASS"]:
                    continue

                if len(rec.alts) != 1:
                    continue

                ref, alt = rec.ref, rec.alts[0]

                if len(ref) == 1 and len(alt) == 1:
                    continue

                if ignore_positions and (rec.chrom, rec.pos) in ignore_positions:
                    continue

                s = rec.samples[0]

                DP = s.get("DP") or 0
                AF = s.get("AF") or s.get("VAF") or (0,)

                if isinstance(AF, (list, tuple)):
                    AF = AF[0]

                if DP < min_dp or AF < min_af:
                    continue

                out[
                    (rec.chrom, rec.pos, ref, alt)
                ] = (
                    s.get("GT"),
                    DP
                )

        return out

    clair3 = _get(clair3_vcf_path)
    deepvariant = _get(deepvariant_vcf_path)

    shared = {
        k: clair3[k]
        for k in clair3.keys() & deepvariant.keys()
    }

    return {
        k: v
        for k, v in shared.items()
        if k not in exclude_keys
    }


def remove_end_variants(
    invcf,
    outvcf,
    bam,
    min_distance_from_read_end=10,
    terminal_variant_proportion=0.5
):
    """Remove variants whose ALT-supporting reads frequently occur near read ends."""

    if outvcf.endswith('.gz'):
        outvcf = outvcf[:-3]

    bamfile = pysam.AlignmentFile(bam, "rb")
    vcf_reader = pysam.VariantFile(invcf, "r")
    vcf_writer = pysam.VariantFile(
        outvcf,
        "w",
        header=vcf_reader.header
    )

    for record in vcf_reader:

        chrom = record.chrom
        pos = record.pos - 1
        ref = record.ref
        alts = record.alts

        keep_record = True

        for alt in alts:

            total_alt_alignments = 0
            end_alt_alignments = 0

            for alignment in bamfile.fetch(
                chrom,
                pos,
                pos + len(alt)
            ):

                if alignment.is_unmapped:
                    continue

                aligned_pairs = alignment.get_aligned_pairs(
                    matches_only=True
                )

                query_pos = None

                for qpos, rpos in aligned_pairs:

                    if rpos == pos:
                        query_pos = qpos
                        break

                if query_pos is None:
                    continue

                query_seq = alignment.query_sequence[
                    query_pos:query_pos + len(alt)
                ]

                if query_seq != alt:
                    continue

                total_alt_alignments += 1

                q_positions = [
                    qpos
                    for qpos, _ in aligned_pairs
                ]

                if not q_positions:
                    continue

                q_start = q_positions[0]
                q_end = q_positions[-1]

                distance_to_end = min(
                    query_pos - q_start,
                    q_end - (
                        query_pos +
                        max(len(ref), len(alt))
                    )
                )

                if distance_to_end < min_distance_from_read_end:
                    end_alt_alignments += 1

            if (
                total_alt_alignments == 0
                or (
                    end_alt_alignments /
                    total_alt_alignments
                ) > terminal_variant_proportion
            ):
                keep_record = False
                break

        if keep_record:
            vcf_writer.write(record)

    bamfile.close()
    vcf_reader.close()
    vcf_writer.close()

    subprocess.run(
        [
            bcftools_exec,
            'sort',
            '-Oz',
            '-o',
            f"{outvcf}.gz",
            outvcf
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )

    os.remove(outvcf)

    subprocess.run(
        [
            tabix_exec,
            '-f',
            '-p',
            'vcf',
            f"{outvcf}.gz"
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )


########################################################################################################################
# Direct allele-based support
########################################################################################################################

def _get_read_allele(read, pos0, ref, alt):
    """
    Determine whether a read directly supports REF, ALT, or neither.

    Returns:
        "REF"
        "ALT"
        None
    """

    if read.is_unmapped:
        return None

    ##################################################################
    # SNV / MNP
    ##################################################################

    if len(ref) == len(alt):

        aligned_pairs = read.get_aligned_pairs(
            matches_only=False
        )

        query_positions = {}

        for qpos, rpos in aligned_pairs:

            if (
                rpos is not None
                and pos0 <= rpos < pos0 + len(ref)
            ):
                query_positions[rpos] = qpos

        if any(
            p not in query_positions
            for p in range(pos0, pos0 + len(ref))
        ):
            return None

        observed = []

        for p in range(pos0, pos0 + len(ref)):

            qpos = query_positions[p]

            if qpos is None:
                return None

            observed.append(
                read.query_sequence[qpos].upper()
            )

        observed = "".join(observed)

        if observed == ref.upper():
            return "REF"

        if observed == alt.upper():
            return "ALT"

        return None

    ##################################################################
    # Indel
    ##################################################################

    return _classify_indel_read(read, pos0, ref, alt)


def _classify_indel_read(read, pos0, ref, alt):
    """
    Classify a read's support at an indel variant (VCF convention: the
    anchor base at pos0 is shared by REF and ALT, followed by the
    inserted/deleted bases), by walking read.cigartuples once and
    inspecting the actual CIGAR operation(s) immediately after the anchor.

    This deliberately does NOT reproduce the original raw-sequence-only
    heuristic (compare query_sequence[qpos:qpos+n] to ALT), which had two
    real accuracy problems, found and confirmed against real CIGAR
    structures:

      - Deletions: the original check only ever returned "ALT" when the
        read was truncated right at the anchor base, and returned "REF"
        for a read carrying the exact deletion CIGAR operation -- an
        inversion, not a rounding choice. This version instead looks at
        whether a real D operation of exactly the right length starts
        immediately after the anchor.

      - Insertions: the original check compared raw query bases to ALT
        without confirming an actual "I" CIGAR operation was present at
        that position, so a read with no insertion at all could
        coincidentally match ALT's sequence and be misclassified. This
        version requires an actual I operation of the right length before
        trusting the sequence comparison.

    Returns "REF", "ALT", or None (ambiguous / not confidently either --
    e.g. the read ends right at the anchor, or a different-length/-type
    event occupies the position instead of a clean match to either allele).
    """

    if pos0 < read.reference_start or pos0 >= read.reference_end:
        return None

    REF_CONSUMING = (0, 2, 3, 7, 8)    # M, D, N, =, X
    QUERY_CONSUMING = (0, 1, 4, 7, 8)  # M, I, S, =, X
    MATCH_OPS = (0, 7, 8)              # M, =, X

    ops = read.cigartuples
    n = len(ops)

    rpos = read.reference_start
    qpos = 0
    op_i = 0

    ##################################################################
    # Locate the CIGAR operation covering the anchor position, pos0.
    ##################################################################

    while op_i < n:

        op, length = ops[op_i]
        op_ref_len = length if op in REF_CONSUMING else 0

        if op_ref_len and rpos <= pos0 < rpos + op_ref_len:
            break

        rpos += op_ref_len
        qpos += (length if op in QUERY_CONSUMING else 0)
        op_i += 1

    else:
        return None

    op, length = ops[op_i]

    if op not in MATCH_OPS:
        # The anchor itself sits inside a deletion/skip -- can't classify.
        return None

    anchor_qpos = qpos + (pos0 - rpos)
    query_base = read.query_sequence[anchor_qpos].upper()

    if query_base != ref[0].upper():
        return None

    is_last_base_of_block = (pos0 == rpos + length - 1)

    ##################################################################
    # Insertion
    ##################################################################

    if len(alt) > len(ref):

        insertion_length = len(alt) - len(ref)

        if not is_last_base_of_block:
            # More matched bases follow before any indel -- no insertion here.
            return "REF"

        if op_i + 1 >= n:
            # Read ends exactly at the anchor -- can't confirm either way.
            return None

        next_op, next_len = ops[op_i + 1]

        if next_op == 1 and next_len == insertion_length:  # I

            observed = read.query_sequence[
                anchor_qpos: anchor_qpos + insertion_length + 1
            ].upper()

            if observed == alt.upper():
                return "ALT"

            return None  # right-length insertion, wrong sequence

        if next_op in MATCH_OPS:
            return "REF"

        return None  # a different-length indel or other event sits here

    ##################################################################
    # Deletion
    ##################################################################

    deletion_length = len(ref) - len(alt)

    def _matched_coverage_confirms_ref(start_op_i, remaining_needed):
        """Walk forward confirming `remaining_needed` more ref-consuming
        bases are covered by plain matches (M/=/X), with no D/N/I
        interrupting -- i.e. the deletion is genuinely absent here."""

        i2 = start_op_i

        while remaining_needed > 0 and i2 < n:

            op2, len2 = ops[i2]

            if op2 in MATCH_OPS:
                remaining_needed -= len2
            elif op2 in (2, 3, 1):  # D, N, or I -- conflicts with a clean REF
                return False

            i2 += 1

        return remaining_needed <= 0

    if not is_last_base_of_block:

        matched_so_far = (rpos + length - 1) - pos0

        if matched_so_far >= deletion_length:
            return "REF"

        if _matched_coverage_confirms_ref(
            op_i + 1,
            deletion_length - matched_so_far
        ):
            return "REF"

        return None

    if op_i + 1 >= n:
        return None

    next_op, next_len = ops[op_i + 1]

    if next_op == 2 and next_len == deletion_length:  # D
        return "ALT"

    if next_op in MATCH_OPS:

        if _matched_coverage_confirms_ref(op_i + 1, deletion_length):
            return "REF"

        return None

    return None


def compute_variant_support(bam_path, variant):
    """
    Compute coverage and direct REF/ALT support for one variant.

    This is the single source of truth for variant-level read support.

    Returns a dictionary containing:
        coverage
        ref_reads
        alt_reads
        unassigned_reads
    """

    chrom, pos, ref, alt = variant
    pos0 = pos - 1

    ref_reads = set()
    alt_reads = set()
    unassigned_reads = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        for read in bam.fetch(
            chrom,
            pos0,
            pos0 + max(1, len(ref))
        ):

            if read.is_unmapped:
                continue

            read_name = read.query_name

            allele = _get_read_allele(
                read,
                pos0,
                ref,
                alt
            )

            if allele == "REF":
                ref_reads.add(read_name)

            elif allele == "ALT":
                alt_reads.add(read_name)

            else:
                unassigned_reads.add(read_name)

    ##################################################################
    # A read name may occur more than once because supplementary
    # alignments are retained. A read supporting either allele should
    # not also be counted as unassigned.
    ##################################################################

    unassigned_reads -= ref_reads
    unassigned_reads -= alt_reads

    overlap = ref_reads & alt_reads

    ref_reads -= overlap
    alt_reads -= overlap

    coverage = (
        len(ref_reads) +
        len(alt_reads) +
        len(unassigned_reads)
    )

    return {
        "coverage": coverage,
        "ref_reads": ref_reads,
        "alt_reads": alt_reads,
        "unassigned_reads": unassigned_reads,
    }


def compute_all_variant_support(bam_path, variants):
    """
    Compute coverage and REF/ALT support for every variant.

    Returns:
        {
            variant_key: {
                "coverage": int,
                "ref_reads": set,
                "alt_reads": set,
                "unassigned_reads": set
            }
        }
    """

    support = {}

    for variant in variants:

        support[variant] = compute_variant_support(
            bam_path,
            variant
        )

    return support


def _is_snv(key):
    """A variant key is a SNV iff both REF and ALT are single bases."""
    _, _, ref, alt = key
    return len(ref) == 1 and len(alt) == 1


def evaluate_snv(
    bam_path,
    variant,
    min_distance_from_read_end,
    terminal_variant_proportion,
    check_end_bias
):
    """
    Single BAM pass for one SNV that computes, together:

      (a) REF/ALT/unassigned read support -- identical semantics to
          compute_variant_support()/_get_read_allele()'s SNV branch, and
      (b) whether the variant should be dropped for having its ALT support
          concentrated too close to read ends -- identical semantics to
          remove_end_variants().

    This only ever runs for SNVs (len(ref) == len(alt) == 1), because that
    is the one case where compute_variant_support()'s fetch window
    (pos0, pos0 + max(1, len(ref))) and remove_end_variants()' fetch window
    (pos, pos + len(alt)) are the same interval (pos0, pos0 + 1), so the two
    original passes can be safely collapsed into one without changing which
    reads either function would have seen. Indels keep the original,
    unmerged code path unchanged (see select_candidate_indels_for_gene /
    compute_all_variant_support / remove_end_variants), since their fetch
    windows differ and their allele classification still requires a real
    CIGAR walk.

    Reads are classified using pysam.pileup(), which resolves the CIGAR at
    the single reference column in pysam's C layer, instead of each read's
    full get_aligned_pairs() being materialized in Python as the original
    _get_read_allele()/remove_end_variants() did.

    Returns:
        (support, drop_for_end_bias)
        support is the same shape compute_variant_support() returns:
            {"coverage", "ref_reads", "alt_reads", "unassigned_reads"}
    """

    chrom, pos, ref, alt = variant
    pos0 = pos - 1
    ref_u, alt_u = ref.upper(), alt.upper()

    ref_reads = set()
    alt_reads = set()
    unassigned_reads = set()

    total_alt_alignments = 0
    end_alt_alignments = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        for pileupcolumn in bam.pileup(
            chrom,
            pos0,
            pos0 + 1,
            truncate=True,
            min_base_quality=0,
            flag_filter=0,
            ignore_overlaps=False,
        ):

            for pread in pileupcolumn.pileups:

                aln = pread.alignment

                if aln.is_unmapped:
                    continue

                read_name = aln.query_name

                if (
                    pread.is_del
                    or pread.is_refskip
                    or pread.query_position is None
                ):
                    unassigned_reads.add(read_name)
                    continue

                observed = aln.query_sequence[
                    pread.query_position
                ].upper()

                if observed == ref_u:
                    ref_reads.add(read_name)
                    continue

                if observed == alt_u:

                    alt_reads.add(read_name)

                    if check_end_bias:

                        total_alt_alignments += 1

                        query_pos = pread.query_position
                        q_start = aln.query_alignment_start
                        q_end = aln.query_alignment_end - 1

                        distance_to_end = min(
                            query_pos - q_start,
                            q_end - (
                                query_pos +
                                max(len(ref), len(alt))
                            )
                        )

                        if distance_to_end < min_distance_from_read_end:
                            end_alt_alignments += 1

                    continue

                unassigned_reads.add(read_name)

    ##################################################################
    # Same overlap-removal / dedup rules as compute_variant_support().
    ##################################################################

    unassigned_reads -= ref_reads
    unassigned_reads -= alt_reads

    overlap = ref_reads & alt_reads

    ref_reads -= overlap
    alt_reads -= overlap

    coverage = (
        len(ref_reads) +
        len(alt_reads) +
        len(unassigned_reads)
    )

    support = {
        "coverage": coverage,
        "ref_reads": ref_reads,
        "alt_reads": alt_reads,
        "unassigned_reads": unassigned_reads,
    }

    drop_for_end_bias = False

    if check_end_bias:

        if (
            total_alt_alignments == 0
            or (
                end_alt_alignments /
                total_alt_alignments
            ) > terminal_variant_proportion
        ):
            drop_for_end_bias = True

    return support, drop_for_end_bias


def write_direct_haplotag(
    gene_outdir,
    ref_reads,
    alt_reads,
    unassigned_reads
):
    """Write a Whatshap-compatible haplotag list from direct allele assignments."""

    haplotag_outfile = os.path.join(
        gene_outdir,
        'whatshap_haplotag.tsv'
    )

    with open(haplotag_outfile, 'w') as f:

        f.write("#read_name\thaplotype\n")

        for read_name in sorted(ref_reads):
            f.write(
                f"{read_name}\tH1\n"
            )

        for read_name in sorted(alt_reads):
            f.write(
                f"{read_name}\tH2\n"
            )

        for read_name in sorted(unassigned_reads):
            f.write(
                f"{read_name}\tnone\n"
            )

    return haplotag_outfile


########################################################################################################################
# Whatshap / coverage helpers
########################################################################################################################

def get_phased_coverage(
    bam_path,
    hap1_bam_path,
    hap2_bam_path,
    region
):
    """Get the max coverage and max phased coverage safely."""

    chrom, positions = region.split(":")
    start, end = map(
        int,
        positions.split("-")
    )

    def max_cov(bam_file):

        with pysam.AlignmentFile(
            bam_file,
            "rb"
        ) as bam:

            A, C, G, T = bam.count_coverage(
                chrom,
                start,
                end,
                quality_threshold=0
            )

            return int(
                np.max(
                    np.add(
                        np.add(A, C),
                        np.add(G, T)
                    )
                )
            )

    total_max = max_cov(bam_path)

    phased_max = (
        max_cov(hap1_bam_path) +
        max_cov(hap2_bam_path)
    )

    return total_max, phased_max


def run_haplotag(
    phased_vcf,
    filtered_bam,
    genome,
    gene_outdir
):
    """Run Whatshap haplotag and return H1/H2/unassigned read names."""

    haplotag_outfile = os.path.join(
        gene_outdir,
        'whatshap_haplotag.tsv'
    )

    subprocess.run(
        [
            whatshap_exec,
            "haplotag",
            "--reference",
            genome,
            "--output",
            "/dev/null",
            "--ignore-read-groups",
            "--output-haplotag-list",
            haplotag_outfile,
            phased_vcf,
            filtered_bam
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )

    hap1_reads = []
    hap2_reads = []
    unassigned_reads = []

    with open(haplotag_outfile) as f:

        for line in f:

            if line.startswith("#") or not line.strip():
                continue

            cols = line.strip().split("\t")

            if len(cols) < 2:
                continue

            hap = cols[1]

            if hap == "H1":
                hap1_reads.append(cols[0])

            elif hap == "H2":
                hap2_reads.append(cols[0])

            else:
                unassigned_reads.append(cols[0])

    return (
        hap1_reads,
        hap2_reads,
        unassigned_reads
    )


def finalize_haplotype_bams(
    filtered_bam,
    gene_outdir,
    name,
    gene,
    region,
    phasing_threshold,
    threads,
    remove_monoexonic,
    phased_vcf,
    hap1_reads,
    hap2_reads,
    unassigned_reads
):
    """Write and QC the final haplotype-specific BAMs."""

    report_message = ""

    haplotag_outfile = os.path.join(
        gene_outdir,
        'whatshap_haplotag.tsv'
    )

    hap1_bam = os.path.join(
        gene_outdir,
        f'{name}_{gene}_hap1.bam'
    )

    hap2_bam = os.path.join(
        gene_outdir,
        f'{name}_{gene}_hap2.bam'
    )

    unassigned_bam = os.path.join(
        gene_outdir,
        f'{name}_{gene}_unassigned.bam'
    )

    subprocess.run(
        [
            whatshap_exec,
            "split",
            "--output-h1",
            hap1_bam,
            "--output-h2",
            hap2_bam,
            "--output-untagged",
            unassigned_bam,
            filtered_bam,
            haplotag_outfile
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True
    )

    for bam_path in (
        hap1_bam,
        hap2_bam,
        unassigned_bam
    ):

        subprocess.run(
            [
                samtools_exec,
                'index',
                '-@',
                str(threads),
                bam_path
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )

    num_phased_reads = (
        len(hap1_reads) +
        len(hap2_reads)
    )

    max_coverage, max_phased_coverage = get_phased_coverage(
        filtered_bam,
        hap1_bam,
        hap2_bam,
        region
    )

    report_message += (
        f"\n    {num_phased_reads} reads can be assigned to haplotypes."
        f"\n      {len(unassigned_reads)} reads can not be assigned to a haplotype."
        f"\n    Max phased coverage: {max_phased_coverage}"
        f"\n    Max overall coverage: {max_coverage}"
    )

    if (
        max_coverage == 0
        or
        max_phased_coverage <
        phasing_threshold * max_coverage
    ):

        for fp in [
            hap1_bam,
            hap1_bam + ".bai",
            hap2_bam,
            hap2_bam + ".bai",
            unassigned_bam,
            unassigned_bam + ".bai"
        ]:

            if os.path.exists(fp):
                os.remove(fp)

        report_message += (
            f"\nCannot phase {phasing_threshold*100}% of the max read coverage "
            f"over {gene} ({region}). Exiting without phasing {name}..."
        )

        return report_message, None

    report_message += (
        f"\n\nSuccessfully phased {num_phased_reads} out of "
        f"{num_phased_reads + len(unassigned_reads)} reads over "
        f"{gene} ({region})."
    )

    if remove_monoexonic:

        multidir = os.path.join(
            gene_outdir,
            "multiexonic_bams"
        )

        for src_bam in (
            hap1_bam,
            hap2_bam
        ):

            dest = os.path.join(
                multidir,
                os.path.splitext(
                    os.path.basename(src_bam)
                )[0] +
                "_multiexonic.bam"
            )

            remove_monoexonic_reads(
                src_bam,
                dest
            )

    summary_row = (
        gene,
        len(hap1_reads),
        len(hap2_reads),
        len(unassigned_reads),
        max_coverage,
        max_phased_coverage
    )

    with open(
        os.path.join(
            gene_outdir,
            f"{name}_{gene}_haplotype_assignment_summary.txt"
        ),
        "w"
    ) as f:

        f.write(
            "gene\thap1_read_count\thap2_read_count\t"
            "unassigned_read_count\tmax_coverage\tmax_phased_coverage\n"
        )

        f.write(
            f"{gene}\t"
            f"{len(hap1_reads)}\t"
            f"{len(hap2_reads)}\t"
            f"{len(unassigned_reads)}\t"
            f"{max_coverage}\t"
            f"{max_phased_coverage}\n"
        )

    return report_message, summary_row


########################################################################################################################
# Main function for phasing reads
########################################################################################################################

def phase_reads(
    bam,
    gene,
    region,
    nanoTS_vcf,
    clair3_vcf,
    deepvariant_vcf,
    min_dp,
    min_af,
    min_distance_from_read_end,
    terminal_variant_proportion,
    ignore_variants_bed,
    snvs_only,
    phasing_threshold,
    name,
    genome,
    outdir,
    threads,
    remove_monoexonic
):

    gene_outdir = os.path.join(
        outdir,
        gene
    )

    tempdir = os.path.join(
        gene_outdir,
        'temp'
    )

    os.makedirs(
        tempdir,
        exist_ok=True
    )

    if remove_monoexonic:

        os.makedirs(
            os.path.join(
                gene_outdir,
                'multiexonic_bams'
            ),
            exist_ok=True
        )

    chrom, positions = region.split(":")
    region_start, region_end = map(
        int,
        positions.split("-")
    )

    ignore_positions = _load_ignore_positions(
        ignore_variants_bed
    )

    report = os.path.join(
        gene_outdir,
        f'{name}_{gene}_phasing_report.txt'
    )

    with open(report, 'w') as report_file:

        ##################################################################
        # Filter BAM
        ##################################################################

        report_file.write(
            f"\nFiltering BAM file for reads aligned to "
            f"{gene} ({region})...\n"
        )
        report_file.flush()

        filtered_bam = os.path.join(
            gene_outdir,
            f'{name}_{gene}.bam'
        )

        filter_bam_by_region(
            bam,
            filtered_bam,
            region,
            threads
        )

        if remove_monoexonic:

            multidir = os.path.join(
                gene_outdir,
                'multiexonic_bams'
            )

            remove_monoexonic_reads(
                filtered_bam,
                os.path.join(
                    multidir,
                    os.path.splitext(
                        os.path.basename(filtered_bam)
                    )[0] +
                    "_multiexonic.bam"
                )
            )

        with pysam.AlignmentFile(filtered_bam, "rb") as _count_bam:
            total_read_count = _count_bam.count(until_eof=True)

        if total_read_count == 0:

            report_file.write(
                f"    No reads aligned to {gene} ({region}). "
                f"Exiting without phasing {name}...\n"
            )

            report_file.write(
                f"\nFinished processing {name} over "
                f"{gene} ({region})."
            )

            report_file.flush()

            os.remove(filtered_bam)
            os.remove(f"{filtered_bam}.bai")
            shutil.rmtree(tempdir)

            return None

        report_file.write(
            f"    {total_read_count} primary/supplementary alignments found."
        )
        report_file.flush()

        ##################################################################
        # Select NanoTS variants
        ##################################################################

        report_file.write(
            f"\n\nSelecting trusted, heterozygous NanoTS variants "
            f"for {gene} ({region})..."
        )
        report_file.flush()

        nanoTS_variants = select_nanoTS_variants_for_gene(
            nanoTS_vcf,
            chrom,
            region_start,
            region_end,
            min_dp,
            snvs_only,
            ignore_positions
        )

        report_file.write(
            f"\n    {len(nanoTS_variants)} trusted NanoTS variant(s) selected "
            f"(overlapping the gene region, or sharing a phase set with one "
            f"that does)."
        )

        if nanoTS_variants:

            write_variant_list(
                report_file,
                nanoTS_variants.keys()
            )

        report_file.flush()

        ##################################################################
        # Select candidate indels
        ##################################################################

        report_file.write(
            f"\nSelecting candidate (Clair3/DeepVariant-shared) indels "
            f"for {gene} ({region})..."
        )
        report_file.flush()

        candidate_indels = select_candidate_indels_for_gene(
            clair3_vcf,
            deepvariant_vcf,
            chrom,
            region_start,
            region_end,
            min_dp,
            min_af,
            set(nanoTS_variants.keys()),
            ignore_positions
        )

        report_file.write(
            f"\n    {len(candidate_indels)} candidate indel(s) selected "
            f"(not already in the trusted NanoTS set)."
        )

        if candidate_indels:

            write_variant_list(
                report_file,
                candidate_indels.keys()
            )

        report_file.flush()

        ##################################################################
        # End-of-read filtering
        #
        # SNVs and indels are now handled by two different code paths that
        # are required to reach the *same* filtering decision as the
        # original single VCF-round-trip implementation:
        #
        #   - indels keep the original remove_end_variants() VCF round trip,
        #     unchanged (their fetch window differs from
        #     compute_variant_support()'s, and their allele classification
        #     genuinely needs a CIGAR walk, so there is nothing safe to
        #     collapse here).
        #   - SNVs are filtered directly against the in-memory dict, using
        #     evaluate_snv() (see its docstring) instead of a temp-VCF +
        #     bcftools sort/tabix + remove_end_variants() + VCF-reread
        #     round trip. This also computes their REF/ALT support in the
        #     same pass, so SNVs skip the separate
        #     compute_all_variant_support() call entirely, below.
        ##################################################################

        end_filter_enabled = (
            min_distance_from_read_end > 0
            and terminal_variant_proportion < 1
        )

        if end_filter_enabled:

            report_file.write(
                f"\nRemoving variants that occur within "
                f"{min_distance_from_read_end}nt of the end of a read "
                f"more than {terminal_variant_proportion*100}% of the time..."
            )
            report_file.flush()

        nanoTS_indel_variants = {
            k: v
            for k, v in nanoTS_variants.items()
            if not _is_snv(k)
        }

        nanoTS_snv_variants = {
            k: v
            for k, v in nanoTS_variants.items()
            if _is_snv(k)
        }

        ##################################################################
        # Indels (NanoTS indel subset + candidate indels): unchanged
        # VCF-round-trip end-filtering.
        ##################################################################

        if end_filter_enabled:

            if nanoTS_indel_variants:

                nanoTS_vcf = os.path.join(
                    tempdir,
                    f'{gene}_nanoTS.vcf.gz'
                )

                write_nanoTS_vcf(
                    nanoTS_vcf,
                    nanoTS_indel_variants,
                    name,
                    [chrom]
                )

                filtered_nanoTS_vcf = os.path.join(
                    tempdir,
                    f'{gene}_nanoTS_filtered.vcf.gz'
                )

                remove_end_variants(
                    nanoTS_vcf,
                    filtered_nanoTS_vcf,
                    filtered_bam,
                    min_distance_from_read_end,
                    terminal_variant_proportion
                )

                nanoTS_indel_variants = read_nanoTS_vcf(
                    filtered_nanoTS_vcf
                )

            if candidate_indels:

                candidate_vcf = os.path.join(
                    tempdir,
                    f'{gene}_candidate_indels.vcf.gz'
                )

                write_variants_vcf(
                    candidate_vcf,
                    candidate_indels,
                    name,
                    [chrom]
                )

                filtered_candidate_vcf = os.path.join(
                    tempdir,
                    f'{gene}_candidate_indels_filtered.vcf.gz'
                )

                remove_end_variants(
                    candidate_vcf,
                    filtered_candidate_vcf,
                    filtered_bam,
                    min_distance_from_read_end,
                    terminal_variant_proportion
                )

                filtered_candidate_indels = {}

                with pysam.VariantFile(
                    filtered_candidate_vcf
                ) as vcf:

                    for rec in vcf.fetch():

                        if not rec.alts:
                            continue

                        key = (
                            rec.chrom,
                            rec.pos,
                            rec.ref,
                            rec.alts[0]
                        )

                        filtered_candidate_indels[key] = (
                            rec.samples[0].get("GT"),
                            rec.samples[0].get("DP") or 0
                        )

                candidate_indels = filtered_candidate_indels

        ##################################################################
        # SNVs (NanoTS SNV subset): end-filtering and REF/ALT support in a
        # single pysam.pileup()-based pass per variant (evaluate_snv()).
        ##################################################################

        snv_support = {}
        surviving_snv_variants = {}

        for key, val in nanoTS_snv_variants.items():

            support, drop = evaluate_snv(
                filtered_bam,
                key,
                min_distance_from_read_end,
                terminal_variant_proportion,
                end_filter_enabled
            )

            if drop:
                continue

            surviving_snv_variants[key] = val
            snv_support[key] = support

        nanoTS_variants = {
            **nanoTS_indel_variants,
            **surviving_snv_variants
        }

        ##################################################################
        # Report surviving variants
        ##################################################################

        report_file.write(
            f"\n    {len(nanoTS_variants)} NanoTS variant(s) remain."
        )

        if nanoTS_variants:

            write_variant_list(
                report_file,
                nanoTS_variants.keys()
            )

        report_file.write(
            f"    {len(candidate_indels)} candidate indel(s) remain."
        )

        if candidate_indels:

            write_variant_list(
                report_file,
                candidate_indels.keys()
            )

        report_file.flush()

        ##################################################################
        # Compute REF/ALT allele support for every surviving variant.
        # SNV support was already computed above (snv_support); only the
        # remaining indel keys (NanoTS indel survivors + candidate indels)
        # still need compute_all_variant_support()'s CIGAR-walk path.
        ##################################################################

        pool_keys = (
            list(nanoTS_variants.keys()) +
            list(candidate_indels.keys())
        )

        indel_pool_keys = [
            k for k in pool_keys if not _is_snv(k)
        ]

        variant_support = dict(snv_support)

        if pool_keys:

            report_file.write(
                f"\nComputing REF/ALT allele support over "
                f"{len(pool_keys)} variant(s) for {gene}..."
            )

            report_file.flush()

            if indel_pool_keys:

                variant_support.update(
                    compute_all_variant_support(
                        filtered_bam,
                        indel_pool_keys
                    )
                )

            if nanoTS_variants:

                report_file.write(
                    "\n    NanoTS variants:"
                )

                for key in sorted(nanoTS_variants):

                    result = variant_support[key]

                    report_file.write(
                        f"\n      {format_variant(key)}"
                        f"  REF={len(result['ref_reads'])}"
                        f"  ALT={len(result['alt_reads'])}"
                    )

            if candidate_indels:

                report_file.write(
                    "\n    Candidate indels:"
                )

                for key in sorted(candidate_indels):

                    result = variant_support[key]

                    report_file.write(
                        f"\n      {format_variant(key)}"
                        f"  REF={len(result['ref_reads'])}"
                        f"  ALT={len(result['alt_reads'])}"
                    )

            ##################################################################
            # Select the variant with the highest number of allele-
            # supporting reads (REF + ALT) from the already-computed data.
            ##################################################################

            best_key = max(
                pool_keys,
                key=lambda k: (
                    len(variant_support[k]["ref_reads"]) +
                    len(variant_support[k]["alt_reads"])
                )
            )

            best_support = variant_support[best_key]
            best_dp = best_support["coverage"]

            source = (
                "NanoTS variant"
                if best_key in nanoTS_variants
                else "candidate indel"
            )

            report_file.write(
                f"\n\n    Highest number of allele-supporting reads: "
                f"a {source} at {format_variant(best_key)} "
                f"(REF={len(best_support['ref_reads'])}, "
                f"ALT={len(best_support['alt_reads'])})."
            )

            report_file.flush()

        else:
            best_key = None
            best_dp = 0

        ##################################################################
        # No variants = nothing to phase
        ##################################################################

        if best_key is None:

            report_file.write(
                f"\nNo NanoTS variants or candidate indels found for "
                f"{gene} ({region}); nothing to phase reads on. "
                f"Exiting without attempting to haplotag..."
            )

            report_file.write(
                f"\nFinished processing {name} over "
                f"{gene} ({region})."
            )

            report_file.flush()

            shutil.rmtree(tempdir)

            return None

        ##################################################################
        # Get already-computed allele support for the best variant
        ##################################################################

        best_support = variant_support[best_key]

        best_ref_reads = best_support["ref_reads"]
        best_alt_reads = best_support["alt_reads"]
        best_unassigned_reads = best_support["unassigned_reads"]

        allele_support = (
            len(best_ref_reads) +
            len(best_alt_reads)
        )

        ##################################################################
        # Determine whether a NanoTS phased block exists
        ##################################################################

        ps_counts = {}

        for _, phased_flag, PS, _ in nanoTS_variants.values():

            if PS is not None:
                ps_counts[PS] = (
                    ps_counts.get(PS, 0) +
                    1
                )

        phased_block_ps = {
            ps
            for ps, n in ps_counts.items()
            if n >= 2
        }

        ##################################################################
        # Helper for constructing a single-variant phased VCF
        ##################################################################

        def _single_variant_vcf():

            out = os.path.join(
                gene_outdir,
                f'{gene}_phased.vcf.gz'
            )

            chrom_b, pos_b, ref_b, alt_b = best_key

            write_nanoTS_vcf(
                out,
                {
                    best_key: (
                        (0, 1),
                        True,
                        str(pos_b),
                        best_dp
                    )
                },
                name,
                [chrom]
            )

            return out

        ##################################################################
        # If NanoTS has a phased block, try Whatshap first
        ##################################################################

        if phased_block_ps:

            block_variants = {
                k: v
                for k, v in nanoTS_variants.items()
                if v[2] in phased_block_ps
            }

            phased_vcf = os.path.join(
                gene_outdir,
                f'{gene}_phased.vcf.gz'
            )

            write_nanoTS_vcf(
                phased_vcf,
                block_variants,
                name,
                [chrom]
            )

            report_file.write(
                f"\nNanoTS phased block found "
                f"({len(block_variants)} variant(s), "
                f"{len(phased_block_ps)} phase set(s))."
            )

            report_file.write(
                "\nRunning Whatshap haplotag..."
            )

            report_file.flush()

            hap1_reads, hap2_reads, unassigned_reads = run_haplotag(
                phased_vcf,
                filtered_bam,
                genome,
                gene_outdir
            )

            n_block_phased = (
                len(hap1_reads) +
                len(hap2_reads)
            )

            report_file.write(
                f"\n\n    Whatshap haplotag result:"
                f"\n      H1: {len(hap1_reads)} reads"
                f"\n      H2: {len(hap2_reads)} reads"
                f"\n      Unassigned: {len(unassigned_reads)} reads"
            )

            report_file.write(
                f"\n\n    Whatshap assigned "
                f"{n_block_phased} reads."
            )

            report_file.write(
                f"\n    Highest-coverage variant has "
                f"{allele_support} allele-supporting reads."
            )

            report_file.flush()

            ##################################################################
            # Whatshap fallback
            ##################################################################

            if n_block_phased < allele_support:

                report_file.write(
                    f"\n\n    Fewer reads were assigned by Whatshap "
                    f"than support the highest-coverage variant."
                    f"\n    Falling back to direct allele-based "
                    f"phasing at {format_variant(best_key)}..."
                )

                report_file.flush()

                write_direct_haplotag(
                    gene_outdir,
                    best_ref_reads,
                    best_alt_reads,
                    best_unassigned_reads
                )

                hap1_reads = sorted(best_ref_reads)
                hap2_reads = sorted(best_alt_reads)
                unassigned_reads = sorted(
                    best_unassigned_reads
                )

                report_file.write(
                    f"\n\n    Direct allele-based phasing:"
                    f"\n      H1 ({best_key[2]}): "
                    f"{len(hap1_reads)} reads"
                    f"\n      H2 ({best_key[3]}): "
                    f"{len(hap2_reads)} reads"
                    f"\n      Unassigned: "
                    f"{len(unassigned_reads)} reads"
                )

                report_file.flush()

        ##################################################################
        # No NanoTS phased block: directly phase using the highest-
        # coverage surviving variant.
        ##################################################################

        else:

            report_file.write(
                "\n\nNo NanoTS phased block found."
            )

            report_file.write(
                f"\n\nAssigning reads directly using the "
                f"highest-coverage variant:"
                f"\n{format_variant(best_key)}"
            )

            report_file.flush()

            ##################################################################
            # Reuse the support already computed above.
            ##################################################################

            write_direct_haplotag(
                gene_outdir,
                best_ref_reads,
                best_alt_reads,
                best_unassigned_reads
            )

            hap1_reads = sorted(best_ref_reads)
            hap2_reads = sorted(best_alt_reads)
            unassigned_reads = sorted(
                best_unassigned_reads
            )

            report_file.write(
                f"\n\n    Direct allele-based phasing:"
                f"\n      H1 ({best_key[2]}): "
                f"{len(hap1_reads)} reads"
                f"\n      H2 ({best_key[3]}): "
                f"{len(hap2_reads)} reads"
                f"\n      Unassigned: "
                f"{len(unassigned_reads)} reads"
            )

            report_file.flush()

        ##################################################################
        # Create final haplotype BAMs
        ##################################################################

        report_file.write(
            f"\n\nCreating haplotype-specific BAM files..."
        )

        report_file.flush()

        report_message, summary_row = finalize_haplotype_bams(
            filtered_bam,
            gene_outdir,
            name,
            gene,
            region,
            phasing_threshold,
            threads,
            remove_monoexonic,
            phased_vcf if phased_block_ps else None,
            hap1_reads,
            hap2_reads,
            unassigned_reads
        )

        report_file.write(
            report_message
        )

        report_file.write(
            f"\nFinished processing {name} over "
            f"{gene} ({region})."
        )

        report_file.flush()

    shutil.rmtree(tempdir)

    return summary_row


########################################################################################################################
# Main script
########################################################################################################################

def main():
    """Main function."""

    print(
        "\n\n\n******************************************************************************************"
    )
    print(
        "Phasing reads from RNA-seq data..."
    )
    print(
        "******************************************************************************************\n"
    )

    args = parse_args()

    print(
        f"Preparing to phase reads for {args.name}...\n"
    )

    global samtools_exec, bcftools_exec, tabix_exec
    global bgzip_exec, minimap2_exec, whatshap_exec

    samtools_exec = args.samtools_exec
    bcftools_exec = args.bcftools_exec
    tabix_exec = args.tabix_exec
    bgzip_exec = args.bgzip_exec
    minimap2_exec = args.minimap2_exec
    whatshap_exec = args.whatshap_exec

    if not os.path.exists(args.bam):
        raise FileNotFoundError(
            f"BAM file {args.bam} not found."
        )

    if not os.path.exists(args.genome):
        raise FileNotFoundError(
            f"Genome file {args.genome} not found."
        )

    if args.bed and not os.path.exists(args.bed):
        raise FileNotFoundError(
            f"BED file {args.bed} not found."
        )

    if not os.path.exists(args.nanoTS_vcf):
        raise FileNotFoundError(
            f"NanoTS VCF file {args.nanoTS_vcf} not found."
        )

    if not os.path.exists(args.clair3_vcf):
        raise FileNotFoundError(
            f"Clair3 VCF file {args.clair3_vcf} not found."
        )

    if not os.path.exists(args.deepvariant_vcf):
        raise FileNotFoundError(
            f"DeepVariant VCF file {args.deepvariant_vcf} not found."
        )

    if (
        args.ignore_variants_bed
        and
        not os.path.exists(args.ignore_variants_bed)
    ):
        raise FileNotFoundError(
            f"BED file containing variants to ignore "
            f"({args.ignore_variants_bed}) not found."
        )

    if not args.region and not args.bed:
        raise ValueError(
            "Either --region or --bed must be specified."
        )

    os.makedirs(
        args.outdir,
        exist_ok=True
    )

    if args.ignore_variants_list:

        ignore_variants_bed = os.path.join(
            args.outdir,
            'ignore_variants.bed'
        )

        with open(
            ignore_variants_bed,
            'w'
        ) as f:

            for variant in args.ignore_variants_list:

                chrom, pos = variant.split(':')

                f.write(
                    f"{chrom}\t{int(pos)-1}\t{pos}\n"
                )

    elif args.ignore_variants_bed:

        ignore_variants_bed = args.ignore_variants_bed

    else:

        ignore_variants_bed = None

    snvs_only = bool(
        args.snvs_only
    )

    summary_rows = []

    if args.bed:

        gene_regions = extract_gene_regions(
            args.bed
        )

        threads = min(
            args.threads,
            len(gene_regions)
        )

        threads_per_process = max(
            1,
            args.threads //
            len(gene_regions)
        )

        print(
            f"\nBegin phasing reads in {len(gene_regions)} "
            f"gene regions of interest in parallel with "
            f"{threads} workers. This may take a while...\n"
        )

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=threads
        ) as executor:

            futures = {
                executor.submit(
                    phase_reads,
                    args.bam,
                    gene,
                    f"{chrom}:{start}-{end}",
                    args.nanoTS_vcf,
                    args.clair3_vcf,
                    args.deepvariant_vcf,
                    args.min_dp,
                    args.min_af,
                    args.min_distance_from_read_end,
                    args.terminal_variant_proportion,
                    ignore_variants_bed,
                    snvs_only,
                    args.phasing_threshold,
                    args.name,
                    args.genome,
                    args.outdir,
                    threads_per_process,
                    args.remove_monoexonic
                ): gene
                for gene, (chrom, start, end)
                in gene_regions.items()
            }

            for future in concurrent.futures.as_completed(
                futures
            ):

                gene = futures[future]

                try:

                    row = future.result(
                        timeout=3600
                    )

                    if row is not None:
                        summary_rows.append(row)

                except Exception as e:

                    print(
                        f"Error encountered while phasing "
                        f"{gene}: {e}"
                    )

                    traceback.print_exc()

    else:

        row = phase_reads(
            args.bam,
            args.name,
            args.region,
            args.nanoTS_vcf,
            args.clair3_vcf,
            args.deepvariant_vcf,
            args.min_dp,
            args.min_af,
            args.min_distance_from_read_end,
            args.terminal_variant_proportion,
            ignore_variants_bed,
            snvs_only,
            args.phasing_threshold,
            args.name,
            args.genome,
            args.outdir,
            args.threads,
            args.remove_monoexonic
        )

        if row is not None:
            summary_rows.append(row)

    ##################################################################
    # Combined haplotype table
    ##################################################################

    ase_infile = os.path.join(
        args.outdir,
        f"{args.name}_phasing_summary.tsv"
    )

    os.makedirs(
        os.path.dirname(ase_infile) or ".",
        exist_ok=True
    )

    with open(
        ase_infile,
        "w"
    ) as f:

        f.write(
            "sample\tgene\thap1_read_count\thap2_read_count\t"
            "unassigned_read_count\tmax_coverage\tmax_phased_coverage\n"
        )

        for (
            gene,
            hap1_count,
            hap2_count,
            unassigned_count,
            max_coverage,
            max_phased_coverage
        ) in sorted(
            summary_rows,
            key=lambda r: r[0]
        ):

            f.write(
                f"{args.name}\t{gene}\t"
                f"{hap1_count}\t{hap2_count}\t"
                f"{unassigned_count}\t"
                f"{max_coverage}\t"
                f"{max_phased_coverage}\n"
            )

    print(
        f"\nWrote combined ASE input table for "
        f"{len(summary_rows)} gene(s) to {ase_infile}"
    )

    ##################################################################
    # Gene/BAM mapping file
    ##################################################################

    if args.bed:

        mapping_file = os.path.join(
            args.outdir,
            f"{args.name}_gene_bam_mapping_file.tsv"
        )

        with open(
            mapping_file,
            "w"
        ) as f:

            f.write(
                "name\tregion\tgene\tbulk_bam\thap1_bam\thap2_bam\n"
            )

            for gene, (chrom, start, end) in sorted(
                gene_regions.items()
            ):

                gene_outdir = os.path.join(
                    args.outdir,
                    gene
                )

                bulk_bam = os.path.join(
                    gene_outdir,
                    f"{args.name}_{gene}.bam"
                )

                hap1_bam = os.path.join(
                    gene_outdir,
                    f"{args.name}_{gene}_hap1.bam"
                )

                hap2_bam = os.path.join(
                    gene_outdir,
                    f"{args.name}_{gene}_hap2.bam"
                )

                region = f"{chrom}:{start}-{end}"

                if (
                    os.path.exists(bulk_bam)
                    and
                    os.path.exists(hap1_bam)
                    and
                    os.path.exists(hap2_bam)
                ):

                    f.write(
                        f"{args.name}\t{region}\t{gene}\t"
                        f"{bulk_bam}\t{hap1_bam}\t{hap2_bam}\n"
                    )

                elif os.path.exists(bulk_bam):

                    f.write(
                        f"{args.name}\t{region}\t{gene}\t"
                        f"{bulk_bam}\t\t\n"
                    )

        print(
            f"Wrote gene/BAM mapping file for "
            f"{len(gene_regions)} gene(s) to {mapping_file}"
        )

    if args.ignore_variants_list:
        os.remove(ignore_variants_bed)

    print(
        "\nFinished phasing reads.\n"
    )


if __name__ == "__main__":
    main()