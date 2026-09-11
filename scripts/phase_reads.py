#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2025.01.24
# Optimized: 2025
# Revised: 2026 -- read phasing now trusts NanoTS's genotypes AND its own
# phase directly, and selects its own per-gene variant set straight from
# --nanoTS-vcf (no more build_vcf_for_phasing.py / separate merged-VCF-
# building rule). Per gene, select_nanoTS_variants_for_gene() takes every
# PASS/DP-filtered, heterozygous NanoTS variant overlapping the gene
# region, PLUS every PASS/DP-filtered heterozygous NanoTS variant anywhere
# else in the VCF sharing a phase set (PS) with one of those -- NanoTS
# phases its own calls internally (its own output is literally named
# phased_predict.pass.vcf), so a shared PS means NanoTS already established
# two variants are on the same haplotype block even when one sits outside
# the gene's own boundaries. No AF filter is applied to NanoTS (FILTER=PASS
# and DP only) -- its genotype is trusted directly.
#
# Whatever VCF ends up used for haplotag is written to (and persists at)
# {outdir}/phased_reads/{gene}/{gene}_phased.vcf.gz -- not a scratch/temp
# file -- so it's always available to inspect after the fact, whichever of
# the cases below actually produced it.
#
# Whatshap is used ONLY to haplotag reads (assign each read to whichever
# haplotype its alleles match), NEVER to phase (`whatshap phase` is not run
# anywhere in this script). Per gene, phase_reads() first computes BAM-
# recomputed coverage (samtools depth over the filtered BAM, not each
# source's self-reported DP) over the WHOLE pool of selected NanoTS
# variants and candidate (Clair3/DeepVariant-shared) indels, tracking
# whichever single site has the highest coverage -- unconditionally, not
# just as a fallback, since it's also used as a sanity check below. Then:
#   - If 2+ of the gene's selected NanoTS variants share a real PS (a
#     genuine NanoTS-derived phased block), that block's variants are
#     written with their exact NanoTS phase (allele order + PS) preserved
#     -- write_nanoTS_vcf() -- and haplotagged directly (run_haplotag()).
#     If fewer reads end up phased this way than the single max-coverage
#     variant's own coverage, the block result is discarded and haplotag
#     is redone using just that one max-coverage variant instead -- a
#     multi-site block isn't assumed to always beat a single well-covered
#     site.
#   - Otherwise (no shared-PS block), haplotag runs directly against the
#     single max-coverage variant from the pool (NanoTS or candidate indel,
#     whichever it is). A candidate indel is only ever used this way, as
#     one option in a coverage-ranked pool -- never assigned to or tested
#     against an already-established NanoTS haplotype. If the pool is
#     empty, the gene simply isn't phased.

import argparse
import os
import subprocess
import pysam
import numpy as np
import shutil
import math
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
             "call -- see select_candidate_indels_for_gene(). Whatshap never sees these directly; they're "
             "only used as fallback single-variant haplotag anchors when a gene has no NanoTS phased "
             "block -- see phase_reads().")
    parser.add_argument("--deepvariant-vcf", type=str, required=True,
        help="DeepVariant VCF, used the same way as --clair3-vcf (indel agreement partner).")
    parser.add_argument("--region", type=str)
    parser.add_argument("--bed", type=str)
    parser.add_argument("--genome", type=str, required=True)
    parser.add_argument("--outdir", type=str, required=True)
    parser.add_argument("--name", type=str, default="SAMPLE")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--snvs-only', action='store_true',
        help="Exclude NanoTS indels from the trusted set (SNVs only). Candidate indels (Clair3/"
             "DeepVariant-shared) are unaffected by this flag.")
    parser.add_argument('--min-dp', type=int, default=20,
        help="Minimum DP for a NanoTS or candidate-indel variant to be selected at all. Default: 20")
    parser.add_argument('--min-af', type=float, default=0.1,
        help="Minimum AF for a candidate indel (Clair3/DeepVariant-shared) to be selected. NOT applied "
             "to NanoTS variants, which are trusted directly regardless of AF. Default: 0.1")
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

def filter_bam_by_region(inbam, outbam, region, threads):
    """Filter a BAM file for primary and supplementary alignments over a region."""
    subprocess.run([samtools_exec, 'view', '-F', '256', '-@', str(threads), '-b', '-o', outbam, inbam, region],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([samtools_exec, 'index', '-@', str(threads), outbam],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


def remove_monoexonic_reads(inbam, outbam, threads=1):
    """Remove reads lacking an 'N' CIGAR operation and write spliced reads to a new BAM."""
    with pysam.AlignmentFile(inbam, "rb") as infile, \
         pysam.AlignmentFile(outbam, "wb", header=infile.header, threads=threads) as outfile:
        for read in infile.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            # cigartuples op 3 == N (intron)
            if read.cigartuples and any(op == 3 for op, _ in read.cigartuples):
                outfile.write(read)

    sorted_bam = outbam.replace(".bam", ".sorted.bam")
    subprocess.run([samtools_exec, "sort", "-@", str(threads), "-o", sorted_bam, outbam], check=True)
    os.replace(sorted_bam, outbam)
    subprocess.run([samtools_exec, "index", "-@", str(threads), outbam], check=True)
    return outbam


def _load_ignore_positions(ignore_variants_bed):
    """Returns a set of (chrom, pos) 1-based positions to exclude from
    variant selection, or None if no bed was given."""
    if not ignore_variants_bed:
        return None
    positions = set()
    with open(ignore_variants_bed) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            for p in range(start + 1, end + 1):  # BED 0-based half-open -> 1-based positions
                positions.add((chrom, p))
    return positions


def write_variants_vcf(outfile, variants, sample_name, contigs):
    """Writes a minimal plain-text VCF for candidate indels
    ({(chrom,pos,ref,alt): (GT, DP)}), bgzip+tabix'd in place. GT is always
    written with '/' (unphased) -- candidate indels' own reported genotype
    is irrelevant downstream: they're only ever used for BAM-recomputed-
    coverage ranking as one option in phase_reads()'s single-variant
    fallback pool, which doesn't look at this GT field at all. See
    write_nanoTS_vcf() for the trusted NanoTS set, which -- unlike this --
    preserves real phase."""
    if outfile.endswith(".gz"):
        outfile = outfile[:-3]
    sorted_variants = sorted(variants.items())
    with open(outfile, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in contigs:
            f.write(f"##contig=<ID={c}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for (chrom, pos, ref, alt), (GT, DP) in sorted_variants:
            gt_str = "/".join(map(str, GT)) if GT else "./."
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP\t{gt_str}:{DP}\n")
    subprocess.run([bcftools_exec, "sort", "-Oz", "-o", f"{outfile}.gz", outfile],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    os.remove(outfile)
    subprocess.run([tabix_exec, "-f", "-p", "vcf", f"{outfile}.gz"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


def write_nanoTS_vcf(outfile, variants, sample_name, contigs):
    """Writes a minimal plain-text VCF for the trusted NanoTS set
    ({(chrom,pos,ref,alt): (GT, phased, PS, DP)}), bgzip+tabix'd in place.
    Unlike write_variants_vcf(), this PRESERVES NanoTS's own phase exactly
    as reported -- allele order ('|' vs '/') and PS -- since whatshap is
    only ever used downstream to haplotag reads against this VCF, never to
    re-derive phase (see this script's module docstring). PS is always
    written (as "." when NanoTS didn't report one) so every record has a
    consistent FORMAT column."""
    if outfile.endswith(".gz"):
        outfile = outfile[:-3]
    sorted_variants = sorted(variants.items())
    with open(outfile, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in contigs:
            f.write(f"##contig=<ID={c}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        f.write('##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for (chrom, pos, ref, alt), (GT, phased, PS, DP) in sorted_variants:
            if not GT:
                gt_str = "./."
            else:
                sep = "|" if phased else "/"
                gt_str = sep.join(map(str, GT))
            ps_str = str(PS) if PS not in (None, "") else "."
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP:PS\t{gt_str}:{DP}:{ps_str}\n")
    subprocess.run([bcftools_exec, "sort", "-Oz", "-o", f"{outfile}.gz", outfile],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    os.remove(outfile)
    subprocess.run([tabix_exec, "-f", "-p", "vcf", f"{outfile}.gz"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


def read_nanoTS_vcf(vcf_path):
    """Reverse of write_nanoTS_vcf(): reads a VCF written by it (or passed
    through remove_end_variants(), which preserves FORMAT fields via the
    same header) back into {(chrom,pos,ref,alt): (GT, phased, PS, DP)}.
    Used after remove_end_variants() to resume working with the filtered
    variant set as a plain dict again."""
    variants = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.fetch():
            if not rec.alts:
                continue
            s = rec.samples[0]
            PS = s.get("PS")
            if PS in (".", "", 0, "0"):
                PS = None
            variants[(rec.chrom, rec.pos, rec.ref, rec.alts[0])] = (
                s.get("GT"), s.phased, PS, s.get("DP") or 0)
    return variants


def select_nanoTS_variants_for_gene(nanoTS_vcf_path, chrom, start, end, min_dp,
                                     snvs_only=False, ignore_positions=None):
    """Selects the trusted, heterozygous NanoTS variant set for one gene:
    every PASS/DP-filtered heterozygous variant overlapping the gene
    region, PLUS every PASS/DP-filtered heterozygous variant anywhere else
    in the VCF that shares a phase set (PS) with one of those. NanoTS
    phases its own calls internally (its output is literally named
    phased_predict.pass.vcf -- see rules/1_call_variants.smk), so a shared
    PS means NanoTS already established two variants are on the same
    haplotype block, even when one of them sits outside the gene's own
    boundaries. No AF filter -- NanoTS's genotypes are trusted directly
    here, not treated as a confidence-scored candidate the way Clair3/
    DeepVariant indels still are in select_candidate_indels_for_gene().

    NanoTS's own phase (GT allele order + PS) is preserved exactly as
    reported, not discarded -- see this script's module docstring for why
    whatshap is only ever used to haplotag reads against this, never to
    re-derive phase itself. Returns {(chrom,pos,ref,alt): (GT, phased, PS, DP)}."""
    def _passes(rec):
        if rec.filter.keys() != ["PASS"]:
            return None
        if len(rec.alts) != 1:
            return None
        ref, alt = rec.ref, rec.alts[0]
        if snvs_only and len(ref) == 1 and len(alt) == 1:
            pass  # SNV, allowed
        elif snvs_only:
            return None  # indel, excluded when snvs_only
        if ignore_positions and (rec.chrom, rec.pos) in ignore_positions:
            return None
        s = rec.samples[0]
        GT = s.get("GT")
        if not GT or None in GT or len(set(GT)) < 2:
            return None  # missing or homozygous -- not heterozygous
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
            key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            PS = s.get("PS")
            variants[key] = (s.get("GT"), s.phased, PS, s.get("DP") or 0)
            if PS not in (None, ".", 0, "0"):
                phase_sets.add(PS)

        # Second pass, unrestricted by region: NanoTS is a targeted-panel
        # caller, so its whole per-sample VCF is small -- a full scan here
        # (once per gene, in parallel with other genes) is cheap, and is
        # the only way to correctly pull in same-PS variants that sit
        # outside this gene's own region.
        if phase_sets:
            for rec in vcf.fetch():
                if not rec.alts:
                    continue
                key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
                if key in variants:
                    continue
                s = _passes(rec)
                if s is None:
                    continue
                if s.get("PS") in phase_sets:
                    variants[key] = (s.get("GT"), s.phased, s.get("PS"), s.get("DP") or 0)

    return variants


def select_candidate_indels_for_gene(clair3_vcf_path, deepvariant_vcf_path, chrom, start, end,
                                      min_dp, min_af, exclude_keys, ignore_positions=None):
    """Selects candidate indels for one gene: biallelic PASS indels called
    by BOTH Clair3 and DeepVariant (PASS/DP/AF-filtered independently in
    each), restricted to this gene's region, excluding anything already in
    `exclude_keys` (the gene's trusted NanoTS set from
    select_nanoTS_variants_for_gene() -- those don't need separate post-hoc
    phasing, they're already trusted and already whatshap-phased). Returns
    {(chrom,pos,ref,alt): (GT, DP)}."""
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
                    continue  # indels only
                if ignore_positions and (rec.chrom, rec.pos) in ignore_positions:
                    continue
                s = rec.samples[0]
                DP = s.get("DP") or 0
                AF = s.get("AF") or s.get("VAF") or (0,)
                if isinstance(AF, (list, tuple)):
                    AF = AF[0]
                if DP < min_dp or AF < min_af:
                    continue
                out[(rec.chrom, rec.pos, ref, alt)] = (s.get("GT"), DP)
        return out

    clair3 = _get(clair3_vcf_path)
    deepvariant = _get(deepvariant_vcf_path)
    shared = {k: clair3[k] for k in clair3.keys() & deepvariant.keys()}
    return {k: v for k, v in shared.items() if k not in exclude_keys}


def remove_end_variants(invcf, outvcf, bam, min_distance_from_read_end=10, terminal_variant_proportion=0.5):
    """Remove variants that frequently occur close to the end of an alignment."""
    if outvcf.endswith('.gz'):
        outvcf = outvcf[:-3]

    bamfile = pysam.AlignmentFile(bam, "rb")
    vcf_reader = pysam.VariantFile(invcf, "r")
    vcf_writer = pysam.VariantFile(outvcf, "w", header=vcf_reader.header)

    for record in vcf_reader:
        chrom = record.chrom
        pos = record.pos - 1  # VCF 1-based → BAM 0-based
        ref = record.ref
        alts = record.alts

        for alt in alts:
            total_alt_alignments, end_alt_alignments = 0, 0
            for alignment in bamfile.fetch(chrom, pos, pos + len(alt)):
                # Use get_aligned_pairs once and cache the result
                aligned_pairs = alignment.get_aligned_pairs(matches_only=True)
                query_pos = None
                for qpos, rpos in aligned_pairs:
                    if rpos == pos:
                        query_pos = qpos
                        break
                if query_pos is None:
                    continue

                query_seq = alignment.query_sequence[query_pos: query_pos + len(alt)]
                if query_seq == alt:
                    total_alt_alignments += 1
                    # Compute bounds from the same aligned_pairs list (already fetched above)
                    q_positions = [qpos for qpos, _ in aligned_pairs]
                    q_start = q_positions[0]
                    q_end = q_positions[-1]
                    distance_to_end = min(query_pos - q_start,
                                          q_end - (query_pos + max(len(ref), len(alt))))
                    if distance_to_end < min_distance_from_read_end:
                        end_alt_alignments += 1

        if total_alt_alignments == 0 or (end_alt_alignments / total_alt_alignments) > terminal_variant_proportion:
            continue
        vcf_writer.write(record)

    bamfile.close()
    vcf_reader.close()
    vcf_writer.close()

    subprocess.run([bcftools_exec, 'sort', '-Oz', '-o', f"{outvcf}.gz", outvcf],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([tabix_exec, '-f', '-p', 'vcf', f"{outvcf}.gz"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


def annotate_read_depth(invcf, bam, outvcf, tempdir):
    """Re-annotates a VCF's DP with BAM-recomputed coverage (samtools depth
    over `bam`), overwriting whatever DP the source caller(s) reported.
    Despite the old name, `invcf` is no longer necessarily whatshap output
    -- it's just any VCF (e.g. the trusted NanoTS set, or the gene-scoped
    candidate indels) needing coverage recomputed on equal footing with
    another VCF being compared against it."""
    with pysam.VariantFile(invcf) as v, \
         open(os.path.join(tempdir, 'variants.bed'), 'w') as bedfile:
        for record in v.fetch():
            bedfile.write(f'{record.chrom}\t{record.pos - 1}\t{record.pos}\n')

    with open(os.path.join(tempdir, 'variant_coverage.tsv'), 'w') as f:
        subprocess.run([samtools_exec, 'depth', '-b', os.path.join(tempdir, 'variants.bed'), bam],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([bgzip_exec, '-f', os.path.join(tempdir, 'variant_coverage.tsv')],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([tabix_exec, '-f', '-s', '1', '-b', '2', '-e', '2',
                    os.path.join(tempdir, 'variant_coverage.tsv.gz')],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    with open(os.path.join(tempdir, 'header.txt'), 'w') as f:
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    subprocess.run(
        [bcftools_exec, 'annotate',
         '-a', os.path.join(tempdir, 'variant_coverage.tsv.gz'),
         '-h', os.path.join(tempdir, 'header.txt'),
         '-c', 'CHROM,POS,FORMAT/DP', invcf, '-Oz', '-o', outvcf],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([tabix_exec, '-f', '-p', 'vcf', outvcf],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


def get_phased_coverage(bam_path, hap1_bam_path, hap2_bam_path, region):
    """Get the max coverage and max phased coverage safely."""
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))

    def max_cov(bam_file):
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            A, C, G, T = bam.count_coverage(chrom, start, end, quality_threshold=0)
            # np.add is faster than element-wise zip sum
            return int(np.max(np.add(np.add(A, C), np.add(G, T))))

    total_max = max_cov(bam_path)
    phased_max = max_cov(hap1_bam_path) + max_cov(hap2_bam_path)
    return total_max, phased_max


def run_haplotag(phased_vcf, filtered_bam, genome, gene_outdir):
    """Runs `whatshap haplotag` against `phased_vcf` and returns
    (hap1_reads, hap2_reads, unassigned_reads) as lists of read names --
    just the read-name-to-haplotype assignment, no BAM writing. Cheap
    enough to call more than once per gene (phase_reads() does, when a
    NanoTS phased block's own haplotag result underperforms the single
    max-coverage-variant fallback and needs redoing with that instead) --
    unlike whatshap split, which actually writes and indexes BAM files and
    is only worth doing once the final phased_vcf choice is settled. No
    `whatshap phase` is run anywhere in this pipeline; NanoTS already
    phases its own calls internally (its own output is literally named
    phased_predict.pass.vcf -- see rules/1_call_variants.smk), so whichever
    phase phase_reads() decided on is what haplotag uses directly."""
    haplotag_outfile = os.path.join(gene_outdir, 'whatshap_haplotag.tsv')
    subprocess.run(
        [whatshap_exec, "haplotag", "--reference", genome, "--output", "/dev/null",
         "--ignore-read-groups", "--output-haplotag-list", haplotag_outfile, phased_vcf, filtered_bam],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    hap1_reads, hap2_reads, unassigned_reads = [], [], []
    with open(haplotag_outfile) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            hap = cols[1]
            if hap == "H1":
                hap1_reads.append(cols[0])
            elif hap == "H2":
                hap2_reads.append(cols[0])
            else:
                unassigned_reads.append(cols[0])

    return hap1_reads, hap2_reads, unassigned_reads


def finalize_haplotype_bams(filtered_bam, gene_outdir, name, gene, region, phasing_threshold, threads,
                             remove_monoexonic, phased_vcf, hap1_reads, hap2_reads, unassigned_reads):
    """Given the final (hap1_reads, hap2_reads, unassigned_reads) from
    run_haplotag() against the settled `phased_vcf`, writes/indexes hap1/
    hap2/unassigned BAMs (`whatshap split`), applies the phasing_threshold
    gate, and writes the per-gene summary."""
    report_message = ""

    haplotag_outfile = os.path.join(gene_outdir, 'whatshap_haplotag.tsv')
    hap1_bam = os.path.join(gene_outdir, f'{name}_{gene}_hap1.bam')
    hap2_bam = os.path.join(gene_outdir, f'{name}_{gene}_hap2.bam')
    unassigned_bam = os.path.join(gene_outdir, f'{name}_{gene}_unassigned.bam')
    subprocess.run(
        [whatshap_exec, "split", "--output-h1", hap1_bam, "--output-h2", hap2_bam,
         "--output-untagged", unassigned_bam, filtered_bam, haplotag_outfile],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    for bam_path in (hap1_bam, hap2_bam, unassigned_bam):
        subprocess.run([samtools_exec, 'index', '-@', str(threads), bam_path],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    num_phased_reads = len(hap1_reads) + len(hap2_reads)
    max_coverage, max_phased_coverage = get_phased_coverage(filtered_bam, hap1_bam, hap2_bam, region)
    report_message += (
        f"\n    {num_phased_reads} reads can be assigned to haplotypes."
        f"\n      {len(unassigned_reads)} reads can not be assigned to a haplotype."
        f"\n    Max phased coverage: {max_phased_coverage}"
        f"\n    Max overall coverage: {max_coverage}"
    )

    if max_coverage == 0 or max_phased_coverage < phasing_threshold * max_coverage:
        for fp in [hap1_bam, hap1_bam + ".bai", hap2_bam, hap2_bam + ".bai",
                   unassigned_bam, unassigned_bam + ".bai"]:
            if os.path.exists(fp):
                os.remove(fp)
        report_message += (f"\nCannot phase {phasing_threshold*100}% of the max read coverage "
                           f"over {gene} ({region}). Exiting without phasing {name}...")
        return report_message, None

    report_message += (f"\nSuccessfully phased {num_phased_reads} out of "
                       f"{num_phased_reads + len(unassigned_reads)} reads over {gene} ({region}).")

    if remove_monoexonic:
        multidir = os.path.join(gene_outdir, "multiexonic_bams")
        for src_bam in (hap1_bam, hap2_bam):
            dest = os.path.join(multidir,
                                os.path.splitext(os.path.basename(src_bam))[0] + "_multiexonic.bam")
            remove_monoexonic_reads(src_bam, dest)

    summary_row = (gene, len(hap1_reads), len(hap2_reads), len(unassigned_reads),
                   max_coverage, max_phased_coverage)

    # Per-gene summary file is kept for human-readable debugging, but is no longer
    # read back in by any downstream rule — build_ase_infile has been retired in
    # favor of phase_reads.py always writing the combined table directly to
    # {outdir}/{name}_phasing_summary.tsv (see main()).
    with open(os.path.join(gene_outdir, f"{name}_{gene}_haplotype_assignment_summary.txt"), "w") as f:
        f.write("gene\thap1_read_count\thap2_read_count\tunassigned_read_count\tmax_coverage\tmax_phased_coverage\n")
        f.write(f"{gene}\t{len(hap1_reads)}\t{len(hap2_reads)}\t{len(unassigned_reads)}\t{max_coverage}\t{max_phased_coverage}\n")

    return report_message, summary_row


########################################################################################################################
# Main function for phasing reads
########################################################################################################################

def phase_reads(bam, gene, region, nanoTS_vcf, clair3_vcf, deepvariant_vcf, min_dp, min_af,
                min_distance_from_read_end, terminal_variant_proportion,
                ignore_variants_bed, snvs_only, phasing_threshold, name, genome, outdir, threads, remove_monoexonic):

    gene_outdir = os.path.join(outdir, gene)
    tempdir = os.path.join(gene_outdir, 'temp')
    os.makedirs(tempdir, exist_ok=True)
    if remove_monoexonic:
        os.makedirs(os.path.join(gene_outdir, 'multiexonic_bams'), exist_ok=True)

    chrom, positions = region.split(":")
    region_start, region_end = map(int, positions.split("-"))
    ignore_positions = _load_ignore_positions(ignore_variants_bed)

    report = os.path.join(gene_outdir, f'{name}_{gene}_phasing_report.txt')
    with open(report, 'w') as report_file:

        report_file.write(f"\nFiltering BAM file for reads aligned to {gene} ({region})...\n")
        report_file.flush()
        filtered_bam = os.path.join(gene_outdir, f'{name}_{gene}.bam')
        filter_bam_by_region(bam, filtered_bam, region, threads)
        if remove_monoexonic:
            multidir = os.path.join(gene_outdir, 'multiexonic_bams')
            remove_monoexonic_reads(
                filtered_bam,
                os.path.join(multidir, os.path.splitext(os.path.basename(filtered_bam))[0] + "_multiexonic.bam"))

        total_read_count = int(
            subprocess.run([samtools_exec, 'view', '-@', str(threads), '-c', filtered_bam],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode())
        if total_read_count == 0:
            report_file.write(f"    No reads aligned to {gene} ({region}). Exiting without phasing {name}...\n")
            report_file.write(f"\nFinished processing {name} over {gene} ({region}).")
            report_file.flush()
            os.remove(filtered_bam)
            os.remove(f"{filtered_bam}.bai")
            shutil.rmtree(tempdir)
            return None
        report_file.write(f"    {total_read_count} primary/supplementary alignments found.")
        report_file.flush()

        report_file.write(f"\nSelecting trusted, heterozygous NanoTS variants for {gene} ({region})...")
        report_file.flush()
        nanoTS_variants = select_nanoTS_variants_for_gene(
            nanoTS_vcf, chrom, region_start, region_end, min_dp, snvs_only, ignore_positions)
        report_file.write(f"\n    {len(nanoTS_variants)} trusted NanoTS variant(s) selected (overlapping "
                          f"the gene region, or sharing a phase set with one that does).")
        report_file.flush()

        report_file.write(f"\nSelecting candidate (Clair3/DeepVariant-shared) indels for {gene} ({region})...")
        report_file.flush()
        candidate_indels = select_candidate_indels_for_gene(
            clair3_vcf, deepvariant_vcf, chrom, region_start, region_end, min_dp, min_af,
            set(nanoTS_variants.keys()), ignore_positions)
        report_file.write(f"\n    {len(candidate_indels)} candidate indel(s) selected (not already in the "
                          f"trusted NanoTS set).")
        report_file.flush()

        if nanoTS_variants:
            gene_vcf = os.path.join(tempdir, f'{gene}.vcf.gz')
            write_nanoTS_vcf(gene_vcf, nanoTS_variants, name, [chrom])

            if min_distance_from_read_end > 0 and terminal_variant_proportion < 1:
                report_file.write(
                    f"\nRemoving variants that occur within {min_distance_from_read_end}nt of the "
                    f"end of a read more than {terminal_variant_proportion*100}% of the time...")
                report_file.flush()
                filtered_gene_vcf = os.path.join(tempdir, f'filtered_{gene}.vcf.gz')
                remove_end_variants(gene_vcf, filtered_gene_vcf, filtered_bam,
                                     min_distance_from_read_end, terminal_variant_proportion)
                nanoTS_variants = read_nanoTS_vcf(filtered_gene_vcf)
                report_file.write(f"\n    {len(nanoTS_variants)} variant(s) remain.")
                report_file.flush()

        # Coverage (BAM-recomputed via samtools depth, not each source's
        # self-reported DP) is computed over the WHOLE pool -- every
        # surviving NanoTS variant plus every candidate indel -- regardless
        # of whether a NanoTS phased block exists, since it's needed either
        # way: as the sole basis for haplotagging when there's no block, or
        # as a sanity check against the block's own haplotag result when
        # there is one (see below).
        pool_keys = list(nanoTS_variants.keys()) + list(candidate_indels.keys())
        best_key, best_dp = None, -1
        if pool_keys:
            report_file.write(f"\nComputing BAM-recomputed coverage over all {len(pool_keys)} "
                              f"NanoTS variant(s) and candidate indel(s) for {gene}...")
            report_file.flush()
            pool_vcf = os.path.join(tempdir, f'{gene}_pool.vcf.gz')
            write_variants_vcf(pool_vcf, {k: (None, 0) for k in pool_keys}, name, [chrom])
            pool_annotated = os.path.join(tempdir, f'{gene}_pool_annotated.vcf.gz')
            annotate_read_depth(pool_vcf, filtered_bam, pool_annotated, tempdir)
            with pysam.VariantFile(pool_annotated) as vcf:
                for rec in vcf.fetch():
                    if not rec.alts:
                        continue
                    dp = rec.samples[0].get("DP") or 0
                    if dp > best_dp:
                        best_dp = dp
                        best_key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            chrom_b, pos_b, ref_b, alt_b = best_key
            source = "NanoTS variant" if best_key in nanoTS_variants else "candidate indel"
            report_file.write(f"\n    Highest-coverage variant: a {source} at "
                              f"{chrom_b}:{pos_b} {ref_b}>{alt_b} (DP={best_dp}).")
            report_file.flush()

        def _single_variant_vcf():
            out = os.path.join(gene_outdir, f'{gene}_phased.vcf.gz')
            chrom_b, pos_b, ref_b, alt_b = best_key
            write_nanoTS_vcf(out, {best_key: ((0, 1), True, str(pos_b), best_dp)}, name, [chrom])
            return out

        # A phased block exists if 2+ of the (surviving) NanoTS variants
        # share a real PS -- that's NanoTS's own phasing, established
        # independently of this BAM.
        ps_counts = {}
        for _, phased_flag, PS, _ in nanoTS_variants.values():
            if PS is not None:
                ps_counts[PS] = ps_counts.get(PS, 0) + 1
        phased_block_ps = {ps for ps, n in ps_counts.items() if n >= 2}

        if phased_block_ps:
            block_variants = {k: v for k, v in nanoTS_variants.items() if v[2] in phased_block_ps}
            phased_vcf = os.path.join(gene_outdir, f'{gene}_phased.vcf.gz')
            write_nanoTS_vcf(phased_vcf, block_variants, name, [chrom])
            report_file.write(f"\nUsing NanoTS's own phased block ({len(block_variants)} variant(s), "
                              f"{len(phased_block_ps)} phase set(s)) for whatshap haplotag...")
            report_file.flush()
            hap1_reads, hap2_reads, unassigned_reads = run_haplotag(phased_vcf, filtered_bam, genome, gene_outdir)
            n_block_phased = len(hap1_reads) + len(hap2_reads)
            report_file.write(f"\n    {n_block_phased} read(s) phased from the block.")
            report_file.flush()

            if best_key is not None and n_block_phased < best_dp:
                report_file.write(f"\n    Fewer reads phased ({n_block_phased}) than the max single-variant "
                                  f"coverage ({best_dp}); falling back to that single variant and redoing "
                                  f"haplotag...")
                report_file.flush()
                phased_vcf = _single_variant_vcf()
                hap1_reads, hap2_reads, unassigned_reads = run_haplotag(phased_vcf, filtered_bam, genome, gene_outdir)
        elif best_key is not None:
            phased_vcf = _single_variant_vcf()
            report_file.write(f"\nNo NanoTS phased block found; phasing reads based on the single "
                              f"highest-(BAM-recomputed)-coverage variant above for whatshap haplotag...")
            report_file.flush()
            hap1_reads, hap2_reads, unassigned_reads = run_haplotag(phased_vcf, filtered_bam, genome, gene_outdir)
        else:
            # Nothing at all to anchor on -- there's no variant to build a
            # haplotype split from, so skip haplotagging entirely rather
            # than running it against an empty VCF (which would always
            # phase 0 reads, and -- since 0 coverage trivially satisfies
            # the phasing_threshold ratio check in finalize_haplotype_bams,
            # 0 < threshold * 0 is False -- could get misreported as a
            # "successful" 0-read phasing instead of being cleanly skipped).
            # filtered_bam (bulk) is kept, same as any other gene that fails
            # to phase -- only hap1/hap2/unassigned BAMs are ever skipped.
            report_file.write(f"\nNo NanoTS variants or candidate indels found for {gene} ({region}); "
                              f"nothing to phase reads on. Exiting without attempting to haplotag...")
            report_file.write(f"\nFinished processing {name} over {gene} ({region}).")
            report_file.flush()
            shutil.rmtree(tempdir)
            return None

        report_file.write(f"\nCreating haplotype-specific BAM files...")
        report_file.flush()
        report_message, summary_row = finalize_haplotype_bams(
            filtered_bam, gene_outdir, name, gene, region, phasing_threshold, threads, remove_monoexonic,
            phased_vcf, hap1_reads, hap2_reads, unassigned_reads)
        report_file.write(report_message)

        report_file.write(f"\nFinished processing {name} over {gene} ({region}).")
        report_file.flush()

    shutil.rmtree(tempdir)
    return summary_row


########################################################################################################################
# Main script
########################################################################################################################

def main():
    """Main function."""

    print(f"\n\n\n******************************************************************************************")
    print(f"Phasing reads from RNA-seq data...")
    print(f"******************************************************************************************\n")

    args = parse_args()
    print(f"Preparing to phase reads for {args.name}...\n")

    global samtools_exec, bcftools_exec, tabix_exec, bgzip_exec, minimap2_exec, whatshap_exec
    samtools_exec = args.samtools_exec
    bcftools_exec = args.bcftools_exec
    tabix_exec = args.tabix_exec
    bgzip_exec = args.bgzip_exec
    minimap2_exec = args.minimap2_exec
    whatshap_exec = args.whatshap_exec

    if not os.path.exists(args.bam):
        raise FileNotFoundError(f"BAM file {args.bam} not found.")
    if not os.path.exists(args.genome):
        raise FileNotFoundError(f"Genome file {args.genome} not found.")
    if args.bed and not os.path.exists(args.bed):
        raise FileNotFoundError(f"BED file {args.bed} not found.")
    if not os.path.exists(args.nanoTS_vcf):
        raise FileNotFoundError(f"NanoTS VCF file {args.nanoTS_vcf} not found.")
    if not os.path.exists(args.clair3_vcf):
        raise FileNotFoundError(f"Clair3 VCF file {args.clair3_vcf} not found.")
    if not os.path.exists(args.deepvariant_vcf):
        raise FileNotFoundError(f"DeepVariant VCF file {args.deepvariant_vcf} not found.")
    if args.ignore_variants_bed and not os.path.exists(args.ignore_variants_bed):
        raise FileNotFoundError(f"BED file containing variants to ignore ({args.ignore_variants_bed}) not found.")
    if not args.region and not args.bed:
        raise ValueError("Either --region or --bed must be specified.")

    os.makedirs(args.outdir, exist_ok=True)

    if args.ignore_variants_list:
        ignore_variants_bed = os.path.join(args.outdir, 'ignore_variants.bed')
        with open(ignore_variants_bed, 'w') as f:
            for variant in args.ignore_variants_list:
                chrom, pos = variant.split(':')
                f.write(f"{chrom}\t{int(pos)-1}\t{pos}\n")
    elif args.ignore_variants_bed:
        ignore_variants_bed = args.ignore_variants_bed
    else:
        ignore_variants_bed = None

    snvs_only = bool(args.snvs_only)

    summary_rows = []

    if args.bed:
        gene_regions = extract_gene_regions(args.bed)
        threads = min(args.threads, len(gene_regions))
        threads_per_process = max(1, args.threads // len(gene_regions))
        print(f"\nBegin phasing reads in {len(gene_regions)} gene regions of interest in parallel "
              f"with {threads} workers. This may take a while...\n")
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(
                    phase_reads, args.bam, gene, f"{chrom}:{start}-{end}",
                    args.nanoTS_vcf, args.clair3_vcf, args.deepvariant_vcf, args.min_dp, args.min_af,
                    args.min_distance_from_read_end, args.terminal_variant_proportion,
                    ignore_variants_bed, snvs_only, args.phasing_threshold, args.name,
                    args.genome, args.outdir, threads_per_process, args.remove_monoexonic): gene
                for gene, (chrom, start, end) in gene_regions.items()
            }
        for future in concurrent.futures.as_completed(futures):
            gene = futures[future]
            try:
                row = future.result(timeout=3600)
                if row is not None:
                    summary_rows.append(row)
            except Exception as e:
                print(f"Error encountered while phasing {gene}: {e}")
                traceback.print_exc()
    else:
        row = phase_reads(args.bam, args.name, args.region,
                    args.nanoTS_vcf, args.clair3_vcf, args.deepvariant_vcf, args.min_dp, args.min_af,
                    args.min_distance_from_read_end, args.terminal_variant_proportion,
                    ignore_variants_bed, snvs_only, args.phasing_threshold, args.name,
                    args.genome, args.outdir, args.threads, args.remove_monoexonic)
        if row is not None:
            summary_rows.append(row)

    # Write the combined per-gene haplotype table directly — this is what
    # detect_ase_outliers.py consumes. Writing it here (rather than reconstructing
    # it downstream via a shell/awk pass over per-gene summary files) guarantees the
    # sample name and column count are always correct and in sync with each other.
    # Always written to a fixed, predictable path — {outdir}/{name}_phasing_summary.tsv —
    # rather than taken as a CLI arg, so callers don't need to know or specify it.
    ase_infile = os.path.join(args.outdir, f"{args.name}_phasing_summary.tsv")
    os.makedirs(os.path.dirname(ase_infile) or ".", exist_ok=True)
    with open(ase_infile, "w") as f:
        f.write("sample\tgene\thap1_read_count\thap2_read_count\t"
                 "unassigned_read_count\tmax_coverage\tmax_phased_coverage\n")
        for gene, hap1_count, hap2_count, unassigned_count, max_coverage, max_phased_coverage in \
                sorted(summary_rows, key=lambda r: r[0]):
            f.write(f"{args.name}\t{gene}\t{hap1_count}\t{hap2_count}\t"
                    f"{unassigned_count}\t{max_coverage}\t{max_phased_coverage}\n")
    print(f"\nWrote combined ASE input table for {len(summary_rows)} gene(s) to {ase_infile}")

    # Index phase_reads' own per-gene BAM outputs (bulk + hap1/hap2, written
    # above by phase_reads()) for get_splice_junction_counts in
    # junction_analysis.smk. Only meaningful in --bed mode (gene_regions is
    # only built above when a BED file was given); a single --region run has
    # no gene to key rows on. Always written to a fixed, predictable path --
    # {outdir}/{name}_gene_bam_mapping_file.tsv -- for the same reason as
    # ase_infile above.
    if args.bed:
        mapping_file = os.path.join(args.outdir, f"{args.name}_gene_bam_mapping_file.tsv")
        with open(mapping_file, "w") as f:
            f.write("name\tregion\tgene\tbulk_bam\thap1_bam\thap2_bam\n")
            for gene, (chrom, start, end) in sorted(gene_regions.items()):
                gene_outdir = os.path.join(args.outdir, gene)
                bulk_bam = os.path.join(gene_outdir, f"{args.name}_{gene}.bam")
                hap1_bam = os.path.join(gene_outdir, f"{args.name}_{gene}_hap1.bam")
                hap2_bam = os.path.join(gene_outdir, f"{args.name}_{gene}_hap2.bam")
                region = f"{chrom}:{start}-{end}"
                if os.path.exists(bulk_bam) and os.path.exists(hap1_bam) and os.path.exists(hap2_bam):
                    f.write(f"{args.name}\t{region}\t{gene}\t{bulk_bam}\t{hap1_bam}\t{hap2_bam}\n")
                elif os.path.exists(bulk_bam):
                    f.write(f"{args.name}\t{region}\t{gene}\t{bulk_bam}\t\t\n")
        print(f"Wrote gene/BAM mapping file for {len(gene_regions)} gene(s) to {mapping_file}")

    if args.ignore_variants_list:
        os.remove(ignore_variants_bed)

    print(f"\nFinished phasing reads.\n")


if __name__ == "__main__":
    main()