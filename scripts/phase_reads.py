#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2025.01.24
# Optimized: 2025

import argparse
import os
import subprocess
import pysam
import numpy as np
import shutil
import math
import concurrent.futures

def parse_args():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Phases reads from an RNA-seq BAM file that cover gene regions of interest.")
    parser.add_argument("--bam", type=str, required=True)
    parser.add_argument("--vcf", type=str, required=True)
    parser.add_argument("--region", type=str)
    parser.add_argument("--bed", type=str)
    parser.add_argument("--genome", type=str, required=True)
    parser.add_argument("--outdir", type=str, required=True)
    parser.add_argument("--name", type=str, default="SAMPLE")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument('--snvs-only', action='store_true')
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


def filter_vcf(invcf, outvcf, region, snvs_only=False, ignore_variants_bed=None):
    """Filter a VCF file for bi-allelic heterozygous variants in a specific region."""
    if outvcf.endswith('.gz'):
        outvcf = outvcf[:-3]

    if ignore_variants_bed:
        tempvcf = os.path.join(os.path.dirname(outvcf), "temp.vcf")
        subprocess.run([bcftools_exec, 'view', '-T', f'^{ignore_variants_bed}', '-Oz', '-o', tempvcf, invcf],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        subprocess.run([bcftools_exec, 'sort', '-Oz', '-o', f"{tempvcf}.gz", tempvcf],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        subprocess.run([tabix_exec, '-f', '-p', 'vcf', f"{tempvcf}.gz"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        invcf = f"{tempvcf}.gz"

    if snvs_only:
        cmd = [bcftools_exec, 'view', '-m2', '-M2', '-v', 'snps', '--regions', region, '-Oz', '-o', outvcf, invcf]
    else:
        cmd = [bcftools_exec, 'view', '-m2', '-M2', '--regions', region, '-Oz', '-o', outvcf, invcf]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([bcftools_exec, 'sort', '-Oz', '-o', f"{outvcf}.gz", outvcf],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([tabix_exec, '-f', '-p', 'vcf', f"{outvcf}.gz"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    os.remove(outvcf)

    if ignore_variants_bed:
        os.remove(tempvcf)
        os.remove(f"{tempvcf}.gz")
        os.remove(f"{tempvcf}.gz.tbi")


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


def phase_variants(bam, genome, vcf, tempdir):
    """Phase variants in a BAM file using WhatsHap."""
    subprocess.run(
        [whatshap_exec, 'phase', '-o', os.path.join(tempdir, 'whatshap.vcf'),
         '--reference', genome, '--mapping-quality', '0', '--ignore-read-groups',
         '--distrust-genotypes', vcf, bam],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run(
        [bcftools_exec, 'sort', '-Oz', '-o', os.path.join(tempdir, 'whatshap.vcf.gz'),
         os.path.join(tempdir, 'whatshap.vcf')],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run(
        [tabix_exec, '-p', 'vcf', os.path.join(tempdir, 'whatshap.vcf.gz')],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)


def annotate_read_depth(whatshapvcf, bam, outvcf, tempdir):
    """Annotate a WhatsHap VCF file with read depth information."""
    with pysam.VariantFile(whatshapvcf) as v, \
         open(os.path.join(tempdir, 'variants.bed'), 'w') as bedfile:
        for record in v.fetch():
            bedfile.write(f'{record.chrom}\t{record.pos - 1}\t{record.pos}\n')

    with open(os.path.join(tempdir, 'variant_coverage.tsv'), 'w') as f:
        subprocess.run([samtools_exec, 'depth', '-b', os.path.join(tempdir, 'variants.bed'), bam],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([bgzip_exec, os.path.join(tempdir, 'variant_coverage.tsv')],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([tabix_exec, '-s', '1', '-b', '2', '-e', '2',
                    os.path.join(tempdir, 'variant_coverage.tsv.gz')],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    with open(os.path.join(tempdir, 'header.txt'), 'w') as f:
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    subprocess.run(
        [bcftools_exec, 'annotate',
         '-a', os.path.join(tempdir, 'variant_coverage.tsv.gz'),
         '-h', os.path.join(tempdir, 'header.txt'),
         '-c', 'CHROM,POS,FORMAT/DP', whatshapvcf, '-Oz', '-o', outvcf],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([tabix_exec, '-p', 'vcf', outvcf],
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


def create_haplotype_specific_bams(filtered_bam, genome, gene_outdir, tempdir, name, gene,
                                    region, phasing_threshold, threads, remove_monoexonic):
    """Haplotag reads using phased VCF."""
    report_message = ""

    phased_vcf = os.path.join(gene_outdir, 'phased.vcf.gz')
    gt_filter = 'GT="0|1" || GT="1|0"'
    subprocess.run([bcftools_exec, "view", "-i", gt_filter, "-Oz", "-o", phased_vcf,
                    os.path.join(gene_outdir, 'whatshap_annotated.vcf.gz')], check=True)
    subprocess.run([tabix_exec, '-f', '-p', 'vcf', phased_vcf],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    result = subprocess.run([bcftools_exec, 'view', '-H', phased_vcf],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode()
    has_phased = any(line.strip() for line in result.splitlines())

    if not has_phased:
        gt_filter = 'GT="0/1" || GT="1/0"'
        het_output = subprocess.run(
            [bcftools_exec, 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\n',
             '-i', gt_filter, os.path.join(gene_outdir, 'whatshap_annotated.vcf.gz')],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode()
        if het_output:
            best_line = max(het_output.splitlines(), key=lambda x: int(x.split('\t')[4]))
            chrom, pos, ref, alt, dp = best_line.split('\t')
            synthetic_vcf = os.path.join(gene_outdir, 'phased.vcf')
            header = subprocess.run([bcftools_exec, 'view', '-h', phased_vcf],
                                     stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode()
            with open(synthetic_vcf, "w") as f:
                f.write(header)
                f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:PS\t0|1:{pos}\n")
            phased_vcf = synthetic_vcf + ".gz"
            subprocess.run([bgzip_exec, '-f', synthetic_vcf],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            subprocess.run([tabix_exec, '-f', '-p', 'vcf', phased_vcf],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            report_message += (f"\n    Phasing reads based on a single heterozygous variant with "
                                f"the highest coverage at {chrom}:{pos} {ref}>{alt} (DP={dp})")

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

    if max_phased_coverage < phasing_threshold * max_coverage:
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

def phase_reads(bam, gene, region, vcf, min_distance_from_read_end, terminal_variant_proportion,
                ignore_variants_bed, snvs_only, phasing_threshold, name, genome, outdir, threads, remove_monoexonic):

    gene_outdir = os.path.join(outdir, gene)
    tempdir = os.path.join(gene_outdir, 'temp')
    os.makedirs(tempdir, exist_ok=True)
    if remove_monoexonic:
        os.makedirs(os.path.join(gene_outdir, 'multiexonic_bams'), exist_ok=True)

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

        report_file.write(f"\nFiltering VCF file for variants in {gene} ({region})...")
        report_file.flush()
        gene_vcf = os.path.join(tempdir, f'{gene}.vcf.gz')
        filter_vcf(vcf, gene_vcf, region, snvs_only, ignore_variants_bed)
        num_variants = int(
            subprocess.run([bcftools_exec, 'index', '-n', gene_vcf],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode())

        report_file.write(f"\n    {num_variants} variants found.")
        report_file.flush()
        if num_variants == 0:
            report_file.write(f"\n    No variants found in {gene} ({region}). Exiting without phasing {name}...\n")
            report_file.write(f"\nFinished processing {name} over {gene} ({region}).")
            report_file.flush()
            shutil.rmtree(tempdir)
            return None

        if min_distance_from_read_end > 0 and terminal_variant_proportion < 1:
            report_file.write(
                f"\nRemoving variants that occur within {min_distance_from_read_end}nt of the "
                f"end of a read more than {terminal_variant_proportion*100}% of the time...")
            report_file.flush()
            filtered_gene_vcf = os.path.join(tempdir, f'filtered_{gene}.vcf.gz')
            remove_end_variants(gene_vcf, filtered_gene_vcf, filtered_bam,
                                 min_distance_from_read_end, terminal_variant_proportion)
            num_variants = int(
                subprocess.run([bcftools_exec, 'index', '-n', filtered_gene_vcf],
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode())
            report_file.write(f"\n    {num_variants} variants remain.")
            report_file.flush()
            gene_vcf = filtered_gene_vcf

        report_file.write(f"\nPhasing variants with whatshap...")
        report_file.flush()
        phase_variants(filtered_bam, genome, gene_vcf, tempdir)

        report_file.write(f"\nAnnotating phased variants with read depth information...")
        report_file.flush()
        annotate_read_depth(os.path.join(tempdir, 'whatshap.vcf.gz'), filtered_bam,
                            os.path.join(gene_outdir, 'whatshap_annotated.vcf.gz'), tempdir)

        report_file.write(f"\nCreating haplotype-specific BAM files...")
        report_file.flush()
        report_message, summary_row = create_haplotype_specific_bams(
            filtered_bam, genome, gene_outdir, tempdir, name, gene, region,
            phasing_threshold, threads, remove_monoexonic)
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
    if not os.path.exists(args.vcf):
        raise FileNotFoundError(f"VCF file {args.vcf} not found.")
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
                    phase_reads, args.bam, gene, f"{chrom}:{start}-{end}", args.vcf,
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
    else:
        row = phase_reads(args.bam, args.name, args.region, args.vcf,
                    args.min_distance_from_read_end, args.terminal_variant_proportion,
                    ignore_variants_bed, snvs_only, args.phasing_threshold, args.name,
                    args.genome, args.outdir, args.threads)
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