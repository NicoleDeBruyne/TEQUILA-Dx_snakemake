#!/usr/bin/env python3

# Author: adapted from Nicole DeBruyne (Lin Lab)
# Optimized: 2025

# Extracts all unique splice junctions and their read counts from an entire BAM file.

import os, argparse, warnings, pysam
import pandas as pd
from pandas.errors import PerformanceWarning

warnings.filterwarnings('ignore', category=PerformanceWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# Pre-compile CIGAR pattern once at module load — shared across all calls
import re
_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Extracts all unique splice junctions (format: chrom_intronStart_intronEnd) '
                    'and their read counts from an entire BAM file.'
    )
    parser.add_argument('--bam', required=True, help='Path to the input BAM file')
    parser.add_argument('--outfile', required=True, help='Path to the output TSV file')
    return parser.parse_args()


def get_splice_junctions_from_CIGAR(cigar, chrom, start):
    """Extract splice junctions from a single alignment.

    Returns a list of junction strings: chrom_intronStart_intronEnd
    """
    pos = start
    sj_list = []

    for length_str, operation in _CIGAR_RE.findall(cigar):
        length = int(length_str)

        if operation == 'P':
            warnings.warn(f"Padding operation 'P' encountered in CIGAR string: {cigar}")
        elif operation in 'SHI':
            pass
        elif operation in 'MXD=':
            pos += length
        elif operation == 'N':
            sj_list.append(f"{chrom}_{pos + 1}_{pos + length}")
            pos += length

    return sj_list


def extract_all_splice_junctions(bamfile):
    """Iterate over every alignment in the BAM file and collect all splice junctions.

    Returns:
        sj_dict             : dict mapping junction string -> read count
        total_alignments    : total number of non-unmapped alignments processed
        junction_alignments : number of alignments that contained at least one junction
    """
    sj_dict = {}
    total_alignments = 0
    junction_alignments = 0

    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        for alignment in bam.fetch(until_eof=True):
            cigar = alignment.cigarstring
            if alignment.is_unmapped or cigar is None:
                continue

            total_alignments += 1
            # Use cigartuples (pre-parsed by pysam) instead of re-parsing the string
            has_junction = False
            pos = alignment.reference_start
            chrom = alignment.reference_name

            for op, length in alignment.cigartuples:
                # pysam cigar op codes: 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X
                if op in (0, 2, 7, 8):   # M, D, =, X — consumes reference
                    pos += length
                elif op == 3:             # N — intron (splice junction)
                    sj = f"{chrom}_{pos + 1}_{pos + length}"
                    sj_dict[sj] = sj_dict.get(sj, 0) + 1
                    pos += length
                    has_junction = True
                # op 1 (I), 4 (S), 5 (H) — do not consume reference; skip
                # op 6 (P) — padding; skip

            if has_junction:
                junction_alignments += 1

    return sj_dict, total_alignments, junction_alignments


def main():
    print("\n\n\n******************************************************************************************")
    print("Extracting all splice junctions from BAM file...")
    print("******************************************************************************************\n")

    args = parse_args()

    if not os.path.isfile(args.bam):
        raise FileNotFoundError(f"BAM file not found: {args.bam}")

    outdir = os.path.dirname(os.path.abspath(args.outfile))
    os.makedirs(outdir, exist_ok=True)

    print(f"Input BAM : {args.bam}")
    print(f"Output    : {args.outfile}\n")

    print("Scanning alignments...")
    sj_dict, total_alignments, junction_alignments = extract_all_splice_junctions(args.bam)

    print(f"  Processed alignments       : {total_alignments:,}")
    print(f"  Alignments with junctions  : {junction_alignments:,}")
    print(f"  Unique junctions found     : {len(sj_dict):,}")

    bam_name = os.path.basename(args.bam)
    df = pd.DataFrame(list(sj_dict.items()), columns=['junction', bam_name])
    df.sort_values('junction', inplace=True)

    df.to_csv(args.outfile, sep='\t', index=False)
    print(f"\nResults saved to {args.outfile}\n")


if __name__ == '__main__':
    main()
