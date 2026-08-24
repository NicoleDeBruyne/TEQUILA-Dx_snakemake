#!/usr/bin/env python

# Split out of compile_variants.py so that changing a final filter
# threshold (gnomAD AF, CLNSIG, CADD, SpliceAI, DP, AF) doesn't require
# re-running annotation (ANNOVAR, gnomAD, ClinVar, CADD, SpliceAI) --
# this script just reads compile_variants.py's
# {outprefix}_compiled_variants.tsv output and re-filters it.

import argparse
import pandas as pd


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Apply final filtering thresholds to a variants tsv produced by "
                    "compile_variants.py ({outprefix}_compiled_variants.tsv), and write "
                    "the filtered result to a new tsv. Runs entirely off that already-"
                    "annotated file -- no VCFs, BED, ANNOVAR, gnomAD, ClinVar, CADD, or "
                    "SpliceAI resources needed, so re-filtering with new thresholds is fast.")
    parser.add_argument("--compiled-tsv", required=True,
        help="Path to the {outprefix}_compiled_variants.tsv file written by compile_variants.py.")
    parser.add_argument("--outfile", required=True,
        help="Path to write the filtered tsv to.")
    parser.add_argument("--final-gnomadAF-threshold", type=float,
        help="Only keep variants with gnomAD allele frequency <= threshold.")
    parser.add_argument("--final-CLNSIG-filter", nargs="+",
        help="Remove variants with specified CLNSIG values (ex: Benign, Likely_benign, Benign/Likely_benign).")
    parser.add_argument("--final-CADD-phred-threshold", type=int,
        help="Only keep variants with CADD phred score >= threshold (or SpliceAI score >= threshold, "
             "if --final-SpliceAI-threshold is set, or CADD_PHRED missing/unscored).")
    parser.add_argument("--final-SpliceAI-threshold", type=float,
        help="Only keep variants with SpliceAI score >= threshold (or CADD phred score >= threshold, "
             "if --final-CADD-phred-threshold is set, or CADD_PHRED missing/unscored). A missing/"
             "unscored SpliceAI value does NOT count as a pass on its own -- it only rescues a "
             "variant when CADD is also missing, or when CADD is present but low and SpliceAI is "
             "high.")
    parser.add_argument("--final-DP-threshold", type=int,
        help="Only keep variants with depth (DP) >= threshold.")
    parser.add_argument("--final-AF-threshold", type=float,
        help="Only keep variants with allele frequency (AF/VAF) >= threshold.")
    parser.add_argument('--keep-CLNSIG', nargs="+",
        default=['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'],
        help='Keep variants with these CLNSIG values in the final output regardless of all other '
             'filters. Default keeps all pathogenic variants: Pathogenic, Likely_pathogenic, '
             'Pathogenic/Likely_pathogenic')
    return parser.parse_args()


def filter_variants(input_df, args):
    """Filter variants based on user-defined criteria."""
    filtered_df = input_df.copy()

    if args.final_gnomadAF_threshold:
        if 'gnomAD_AF' not in filtered_df.columns:
            print(f"\nWARNING: gnomAD_AF column not found in DataFrame. Skipping filter.")
        else:
            gnomad_num = pd.to_numeric(filtered_df['gnomAD_AF'], errors='coerce')
            filtered_df = filtered_df[
                (gnomad_num <= args.final_gnomadAF_threshold) | filtered_df['gnomAD_AF'].isin(['.', 'fail'])
            ]

    if args.final_CLNSIG_filter:
        if 'CLNSIG' not in filtered_df.columns:
            print(f"\nWARNING: CLNSIG column not found in DataFrame. Skipping filter.")
        else:
            filtered_df = filtered_df[~filtered_df['CLNSIG'].isin(args.final_CLNSIG_filter)]

    if args.final_CADD_phred_threshold or args.final_SpliceAI_threshold:
        # A variant is kept if: its CADD score is missing/unscored (we never
        # want to drop something CADD couldn't evaluate), OR its CADD score
        # is high enough, OR its SpliceAI score is high enough. A missing
        # SpliceAI score is NOT treated as a pass on its own -- it only
        # matters as a tiebreaker when CADD is present but low. So a variant
        # only gets dropped when CADD was actually computed and is low, and
        # SpliceAI is either also low or missing.
        cadd_missing = pd.Series(False, index=filtered_df.index)
        cadd_high = pd.Series(False, index=filtered_df.index)
        spliceai_high = pd.Series(False, index=filtered_df.index)

        if args.final_CADD_phred_threshold:
            if 'CADD_PHRED' not in filtered_df.columns:
                print(f"\nWARNING: CADD_PHRED column not found in DataFrame. Skipping CADD phred filter.")
            else:
                cadd_num = pd.to_numeric(filtered_df['CADD_PHRED'], errors='coerce')
                cadd_missing = filtered_df['CADD_PHRED'].isin(['.', 'fail']) | cadd_num.isna()
                cadd_high = cadd_num >= args.final_CADD_phred_threshold

        if args.final_SpliceAI_threshold:
            if 'SpliceAI' not in filtered_df.columns:
                print(f"\nWARNING: SpliceAI column not found in DataFrame. Skipping SpliceAI filter.")
            else:
                thr = args.final_SpliceAI_threshold
                def _spliceai_score(x):
                    """Max of the four DS_* SpliceAI delta scores, or NaN if
                    missing/unscored -- NaN compares False against >= thr,
                    so a missing SpliceAI score never counts as 'high'."""
                    if x in ('.', 'fail', '', None):
                        return float('nan')
                    vals = [pd.to_numeric(s, errors='coerce') for s in str(x).split('|')[2:6]]
                    vals = [v for v in vals if pd.notna(v)]
                    return max(vals) if vals else float('nan')
                spliceai_num = filtered_df['SpliceAI'].apply(_spliceai_score)
                spliceai_high = spliceai_num >= thr

        if args.final_CADD_phred_threshold and args.final_SpliceAI_threshold:
            keep_mask = cadd_missing | cadd_high | spliceai_high
        elif args.final_CADD_phred_threshold:
            keep_mask = cadd_missing | cadd_high
        else:
            keep_mask = spliceai_high

        filtered_df = filtered_df[keep_mask]

    if args.final_DP_threshold:
        if 'format' not in filtered_df.columns or 'value' not in filtered_df.columns:
            print(f"\nWARNING: format or value column not found in DataFrame. Skipping filter.")
        else:
            # Vectorised DP extraction
            def _dp_ok(row):
                fmt = row['format'].split(':')
                if 'DP' not in fmt:
                    return True
                val = row['value'].split(':')[fmt.index('DP')]
                dp = pd.to_numeric(val, errors='coerce')
                return pd.isna(dp) or dp >= args.final_DP_threshold
            filtered_df = filtered_df[filtered_df.apply(_dp_ok, axis=1)]

    if args.final_AF_threshold:
        if 'format' not in filtered_df.columns or 'value' not in filtered_df.columns:
            print(f"\nWARNING: format or value column not found in DataFrame. Skipping filter.")
        else:
            def _af_ok(row):
                fmt = row['format'].split(':')
                vals = row['value'].split(':')
                for tag in ('AF', 'VAF'):
                    if tag in fmt:
                        v = pd.to_numeric(vals[fmt.index(tag)], errors='coerce')
                        if not (pd.isna(v) or v >= args.final_AF_threshold):
                            return False
                return True
            filtered_df = filtered_df[filtered_df.apply(_af_ok, axis=1)]

    # Ensure variants with keep-CLNSIG values are preserved regardless of other filters
    def extract_CLNSIG_from_CLNSIGCONF(CLNSIG):
        if pd.isna(CLNSIG):
            return []
        if CLNSIG.startswith("Conflicting_classifications_of_pathogenicity:"):
            return [t.strip() for t in [x.split("(")[0] for x in CLNSIG.split(":", 1)[1].split("|")] if t.strip()]
        return [CLNSIG]

    if args.keep_CLNSIG:
        if 'CLNSIG' not in filtered_df.columns:
            print(f"\nWARNING: CLNSIG column not found in DataFrame. Skipping filter.")
        else:
            keep_set = set(args.keep_CLNSIG)
            keep_mask = input_df['CLNSIG'].apply(
                lambda x: any(term in keep_set for term in extract_CLNSIG_from_CLNSIGCONF(x))
            )
            keep_df = input_df[keep_mask].copy()
            filtered_df = pd.concat([filtered_df, keep_df]).drop_duplicates()

    return filtered_df


def main():
    """Main function."""

    print(f"\n\n\n******************************************************************************************")
    print(f"Filtering variants...")
    print(f"******************************************************************************************\n")

    args = parse_args()

    df = pd.read_csv(args.compiled_tsv, sep='\t', dtype=str)

    if (args.final_gnomadAF_threshold or args.final_CLNSIG_filter or
            args.final_CADD_phred_threshold or args.final_SpliceAI_threshold or
            args.final_DP_threshold or args.final_AF_threshold):
        filtered_df = filter_variants(df, args)
    else:
        # No thresholds configured at all -- nothing to filter, pass
        # everything through untouched (keep-CLNSIG's default value alone
        # shouldn't be treated as "no filters configured", but with no
        # other threshold active there's nothing for it to rescue from,
        # so this is equivalent).
        print(f"\nNo final filter thresholds configured. Writing input through unfiltered.")
        filtered_df = df

    filtered_df.to_csv(args.outfile, sep='\t', index=False)
    print(f"\nFiltered results written to: {args.outfile}")


if __name__ == '__main__':
    main()
