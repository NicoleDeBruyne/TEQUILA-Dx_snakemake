#! /usr/bin/env python3

# scripts/merge_hits.py
# Builds one sample's ranked candidate-hits table from that sample's
# variant/ASE/junction/cohort-junction rows. Exposes load_*()/build_hit_table()
# as a library (imported by scripts/merge_group_hits.py, which calls
# build_hit_table() once per sample in a group and concatenates the results)
# as well as a standalone single-sample CLI (below) for manual use/debugging.

import argparse
import os
import re
import pandas as pd
import numpy as np
import warnings

# Ignore FutureWarnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)

# Maps each raw delta source column to its output-column suffix. The
# GTEx-comparison junction_df only ever has 'delta_PSI' (one metric); the
# cohort-comparison cohort_junction_df has all six (see
# scripts/merge_group_hits.py's normalize_cohort_junction_df).
_DELTA_COL_TO_SUFFIX = {
    'delta_PSI':        'deltaPSI',
    'delta_PSI_approx': 'deltaPSIapprox',
    'delta_5ss_IR':     'delta5ssIR',
    'delta_3ss_IR':     'delta3ssIR',
    'delta_full_IR':    'deltaFullIR',
    'delta_IPA':        'deltaIPA',
}

# Which delta-metric suffixes exist per prefix -- '' (GTEx-comparison) only
# ever has deltaPSI; 'cohort_' has all six. max_deltas() (inside
# build_hit_table) takes the max magnitude across whichever of these are
# present for a given prefix: no per-metric weighting, just the single
# largest |delta| across all populated metrics.
_DELTA_SUFFIXES_BY_PREFIX = {
    '':        ('deltaPSI',),
    'cohort_': ('deltaPSI', 'deltaPSIapprox', 'delta5ssIR', 'delta3ssIR', 'deltaFullIR', 'deltaIPA'),
}


def parse_args():
    """ Parse command line arguments """

    parser = argparse.ArgumentParser(description='Merge hits from variant, ASE, and junction analyses for a single sample.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file')
    parser.add_argument('--sample-name', required=True, type=str, help='Name of the sample')
    parser.add_argument('--variant-hits', required=True, type=str, help='Path to variant hit results. At a minimum, file should contain columns: \
                                                                            sample, chrom, pos, ref, alt, GT, gnomAD_AF, CLNSIG, gene, CADD_PHRED, SpliceAI, num_callers, sample_count, \
                                                                            ANNOVAR_AAChange.refGene, ANNOVAR_GeneDetail.refGene')
    parser.add_argument('--ase-hits', required=True, type=str, help='Path to ASE hit results. At a minimum, file should contain columns: sample, gene, ratio')
    parser.add_argument('--junction-hits', required=True, type=str, help='Path to junction hit results. At a minimum, file should contain columns: \
                                                                            sample, gene, phasing, junction, delta_PSI, sample_count, annotation, event')
    parser.add_argument('--cohort-junction-hits', required=False, default=None, type=str, help='Path to cohort-comparison junction hit results (rules/7_cohort_junction_analysis.smk), \
                                                                            same schema as --junction-hits but the sample compared against the rest of its cohort rather than GTEx. \
                                                                            If omitted, the cohort_* junction columns are filled with "." rather than annotated.')
    parser.add_argument('--omim', required=False, default=None, type=str, help='Path to OMIM data. At a minimum, file should contain columns: approved_gene_symbol, phenotypes, inheritance_patterns, haploinsufficient (TRUE/FALSE/empty; empty treated as FALSE). \
                                                                                  If omitted, phenotypes/inheritance_patterns are filled with "." and haploinsufficient is treated as False.')
    
    return parser.parse_args()


# Caller priority for phase resolution, highest-trust first. Only these
# three callers ever phase (see rules/1_call_variants.smk: nanoTS has a
# model_phased option, clair3_rna runs --enable_phasing_model via
# whatshap, longcallR does native joint calling+phasing) -- deepvariant
# never phases and is intentionally absent here, so a deepvariant-only
# call for a variant simply can't contribute phase evidence.
_CALLER_PHASE_PRIORITY = ('nanoTS', 'clair3_rna', 'longcallR')


def _caller_from_name(name):
    """compile_variants.py's --vcf-file-names are '{sample}_{caller}'
    (e.g. 'IDT-160_PID843_nanoTS') -- match against the known caller
    tokens rather than splitting on '_', since sample names themselves
    may contain underscores."""
    for caller in _CALLER_PHASE_PRIORITY + ('deepvariant',):
        if name.endswith('_' + caller):
            return caller
    return None


def resolve_phase(variant_a, variant_b):
    """Determine the trans/cis/unclear relationship between two variants
    in the same gene, using each caller's raw (unmodified) GT + PS.

    variant_a, variant_b: each a dict/Series-like with '_caller', '_raw_GT',
    'PS' for one (caller, variant) row. Callers are checked in priority
    order (nanoTS > clair3_rna > longcallR, see _CALLER_PHASE_PRIORITY) --
    the first caller that phased BOTH variants against the SAME PS wins;
    no fallback to a lower-priority caller once a higher one has already
    given an answer, even if that answer conflicts with what a
    lower-priority caller would have said.

    A caller's phasing for a variant is only trusted if GT is phased
    ('|' present) AND PS is present and not '.' -- PS is meaningless on
    an unphased GT, and comparing '.' PS values would be a false match.
    PS is only comparable within the same caller (see compile_variants.py's
    process_vcf); this function only ever compares same-caller PS pairs
    by construction, since it iterates one caller at a time.

    Returns: 'trans', 'cis', or 'unclear'.
    """
    def _phased_calls(variant, caller):
        """All (raw_GT, PS) pairs for `variant` reported by `caller`,
        restricted to genuinely phased, non-'.' PS entries. A variant can
        have >1 row from the same caller only in unusual multi-record
        VCF situations; in the normal case this is 0 or 1 entries."""
        rows = variant if isinstance(variant, list) else [variant]
        out = []
        for row in rows:
            if row.get('_caller') != caller:
                continue
            gt = row.get('_raw_GT', '.')
            ps = row.get('PS', '.')
            if '|' in str(gt) and ps not in (None, '.', ''):
                out.append((gt, ps))
        return out

    for caller in _CALLER_PHASE_PRIORITY:
        calls_a = _phased_calls(variant_a, caller)
        calls_b = _phased_calls(variant_b, caller)
        for gt_a, ps_a in calls_a:
            for gt_b, ps_b in calls_b:
                if ps_a != ps_b:
                    continue
                # Same phase set from the same caller -- allele order is
                # directly comparable. hap-1-allele is the first GT digit.
                hap1_a = gt_a.split('|')[0]
                hap1_b = gt_b.split('|')[0]
                return 'cis' if hap1_a == hap1_b else 'trans'
        # This caller phased neither variant (or phased them into
        # different/unrelated phase sets) -- fall through to the next
        # caller in priority order rather than deciding 'unclear' yet.

    return 'unclear'


def load_variant_df(path):
    """Read + preprocess a variant-hits tsv (single-sample or group-level --
    same schema either way, just more rows for the latter)."""
    variant_df = pd.read_csv(path, sep='\t', usecols=[
        'sample', 'chrom', 'pos', 'ref', 'alt', 'name', 'GT', 'PS', 'gnomAD_AF', 'CLNSIG', 'gene', 'CADD_PHRED', 'SpliceAI',
        'num_callers', 'sample_count', 'ANNOVAR_AAChange.refGene', 'ANNOVAR_GeneDetail.refGene',
    ]).drop_duplicates()
    variant_df = variant_df.rename(columns={'sample_count': 'variant_nsamples'})
    variant_df = variant_df[variant_df['gene'] != '.']
    # 'gene' lists every BED-panel gene a variant overlaps (comma-separated
    # if more than one) -- explode so a variant overlapping multiple genes
    # is grouped into each of those genes' hit rows individually.
    variant_df['gene'] = variant_df['gene'].str.split(',')
    variant_df = variant_df.explode('gene')
    # Raw, phase-preserving copy -- 'name' is "{sample}_{caller}" (see
    # compile_variants.py's --vcf-file-names), so this identifies which
    # caller reported this specific GT/PS pair. resolve_phase() below
    # needs the untouched '|' separator and caller identity to determine
    # trans/cis; keep this BEFORE the cosmetic canonicalization that
    # follows, which is for the human-readable variant_GT display column
    # only and would otherwise destroy phase information.
    variant_df['_caller'] = variant_df['name'].apply(_caller_from_name)
    variant_df['_raw_GT'] = variant_df['GT']
    variant_df['GT'] = variant_df['GT'].str.replace('|', '/', regex=False).str.replace('1/0', '0/1', regex=False)
    variant_df['variant_ID'] = variant_df.apply(lambda x: f"{x.chrom}-{x.pos}-{x.ref}-{x.alt}", axis=1)
    # ANNOVAR only ever populates one of these two per variant (AAChange.refGene
    # for exonic variants with a codon change to report, GeneDetail.refGene for
    # everything else it has position detail for -- splicing/UTR/intronic/etc.)
    # -- coalesce into a single human-readable consequence column rather than
    # carrying two mostly-empty columns through to the final hits table.
    aachange = variant_df['ANNOVAR_AAChange.refGene']
    genedetail = variant_df['ANNOVAR_GeneDetail.refGene']
    variant_df['variant_consequence'] = aachange.where(aachange.notna() & (aachange != '.'), genedetail)
    variant_df['variant_consequence'] = variant_df['variant_consequence'].fillna('.')
    return variant_df


def load_ase_df(path):
    """Read + preprocess an ASE-hits tsv (single-sample or group-level)."""
    ase_df = pd.read_csv(path, sep='\t', usecols=['sample', 'gene', 'ratio', 'sample_count']).drop_duplicates()
    return ase_df.rename(columns={'ratio': 'ASE_ratio', 'sample_count': 'ASE_nsamples'})


def load_junction_df(path):
    """Read a GTEx-comparison junction-hits tsv (single-sample or group-level) -- one metric (delta_PSI)."""
    return pd.read_csv(path, sep='\t', usecols=[
        'sample', 'gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event',
    ]).drop_duplicates()


def load_cohort_junction_df(path):
    """Read a cohort-comparison junction-hits tsv (already normalized to one
    delta column per metric by scripts/merge_group_hits.py's
    normalize_cohort_junction_df -- six metrics, unlike load_junction_df's one)."""
    return pd.read_csv(path, sep='\t', usecols=[
        'sample', 'gene', 'phasing', 'junction', 'jxn_coverage',
        'delta_PSI', 'delta_PSI_approx', 'delta_5ss_IR', 'delta_3ss_IR', 'delta_full_IR', 'delta_IPA',
        'sample_count', 'annotation', 'event',
    ]).drop_duplicates()


def load_omim_df(path):
    omim_df = pd.read_csv(path, sep='\t', usecols=[
        'approved_gene_symbol', 'phenotypes', 'inheritance_patterns', 'haploinsufficient',
    ])
    omim_df = omim_df.rename(columns={'approved_gene_symbol': 'gene'})
    # TRUE/FALSE (any case) or empty/missing -- empty is treated the same
    # as FALSE (see setup.sh's OMIM section for the full column spec).
    omim_df['haploinsufficient'] = (
        omim_df['haploinsufficient'].astype(str).str.strip().str.upper() == 'TRUE'
    )
    return omim_df


def build_phased_junction_df(df, prefix, delta_cols=('delta_PSI',)):
    """Aggregate a junction_df-shaped table (gene, phasing, junction,
    jxn_coverage, one-or-more delta columns, sample_count, annotation,
    event) into one gene-level row per phasing tier, with columns named
    '{prefix}bulk_jxns', '{prefix}bulk_jxn_coverage', one
    '{prefix}bulk_{suffix}' per entry in delta_cols, etc. Used for both
    the GTEx-comparison junction_df (prefix='', delta_cols=('delta_PSI',))
    and the cohort-comparison cohort_junction_df (prefix='cohort_',
    delta_cols=all six metrics)."""
    tiers = {}
    for phasing, sep in (('bulk', ';'), ('hap1', ','), ('hap2', ',')):
        agg_kwargs = {
            prefix + phasing + '_jxns':         ('junction', lambda x, sep=sep: sep.join(map(str, x))),
            prefix + phasing + '_jxn_coverage': ('jxn_coverage', lambda x, sep=sep: sep.join(map(str, x))),
        }
        for col in delta_cols:
            suffix = _DELTA_COL_TO_SUFFIX[col]
            agg_kwargs[prefix + phasing + '_' + suffix] = (col, lambda x, sep=sep: sep.join(map(str, x)))
        agg_kwargs[prefix + phasing + '_jxn_annotation'] = ('annotation', lambda x, sep=sep: sep.join(map(str, x)))
        agg_kwargs[prefix + phasing + '_jxn_event']      = ('event', lambda x, sep=sep: sep.join(map(str, x)))
        agg_kwargs[prefix + phasing + '_jxn_nsamples']   = ('sample_count', lambda x, sep=sep: sep.join(map(str, x)))
        tiers[phasing] = (
            df[df['phasing'] == phasing]
                .sort_values('junction')
                .groupby('gene')
                .agg(**agg_kwargs)
                .reset_index()
        )
    merged = pd.merge(tiers['bulk'], tiers['hap1'], on='gene', how='outer')
    merged = pd.merge(merged, tiers['hap2'], on='gene', how='outer')
    return merged.drop_duplicates()


def load_gene_expression_df(path):
    """Read the CPTM matrix from rule _10D (quantify_gene_by_assignment.py's
    <outprefix>_matrix.tsv) -- one row per targeted-panel gene, one column
    per sample in that (bed_id, sample_type) cohort. Since this matrix is
    already scoped to exactly one (panel, sample_type) group, every
    sample column in it IS the cohort for percentile/n_cohort purposes --
    no separate cohort-membership lookup needed."""
    return pd.read_csv(path, sep='\t', index_col=0)


def build_hit_table(variant_df, ase_df, junction_df, cohort_junction_df, sample_name, omim_df=None, gene_expression_df=None):
    """Build one sample's ranked candidate-hits table. All four input
    DataFrames are assumed already loaded (via load_*() above) and already
    filtered down to this sample only -- this function itself is agnostic
    to whether they came from single-sample files or were filtered out of
    group-level tables. cohort_junction_df/omim_df may be None (same
    fallback behavior as omitting --cohort-junction-hits/--omim on the CLI:
    cohort_* columns filled with '.', phenotypes/inheritance_patterns
    filled with '.')."""

    # Create modified DataFrames
    # mod_variant_df should contain genes, variant_pos, variant_GT, variant_CLNSIG, variant_nsamples
    mod_variant_df = (
        variant_df
            .sort_values(['gene', 'variant_ID'])
            .groupby(['gene', 'variant_ID'], sort=False)
            .agg(
                GT=('GT', lambda x: ','.join(x)),
                gnomAD_AF=('gnomAD_AF', lambda x: ','.join(x.dropna().astype(str))),
                CLNSIG=('CLNSIG', lambda x: ','.join(x.dropna().astype(str))),
                CADD_PHRED=('CADD_PHRED', lambda x: ','.join(x.dropna().astype(str))),
                SpliceAI=('SpliceAI', lambda x: ','.join(x.dropna().astype(str))),
                variant_consequence=('variant_consequence', lambda x: ','.join(x.dropna().astype(str))),
                num_callers=('num_callers', lambda x: ','.join(x.astype(str))),
                variant_nsamples=('variant_nsamples', lambda x: ','.join(x.astype(str))),
            )
            .reset_index()
            .groupby('gene', sort=False)
            .agg(
                variant_ID=('variant_ID', ';'.join),
                variant_GT=('GT', ';'.join),
                variant_gnomAD_AF=('gnomAD_AF', ';'.join),
                variant_CLNSIG=('CLNSIG', ';'.join),
                variant_CADD_PHRED=('CADD_PHRED', ';'.join),
                variant_SpliceAI=('SpliceAI', ';'.join),
                variant_consequence=('variant_consequence', ';'.join),
                variant_num_callers=('num_callers', ';'.join),
                variant_nsamples=('variant_nsamples', ';'.join),
            )
            .reset_index()
    )
    # mod_junction_df should contain genes, bulk_jxns, bulk_jxn_coverage, bulk_deltaPSI, bulk_jxn_nsamples, 
    # hap1_jxns, hap1_jxn_coverage, hap1_deltaPSI, hap1_jxn_nsamples, hap2_jxns, hap2_jxn_coverage, hap2_delta_PSI, hap2_jxn_nsamples
    mod_junction_df = build_phased_junction_df(junction_df, '')
    if cohort_junction_df is not None:
        mod_cohort_junction_df = build_phased_junction_df(
            cohort_junction_df, 'cohort_',
            delta_cols=('delta_PSI', 'delta_PSI_approx', 'delta_5ss_IR', 'delta_3ss_IR', 'delta_full_IR', 'delta_IPA'),
        )
    else:
        mod_cohort_junction_df = None

    # Merge hits
    hit_df = pd.merge(mod_variant_df, ase_df, on='gene', how='outer')
    hit_df = pd.merge(hit_df, mod_junction_df, on='gene', how='outer')
    if mod_cohort_junction_df is not None:
        hit_df = pd.merge(hit_df, mod_cohort_junction_df, on='gene', how='outer')
    else:
        # No cohort-comparison junction data provided -- keep the same
        # output schema, just unannotated, rather than dropping these
        # columns entirely (mirrors the phenotypes/inheritance_patterns
        # fallback below when omim_df is omitted).
        for phasing in ('bulk', 'hap1', 'hap2'):
            for suffix in ('jxns', 'jxn_coverage', 'deltaPSI', 'deltaPSIapprox', 'delta5ssIR',
                           'delta3ssIR', 'deltaFullIR', 'deltaIPA', 'jxn_annotation', 'jxn_event', 'jxn_nsamples'):
                hit_df['cohort_' + phasing + '_' + suffix] = '.'
    if omim_df is not None:
        hit_df = pd.merge(hit_df, omim_df, on='gene', how='left')
        # Genes not present in the OMIM table get NaN from the left merge --
        # treat missing haploinsufficiency data the same as FALSE, same
        # convention as an explicitly empty cell (see setup.sh's OMIM
        # column spec).
        hit_df['haploinsufficient'] = hit_df['haploinsufficient'].fillna(False)
    else:
        # No OMIM data provided -- keep the same output schema, just
        # unannotated, rather than dropping these columns entirely.
        hit_df['phenotypes'] = '.'
        hit_df['inheritance_patterns'] = '.'
        hit_df['haploinsufficient'] = False
    hit_df = hit_df.drop_duplicates()

    # Create boolean columns to indicate whether there is a candidate variant or allele-specific expression
    hit_df["variant"] = hit_df["variant_ID"].notna()
    hit_df["pathogenic_variant"] = (
        hit_df["variant_CLNSIG"]
            .astype(str)
            .str.contains(r"Pathogenic|Likely_pathogenic", regex=True)
    )
    hit_df["ASE"] = pd.to_numeric(hit_df["ASE_ratio"], errors="coerce").notna()

    # Add a column to indicate whether a gene has strong, moderate, weak, splicing dysregulation

    def max_deltas(row, prefix=''):
        """ Extract the max |delta| across all populated delta metrics for a gene, per phasing tier """
        def parse_vals(s):
            if not pd.notna(s):
                return []
            vals = []
            for v in re.split('[,;]', s):
                v = v.strip()
                if not v:
                    continue
                try:
                    vals.append(float(v))
                except ValueError:
                    # Non-numeric sentinel from rules/7_cohort_junction_analysis.smk's
                    # identify_cohort_junction_outliers.py -- "low_n" (too few
                    # cohort samples with good coverage), "error" (beta_binomial
                    # fit failed), or "no_variance" (modified_zscore: zero
                    # variance in the reference distribution). These mean "not
                    # statistically testable", not a numeric value -- skip
                    # rather than crash, same as a blank/missing entry.
                    continue
            return vals

        def max_for_phasing(phasing):
            all_vals = []
            for suffix in _DELTA_SUFFIXES_BY_PREFIX[prefix]:
                all_vals.extend(parse_vals(row.get(prefix + phasing + '_' + suffix)))
            return max(abs(v) for v in all_vals) if all_vals else np.nan

        max_bulk = max_for_phasing('bulk')
        max_hap1 = max_for_phasing('hap1')
        max_hap2 = max_for_phasing('hap2')

        return max_bulk, max_hap1, max_hap2

    def inspect_row(row, prefix=''):
        """ Determine if a gene has splicing dysregulation """
        max_bulk, max_hap1, max_hap2 = max_deltas(row, prefix)
        dominant = any(x in row['inheritance_patterns'] for x in ['AD', 'XLD']) if pd.notna(row.get('inheritance_patterns')) else False
        if (
            (max_bulk >= 0.5)
            or (max_hap1 >= 0.5 and max_hap2 >= 0.5)
            or (dominant and (max_bulk >= 0.2 or max_hap1 >= 0.5 or max_hap2 >= 0.5))
        ):
            return "Strong"
        elif (max_bulk >= 0.2 or max_hap1 >= 0.5 or max_hap2 >= 0.5):
            return "Moderate"
        elif (max_bulk > 0 or max_hap1 > 0 or max_hap2 > 0):
            return "Weak"
        else:
            return "None"
    hit_df['outlier_junction'] = hit_df.apply(inspect_row, axis=1)
    if mod_cohort_junction_df is not None:
        # NOTE: cohort_bulk_deltaPSI (etc.) is delta_junction_PSI for groups
        # analyzed with the beta_binomial method (a PSI-scale value in
        # [-1, 1], same units the 0.2/0.5 thresholds below assume) but
        # modz_junction_PSI for groups analyzed with modified_zscore (an
        # unbounded z-score) -- rules/7_cohort_junction_analysis.smk picks
        # the method per group based on its sample count. Strong/Moderate/
        # Weak is therefore not on a consistent scale across groups that
        # used different methods; treat cross-group comparisons of this
        # column with that in mind.
        hit_df['cohort_outlier_junction'] = hit_df.apply(lambda row: inspect_row(row, prefix='cohort_'), axis=1)
    else:
        hit_df['cohort_outlier_junction'] = '.'

    # Fill missing values
    # Cast to object before fillna: on pandas >=3.0, filling a still-numeric
    # column (e.g. ASE_ratio, which is never string-joined the way the
    # variant/junction columns are) with the string "." raises
    # LossySetitemError instead of silently upcasting like older pandas did.
    hit_df = hit_df.astype(object)
    hit_df.fillna(".", inplace=True)

    # ------------------------------------------------------------------
    # Tier assignment (replaces the old points-based 'score').
    #
    # Tiers are assigned by the decision tree below, branching first on
    # the gene's inheritance pattern (bucket 1/2/3), then walking down a
    # fixed, ordered list of conditions per bucket -- the first condition
    # that matches wins; lower-priority conditions are never reached once
    # a higher one has already matched. Lower tier number = higher
    # priority / more likely to be diagnostic.
    #
    # Bucket 1 -- AD or XLD only (inheritance_patterns has AD/XLD, no AR/XLR):
    #   1: any pathogenic variant present
    #   2: ASE, or a Strong junction outlier (GTEx- or cohort-comparison)
    #   3: a Moderate junction outlier (not currently reachable -- see
    #      inspect_row()'s dominant-gene Strong condition, which subsumes
    #      the Moderate condition entirely for dominant genes; kept as an
    #      explicit branch rather than silently dropped, in case that
    #      ever changes)
    #   4: a Weak junction outlier
    #   5: at least one VUS present (nothing else above matched)
    #
    # Bucket 2 -- AD/AR or XLD/XLR combined (has both a dominant AND a
    # recessive marker in inheritance_patterns):
    #   1: 2 pathogenic variants in trans, or 1 homozygous pathogenic variant
    #   2: any other pathogenic variant present
    #   3: ASE or Strong junction
    #   4: Moderate junction (same reachability caveat as bucket 1)
    #   5: Weak junction
    #   6: at least one VUS present
    #
    # Bucket 3 -- anything else (no dominant marker at all: pure AR/XLR,
    # unknown, or no inheritance_patterns data):
    #   1: 2 pathogenic variants in trans, or 1 homozygous pathogenic variant
    #   2: (1 pathogenic + 1 VUS in trans) or (2 pathogenic, unclear phasing)
    #   3-5: (1 pathogenic alone) or (2 VUS in trans) or (1 homozygous VUS)
    #        -- 3 if ASE/Strong, 4 if Moderate, 5 if Weak, else base tier 5... 
    #        wait: base case here is tier 2 if none of ASE/moderate/weak,
    #        see _assign_tier's bucket-3 branch below for the exact tiers
    #   4-6: 2 VUS, unclear phasing -- 3 if ASE/Strong, 4 if Moderate, 5 if
    #        Weak, else tier 6
    #   4-7: anything else with at least one VUS or orthogonal evidence --
    #        4 if ASE/Strong, 5 if Moderate, 6 if Weak, else 7 if >=1 VUS
    #
    # See _assign_tier() for the exact tier numbers per branch -- the
    # prose above is a summary, the code is authoritative.
    def _clean_inheritance(s):
        return str(s).strip().strip('"') if pd.notna(s) else ''

    def _inheritance_bucket(inheritance_patterns):
        s = _clean_inheritance(inheritance_patterns)
        has_dominant = ('AD' in s) or ('XLD' in s)
        has_recessive = ('AR' in s) or ('XLR' in s)
        if has_dominant and has_recessive:
            return 2
        elif has_dominant:
            return 1
        else:
            return 3

    def _classify_gene_variants(gene_variant_rows):
        """gene_variant_rows: raw (pre-aggregation) variant_df rows for
        ONE gene -- one row per (caller, physical variant). Groups those
        rows by variant_ID (chrom-pos-ref-alt) to get one entry per
        distinct physical variant, classifies each as pathogenic/VUS and
        homozygous/het, then checks every pathogenic/VUS pair for
        trans/cis/unclear phasing via resolve_phase() (which itself
        applies the nanoTS > clair3_rna > longcallR caller-priority
        chain). Checking every pair (not just the two highest-scoring
        variants individually) is what lets e.g. a phased pathogenic+VUS
        pair beat two pathogenic variants that turned out to be in cis."""
        empty = dict(
            n_pathogenic_variants=0, has_pathogenic=False, has_vus=False,
            has_homozygous_pathogenic=False, has_homozygous_vus=False,
            has_2_pathogenic_trans=False, has_2_pathogenic_unclear=False,
            has_1_pathogenic_1_vus_trans=False,
            has_2_vus_trans=False, has_2_vus_unclear=False,
        )
        if gene_variant_rows.empty:
            return empty

        variants = {}
        for variant_id, sub in gene_variant_rows.groupby('variant_ID'):
            is_pathogenic = sub['CLNSIG'].astype(str).str.contains(r'Pathogenic|Likely_pathogenic', regex=True).any()
            is_homozygous = (sub['GT'] == '1/1').any()
            variants[variant_id] = dict(
                is_pathogenic=is_pathogenic,
                is_homozygous=is_homozygous,
                rows=sub[['_caller', '_raw_GT', 'PS']].to_dict('records'),
            )

        path_ids = [v for v, d in variants.items() if d['is_pathogenic']]
        vus_ids = [v for v, d in variants.items() if not d['is_pathogenic']]

        has_2_pathogenic_trans = False
        has_2_pathogenic_unclear = False
        for i in range(len(path_ids)):
            for j in range(i + 1, len(path_ids)):
                phase = resolve_phase(variants[path_ids[i]]['rows'], variants[path_ids[j]]['rows'])
                has_2_pathogenic_trans |= (phase == 'trans')
                has_2_pathogenic_unclear |= (phase == 'unclear')

        has_1_pathogenic_1_vus_trans = False
        for p in path_ids:
            for v in vus_ids:
                if resolve_phase(variants[p]['rows'], variants[v]['rows']) == 'trans':
                    has_1_pathogenic_1_vus_trans = True

        has_2_vus_trans = False
        has_2_vus_unclear = False
        for i in range(len(vus_ids)):
            for j in range(i + 1, len(vus_ids)):
                phase = resolve_phase(variants[vus_ids[i]]['rows'], variants[vus_ids[j]]['rows'])
                has_2_vus_trans |= (phase == 'trans')
                has_2_vus_unclear |= (phase == 'unclear')

        return dict(
            n_pathogenic_variants=len(path_ids),
            has_pathogenic=len(path_ids) > 0,
            has_vus=len(vus_ids) > 0,
            has_homozygous_pathogenic=any(variants[v]['is_homozygous'] for v in path_ids),
            has_homozygous_vus=any(variants[v]['is_homozygous'] for v in vus_ids),
            has_2_pathogenic_trans=has_2_pathogenic_trans,
            has_2_pathogenic_unclear=has_2_pathogenic_unclear,
            has_1_pathogenic_1_vus_trans=has_1_pathogenic_1_vus_trans,
            has_2_vus_trans=has_2_vus_trans,
            has_2_vus_unclear=has_2_vus_unclear,
        )

    variant_class_by_gene = {
        gene: _classify_gene_variants(sub)
        for gene, sub in variant_df.groupby('gene')
    }
    _empty_variant_class = dict(
        n_pathogenic_variants=0, has_pathogenic=False, has_vus=False,
        has_homozygous_pathogenic=False, has_homozygous_vus=False,
        has_2_pathogenic_trans=False, has_2_pathogenic_unclear=False,
        has_1_pathogenic_1_vus_trans=False,
        has_2_vus_trans=False, has_2_vus_unclear=False,
    )

    def _assign_tier(row):
        vc = variant_class_by_gene.get(row['gene'], _empty_variant_class)
        ase = row['ASE'] is True
        haploinsufficient = row['haploinsufficient'] is True
        strong = (row['outlier_junction'] == 'Strong') or (row['cohort_outlier_junction'] == 'Strong')
        moderate = (row['outlier_junction'] == 'Moderate') or (row['cohort_outlier_junction'] == 'Moderate')
        weak = (row['outlier_junction'] == 'Weak') or (row['cohort_outlier_junction'] == 'Weak')
        ase_or_strong = ase or strong
        bucket = _inheritance_bucket(row['inheritance_patterns'])

        # For AD/XLD genes specifically, ASE means something different
        # depending on whether the gene is known to be haploinsufficient:
        # for a haploinsufficient gene, one fully-silenced allele (ASE) is
        # itself near-diagnostic, on par with a Strong junction outlier.
        # For a non-haploinsufficient gene, the same ASE call is weaker
        # evidence (the gene tolerates one silenced allele), so it's
        # ranked down with Moderate junction instead. Strong/Moderate
        # junction themselves are unconditional either way -- only ASE's
        # weight depends on haploinsufficient.
        ase_hapi_or_strong = (ase and haploinsufficient) or strong
        ase_not_hapi_or_moderate = (ase and not haploinsufficient) or moderate

        if bucket == 1:  # AD or XLD only
            if vc['has_pathogenic']:
                return 1
            if ase_hapi_or_strong:
                return 2
            if ase_not_hapi_or_moderate:
                return 3
            if weak:
                return 4
            if vc['has_vus']:
                return 5
            return None

        if bucket == 2:  # AD/AR or XLD/XLR combined
            if vc['has_2_pathogenic_trans'] or vc['has_homozygous_pathogenic']:
                return 1
            if vc['has_pathogenic']:
                return 2
            if ase_hapi_or_strong:
                return 3
            if ase_not_hapi_or_moderate:
                return 4
            if weak:
                return 5
            if vc['has_vus']:
                return 6
            return None

        # bucket == 3: anything else (AR/XLR, unknown, or no inheritance data)
        if vc['has_2_pathogenic_trans'] or vc['has_homozygous_pathogenic']:
            return 1
        if vc['has_1_pathogenic_1_vus_trans'] or vc['has_2_pathogenic_unclear']:
            return 2
        if vc['has_pathogenic'] or vc['has_2_vus_trans'] or vc['has_homozygous_vus']:
            if ase_or_strong:
                return 2
            if moderate:
                return 3
            if weak:
                return 4
            return 5
        if vc['has_2_vus_unclear']:
            if ase_or_strong:
                return 3
            if moderate:
                return 4
            if weak:
                return 5
            return 6
        if ase_or_strong:
            return 4
        if moderate:
            return 5
        if weak:
            return 6
        if vc['has_vus']:
            return 7
        return None

    hit_df['tier'] = hit_df.apply(_assign_tier, axis=1)
    hit_df['_tb_n_pathogenic'] = hit_df['gene'].map(
        lambda g: variant_class_by_gene.get(g, _empty_variant_class)['n_pathogenic_variants']
    )

    def _max_bulk_delta(row):
        a = max_deltas(row, '')[0]
        b = max_deltas(row, 'cohort_')[0]
        vals = [v for v in (a, b) if pd.notna(v)]
        return max(vals) if vals else -1

    hit_df['_tb_max_bulk_delta'] = hit_df.apply(_max_bulk_delta, axis=1)

    # Tiebreak chain (applied only when two genes have the same tier),
    # in priority order: # pathogenic variants (descending), presence of
    # ASE (True first), max junction delta magnitude (descending, across
    # both GTEx- and cohort-comparison junctions), then gene name
    # alphabetically as a final stable tiebreak.
    hit_df.sort_values(
        by=['tier', '_tb_n_pathogenic', 'ASE', '_tb_max_bulk_delta', 'gene'],
        ascending=[True, False, False, False, True],
        na_position='last',
        inplace=True
    )
    hit_df.drop(columns=['_tb_n_pathogenic', '_tb_max_bulk_delta'], inplace=True)
    hit_df['ranking'] = np.arange(1, len(hit_df) + 1)

    # Gene expression (rule _10D's targeted-panel CPTM matrix, one row per
    # gene / one column per sample in this exact (bed_id, sample_type)
    # cohort -- see load_gene_expression_df()). relative_gene_expression is
    # this sample's own CPTM value for the gene; cohort_relative_gene_expression
    # is that gene's (1st, 50th, 99th) percentile across every sample in the
    # matrix (i.e. the whole cohort, since the matrix is already scoped to
    # this one group); n_cohort is the total sample count backing those
    # percentiles (same for every gene/row -- the matrix's column count).
    if gene_expression_df is not None:
        n_cohort = gene_expression_df.shape[1]

        def _relative_expression(gene):
            if gene in gene_expression_df.index and sample_name in gene_expression_df.columns:
                return gene_expression_df.loc[gene, sample_name]
            return '.'

        def _cohort_relative_expression(gene):
            if gene not in gene_expression_df.index:
                return '.'
            vals = gene_expression_df.loc[gene].astype(float).values
            p1, p50, p99 = np.percentile(vals, [1, 50, 99])
            return f"({p1:.2f}, {p50:.2f}, {p99:.2f})"

        hit_df['relative_gene_expression'] = hit_df['gene'].apply(_relative_expression)
        hit_df['cohort_relative_gene_expression'] = hit_df['gene'].apply(_cohort_relative_expression)
        hit_df['n_cohort'] = n_cohort
    else:
        # No gene-expression matrix provided -- keep the same output
        # schema, just unannotated (same fallback convention as
        # omim_df/cohort_junction_df being omitted elsewhere in this
        # function).
        hit_df['relative_gene_expression'] = '.'
        hit_df['cohort_relative_gene_expression'] = '.'
        hit_df['n_cohort'] = '.'

    # Reorder
    hit_df['sample'] = sample_name
    hit_df = hit_df[[
        'sample', 'gene', 'phenotypes', 'inheritance_patterns', 'haploinsufficient', 'ranking', 'tier', 'variant', 'pathogenic_variant', 'ASE', 'outlier_junction', 'cohort_outlier_junction',
        'relative_gene_expression', 'cohort_relative_gene_expression', 'n_cohort',
        'variant_ID', 'variant_GT', 'variant_gnomAD_AF',  'variant_CLNSIG', 'variant_CADD_PHRED', 'variant_SpliceAI', 'variant_consequence', 'variant_num_callers', 'variant_nsamples',
        'ASE_ratio', 'ASE_nsamples', 
        'bulk_jxns', 'bulk_jxn_coverage', 'bulk_deltaPSI', 'bulk_jxn_annotation', 'bulk_jxn_event', 'bulk_jxn_nsamples',
        'hap1_jxns', 'hap1_jxn_coverage', 'hap1_deltaPSI', 'hap1_jxn_annotation', 'hap1_jxn_event', 'hap1_jxn_nsamples',
        'hap2_jxns', 'hap2_jxn_coverage', 'hap2_deltaPSI', 'hap2_jxn_annotation', 'hap2_jxn_event', 'hap2_jxn_nsamples',
        'cohort_bulk_jxns', 'cohort_bulk_jxn_coverage',
        'cohort_bulk_deltaPSI', 'cohort_bulk_deltaPSIapprox', 'cohort_bulk_delta5ssIR', 'cohort_bulk_delta3ssIR', 'cohort_bulk_deltaFullIR', 'cohort_bulk_deltaIPA',
        'cohort_bulk_jxn_annotation', 'cohort_bulk_jxn_event', 'cohort_bulk_jxn_nsamples',
        'cohort_hap1_jxns', 'cohort_hap1_jxn_coverage',
        'cohort_hap1_deltaPSI', 'cohort_hap1_deltaPSIapprox', 'cohort_hap1_delta5ssIR', 'cohort_hap1_delta3ssIR', 'cohort_hap1_deltaFullIR', 'cohort_hap1_deltaIPA',
        'cohort_hap1_jxn_annotation', 'cohort_hap1_jxn_event', 'cohort_hap1_jxn_nsamples',
        'cohort_hap2_jxns', 'cohort_hap2_jxn_coverage',
        'cohort_hap2_deltaPSI', 'cohort_hap2_deltaPSIapprox', 'cohort_hap2_delta5ssIR', 'cohort_hap2_delta3ssIR', 'cohort_hap2_deltaFullIR', 'cohort_hap2_deltaIPA',
        'cohort_hap2_jxn_annotation', 'cohort_hap2_jxn_event', 'cohort_hap2_jxn_nsamples'
    ]]
    return hit_df


def main():
    """ Main function -- standalone single-sample CLI. """

    # Parse command line arguments
    args = parse_args()
    for attr, value in vars(args).items():
        if value == "None":
            setattr(args, attr, None)
    print(f"\nMerging variant, ASE, and outlier junction hits for sample: {args.sample_name}")

    variant_df = load_variant_df(args.variant_hits)
    ase_df = load_ase_df(args.ase_hits)
    junction_df = load_junction_df(args.junction_hits)
    cohort_junction_df = load_cohort_junction_df(args.cohort_junction_hits) if args.cohort_junction_hits else None
    omim_df = load_omim_df(args.omim) if args.omim else None

    hit_df = build_hit_table(variant_df, ase_df, junction_df, cohort_junction_df, args.sample_name, omim_df)

    # Make output directory if it doesn't exist
    if args.outfile:
        outdir = os.path.dirname(args.outfile)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir)

    # Save hits
    hit_df.to_csv(args.outfile, sep='\t', index=False)
    print(f"Saved merged hits to {args.outfile}")

if __name__ == "__main__":
    main()
