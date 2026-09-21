

import argparse
import os
import re
import pandas as pd
import numpy as np
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

_DELTA_COL_TO_SUFFIX = {
    'delta_PSI':        'deltaPSI',
    'delta_PSI_approx': 'deltaPSIapprox',
    'delta_5ss_IR':     'delta5ssIR',
    'delta_3ss_IR':     'delta3ssIR',
    'delta_full_IR':    'deltaFullIR',
    'delta_IPA':        'deltaIPA',
}

_DELTA_SUFFIXES_BY_PREFIX = {
    '':        ('deltaPSI',),
    'cohort_': ('deltaPSI', 'deltaPSIapprox', 'delta5ssIR', 'delta3ssIR', 'deltaFullIR', 'deltaIPA'),
}


def parse_args():

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


_CALLER_PHASE_PRIORITY = ('nanoTS', 'clair3_rna', 'longcallR')


def _clnsig_category_rank(clnsig):
    s = str(clnsig)
    base = s.split(':', 1)[0]
    if base == 'Pathogenic':
        return 0
    if base == 'Pathogenic/Likely_pathogenic':
        return 1
    if base == 'Likely_pathogenic':
        return 2
    if base == 'Conflicting_classifications_of_pathogenicity':
        if re.search(r'Pathogenic|Likely_pathogenic', s):
            return 3
        return None
    if re.search(r'Pathogenic|Likely_pathogenic', s):
        return 4
    return None


def _caller_from_name(name):
    for caller in _CALLER_PHASE_PRIORITY + ('deepvariant',):
        if name.endswith('_' + caller):
            return caller
    return None


_CALLER_GT_DISPLAY = {
    'nanoTS':     'nanoTS',
    'longcallR':  'longcallR',
    'clair3_rna': 'clair3-RNA',
    'deepvariant': 'deepvariant',
}


def _build_gt_by_caller(variant_df):
    gt_long = (
        variant_df
            .groupby(['gene', 'variant_ID', '_caller'])['GT']
            .agg(lambda x: ','.join(x))
            .reset_index()
    )
    gt_wide = (
        gt_long
            .pivot(index=['gene', 'variant_ID'], columns='_caller', values='GT')
            .reset_index()
    )
    for caller in _CALLER_GT_DISPLAY:
        if caller not in gt_wide.columns:
            gt_wide[caller] = pd.NA
    gt_wide = gt_wide.rename(
        columns={caller: f'variant_GT_{suffix}' for caller, suffix in _CALLER_GT_DISPLAY.items()}
    )
    gt_cols = [f'variant_GT_{suffix}' for suffix in _CALLER_GT_DISPLAY.values()]
    gt_wide[gt_cols] = gt_wide[gt_cols].fillna('.')
    return gt_wide[['gene', 'variant_ID'] + gt_cols]


def resolve_phase(variant_a, variant_b):
    def _phased_calls(variant, caller):
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
                hap1_a = gt_a.split('|')[0]
                hap1_b = gt_b.split('|')[0]
                return 'cis' if hap1_a == hap1_b else 'trans'

    return 'unclear'


def load_variant_df(path):
    variant_df = pd.read_csv(path, sep='\t', usecols=[
        'sample', 'chrom', 'pos', 'ref', 'alt', 'name', 'GT', 'PS', 'gnomAD_AF', 'CLNSIG', 'gene', 'CADD_PHRED', 'SpliceAI',
        'num_callers', 'sample_count', 'ANNOVAR_AAChange.refGene', 'ANNOVAR_GeneDetail.refGene',
    ]).drop_duplicates()
    variant_df = variant_df.rename(columns={'sample_count': 'variant_nsamples'})
    variant_df = variant_df[variant_df['gene'] != '.']
    variant_df['gene'] = variant_df['gene'].str.split(',')
    variant_df = variant_df.explode('gene')
    variant_df['_caller'] = variant_df['name'].apply(_caller_from_name)
    variant_df['_raw_GT'] = variant_df['GT']
    variant_df['GT'] = variant_df['GT'].where(
        variant_df['GT'].str.contains('|', regex=False),
        variant_df['GT'].str.replace('1/0', '0/1', regex=False),
    )
    variant_df['variant_ID'] = variant_df.apply(lambda x: f"{x.chrom}-{x.pos}-{x.ref}-{x.alt}", axis=1)
    aachange = variant_df['ANNOVAR_AAChange.refGene']
    genedetail = variant_df['ANNOVAR_GeneDetail.refGene']
    variant_df['variant_consequence'] = aachange.where(aachange.notna() & (aachange != '.'), genedetail)
    variant_df['variant_consequence'] = variant_df['variant_consequence'].fillna('.')
    return variant_df


def load_ase_df(path):
    ase_df = pd.read_csv(path, sep='\t', usecols=['sample', 'gene', 'ratio', 'sample_count']).drop_duplicates()
    return ase_df.rename(columns={'ratio': 'ASE_ratio', 'sample_count': 'ASE_nsamples'})


def load_junction_df(path):
    return pd.read_csv(path, sep='\t', usecols=[
        'sample', 'gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event',
    ]).drop_duplicates()


def load_cohort_junction_df(path):
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
    omim_df['haploinsufficient'] = (
        omim_df['haploinsufficient'].astype(str).str.strip().str.upper() == 'TRUE'
    )
    return omim_df


def build_phased_junction_df(df, prefix, delta_cols=('delta_PSI',)):
    has_tissue = 'gtex_tissue' in df.columns
    tiers = {}
    for phasing, sep in (('bulk', ';'), ('hap1', ','), ('hap2', ',')):
        sub = df[df['phasing'] == phasing]

        if has_tissue:
            if len(sub):
                collapsed_rows = []
                for (gene, junction), g in sub.groupby(['gene', 'junction'], sort=False):
                    g = g.sort_values('gtex_tissue')
                    tissues = g['gtex_tissue'].astype(str)

                    if g['jxn_coverage'].astype(str).nunique() > 1:
                        print(f"WARNING: jxn_coverage disagrees across GTEx tissues for "
                              f"{gene} {junction} ({dict(zip(tissues, g['jxn_coverage']))}) -- using the first value.")
                    if g['annotation'].astype(str).nunique() > 1:
                        print(f"WARNING: annotation disagrees across GTEx tissues for "
                              f"{gene} {junction} ({dict(zip(tissues, g['annotation']))}) -- using the first value.")
                    row = {
                        'gene': gene,
                        'junction': junction,
                        'jxn_coverage': g['jxn_coverage'].iloc[0],
                        'annotation': g['annotation'].iloc[0],
                    }
                    for col in delta_cols:
                        row[col] = ' '.join(f"{v} ({t})" for v, t in zip(g[col], tissues))
                    row['event'] = ' '.join(f"{v} ({t})" for v, t in zip(g['event'], tissues))
                    row['sample_count'] = ' '.join(f"{v} ({t})" for v, t in zip(g['sample_count'], tissues))
                    collapsed_rows.append(row)
                sub = pd.DataFrame(collapsed_rows)
            else:
                sub = sub.drop(columns=['gtex_tissue'])

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
            sub.sort_values('junction')
                .groupby('gene')
                .agg(**agg_kwargs)
                .reset_index()
        )
    merged = pd.merge(tiers['bulk'], tiers['hap1'], on='gene', how='outer')
    merged = pd.merge(merged, tiers['hap2'], on='gene', how='outer')
    return merged.drop_duplicates()


def load_gene_expression_df(path):
    return pd.read_csv(path, sep='\t', index_col=0)


def load_gene_expression_zscore_df(path):
    return pd.read_csv(path, sep='\t', index_col=0)


def build_hit_table(variant_df, ase_df, junction_df, cohort_junction_df, sample_name, omim_df=None,
                     gene_expression_df=None, gene_expression_motr_df=None,
                     gene_expression_zscore_df=None, gene_expression_outlier_threshold=3.0):

    gt_by_caller = _build_gt_by_caller(variant_df)
    gt_caller_cols = [f'variant_GT_{suffix}' for suffix in _CALLER_GT_DISPLAY.values()]
    mod_variant_df = (
        variant_df
            .sort_values(['gene', 'variant_ID'])
            .groupby(['gene', 'variant_ID'], sort=False)
            .agg(
                gnomAD_AF=('gnomAD_AF', lambda x: ','.join(dict.fromkeys(x.dropna().astype(str)))),
                CLNSIG=('CLNSIG', lambda x: ','.join(dict.fromkeys(x.dropna().astype(str)))),
                CADD_PHRED=('CADD_PHRED', lambda x: ','.join(dict.fromkeys(x.dropna().astype(str)))),
                SpliceAI=('SpliceAI', lambda x: ','.join(dict.fromkeys(x.dropna().astype(str)))),
                variant_consequence=('variant_consequence', lambda x: ','.join(dict.fromkeys(x.dropna().astype(str)))),
                num_callers=('num_callers', lambda x: ','.join(dict.fromkeys(x.astype(str)))),
                variant_nsamples=('variant_nsamples', lambda x: ','.join(dict.fromkeys(x.astype(str)))),
            )
            .reset_index()
            .merge(gt_by_caller, on=['gene', 'variant_ID'], how='left')
            .groupby('gene', sort=False)
            .agg(
                variant_ID=('variant_ID', ';'.join),
                **{col: (col, ';'.join) for col in gt_caller_cols},
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
    mod_junction_df = build_phased_junction_df(junction_df, '')
    if cohort_junction_df is not None:
        mod_cohort_junction_df = build_phased_junction_df(
            cohort_junction_df, 'cohort_',
            delta_cols=('delta_PSI', 'delta_PSI_approx', 'delta_5ss_IR', 'delta_3ss_IR', 'delta_full_IR', 'delta_IPA'),
        )
    else:
        mod_cohort_junction_df = None

    hit_df = pd.merge(mod_variant_df, ase_df, on='gene', how='outer')
    hit_df = pd.merge(hit_df, mod_junction_df, on='gene', how='outer')
    if mod_cohort_junction_df is not None:
        hit_df = pd.merge(hit_df, mod_cohort_junction_df, on='gene', how='outer')
    else:
        for phasing in ('bulk', 'hap1', 'hap2'):
            for suffix in ('jxns', 'jxn_coverage', 'deltaPSI', 'deltaPSIapprox', 'delta5ssIR',
                           'delta3ssIR', 'deltaFullIR', 'deltaIPA', 'jxn_annotation', 'jxn_event', 'jxn_nsamples'):
                hit_df['cohort_' + phasing + '_' + suffix] = '.'
    if omim_df is not None:
        hit_df = pd.merge(hit_df, omim_df, on='gene', how='left')
        hit_df['haploinsufficient'] = hit_df['haploinsufficient'].fillna(False)
    else:
        hit_df['phenotypes'] = '.'
        hit_df['inheritance_patterns'] = '.'
        hit_df['haploinsufficient'] = False
    hit_df = hit_df.drop_duplicates()

    hit_df["variant"] = hit_df["variant_ID"].notna()
    hit_df["pathogenic_variant"] = (
        hit_df["variant_CLNSIG"]
            .astype(str)
            .str.contains(r"Pathogenic|Likely_pathogenic", regex=True)
    )
    hit_df["ASE"] = pd.to_numeric(hit_df["ASE_ratio"], errors="coerce").notna()


    _FLOAT_TOKEN_RE = re.compile(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?')

    def max_deltas(row, prefix=''):
        def parse_vals(s):
            if not pd.notna(s):
                return []
            vals = []
            for v in re.split('[,;]', s):
                v = v.strip()
                if not v:
                    continue
                for token in _FLOAT_TOKEN_RE.findall(v):
                    try:
                        vals.append(float(token))
                    except ValueError:
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
        hit_df['cohort_outlier_junction'] = hit_df.apply(lambda row: inspect_row(row, prefix='cohort_'), axis=1)
    else:
        hit_df['cohort_outlier_junction'] = 'None'

    hit_df = hit_df.astype(object)
    hit_df.fillna(".", inplace=True)

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
        empty = dict(
            n_pathogenic_variants=0, has_pathogenic=False, has_vus=False,
            has_homozygous_pathogenic=False, has_homozygous_vus=False,
            has_2_pathogenic_trans=False, has_2_pathogenic_unclear=False,
            has_1_pathogenic_1_vus_trans=False,
            has_2_vus_trans=False, has_2_vus_unclear=False,
            best_pathogenic_clnsig_rank=None,
        )
        if gene_variant_rows.empty:
            return empty

        variants = {}
        for variant_id, sub in gene_variant_rows.groupby('variant_ID'):
            clnsig_values = sub['CLNSIG'].astype(str)
            is_pathogenic = clnsig_values.str.contains(r'Pathogenic|Likely_pathogenic', regex=True).any()
            is_homozygous = sub['GT'].isin(['1/1', '1|1']).any()
            clnsig_rank = _clnsig_category_rank(clnsig_values.iloc[0]) if is_pathogenic else None
            variants[variant_id] = dict(
                is_pathogenic=is_pathogenic,
                is_homozygous=is_homozygous,
                clnsig_rank=clnsig_rank,
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
            best_pathogenic_clnsig_rank=(
                min((variants[v]['clnsig_rank'] for v in path_ids if variants[v]['clnsig_rank'] is not None), default=None)
            ),
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
        best_pathogenic_clnsig_rank=None,
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

        ase_hapi_or_strong = (ase and haploinsufficient) or strong
        ase_not_hapi_or_moderate = (ase and not haploinsufficient) or moderate

        if bucket == 1:
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

        if bucket == 2:
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
    hit_df['_tb_clnsig_rank'] = hit_df['gene'].map(
        lambda g: variant_class_by_gene.get(g, _empty_variant_class)['best_pathogenic_clnsig_rank']
    )

    def _max_bulk_delta(row):
        a = max_deltas(row, '')[0]
        b = max_deltas(row, 'cohort_')[0]
        vals = [v for v in (a, b) if pd.notna(v)]
        return max(vals) if vals else -1

    hit_df['_tb_max_bulk_delta'] = hit_df.apply(_max_bulk_delta, axis=1)

    def _max_cadd(row):
        s = row.get('variant_CADD_PHRED')
        if not pd.notna(s):
            return -1
        vals = []
        for v in str(s).split(';'):
            v = v.strip()
            if not v or v == '.':
                continue
            try:
                vals.append(float(v))
            except ValueError:
                continue
        return max(vals) if vals else -1

    hit_df['_tb_max_cadd'] = hit_df.apply(_max_cadd, axis=1)

    hit_df.sort_values(
        by=['tier', '_tb_n_pathogenic', '_tb_clnsig_rank', 'ASE', '_tb_max_bulk_delta', '_tb_max_cadd', 'gene'],
        ascending=[True, False, True, False, False, False, True],
        na_position='last',
        inplace=True
    )
    hit_df.drop(columns=['_tb_n_pathogenic', '_tb_clnsig_rank', '_tb_max_bulk_delta', '_tb_max_cadd'], inplace=True)
    hit_df['ranking'] = np.arange(1, len(hit_df) + 1)


    def _annotate_relative_expression(hit_df, expr_df, value_col, cohort_col, n_cohort_col=None):
        if expr_df is None:
            hit_df[value_col] = '.'
            hit_df[cohort_col] = '.'
            if n_cohort_col:
                hit_df[n_cohort_col] = '.'
            return

        if n_cohort_col:
            hit_df[n_cohort_col] = expr_df.shape[1]

        def _relative_expression(gene):
            if gene in expr_df.index and sample_name in expr_df.columns:
                value = expr_df.loc[gene, sample_name]
                if pd.isna(value):
                    return '.'
                return value
            return '.'

        hit_df[value_col] = hit_df['gene'].apply(_relative_expression)

        def _cohort_expression_summary(gene):
            if gene not in expr_df.index:
                return '.'
            values = pd.to_numeric(expr_df.loc[gene], errors='coerce').dropna()
            if values.empty:
                return '.'
            summary = np.percentile(values.to_numpy(), [0, 25, 50, 75, 100])
            return '[' + ','.join(str(v) for v in summary) + ']'

        hit_df[cohort_col] = hit_df['gene'].apply(_cohort_expression_summary)

    _annotate_relative_expression(
        hit_df, gene_expression_df,
        'relative_gene_expression', 'cohort_relative_gene_expression',
        n_cohort_col='n_cohort',
    )
    _annotate_relative_expression(
        hit_df, gene_expression_motr_df,
        'relative_gene_expression_motr', 'cohort_relative_gene_expression_motr',
    )

    if gene_expression_zscore_df is not None:
        def _expression_zscore(gene):
            if gene in gene_expression_zscore_df.index and sample_name in gene_expression_zscore_df.columns:
                z = gene_expression_zscore_df.loc[gene, sample_name]
                return z if pd.notna(z) else '.'
            return '.'

        def _expression_outlier(gene):
            z = _expression_zscore(gene)
            if z == '.':
                return '.'
            return bool(float(z) <= -abs(gene_expression_outlier_threshold))

        hit_df['gene_expression_zscore'] = hit_df['gene'].apply(_expression_zscore)
        hit_df['gene_expression_outlier'] = hit_df['gene'].apply(_expression_outlier)
    else:
        hit_df['gene_expression_zscore'] = '.'
        hit_df['gene_expression_outlier'] = '.'

    hit_df['sample'] = sample_name
    hit_df = hit_df[[
        'sample', 'gene', 'phenotypes', 'inheritance_patterns', 'haploinsufficient', 'ranking', 'tier', 'variant', 'pathogenic_variant', 'ASE', 'outlier_junction', 'cohort_outlier_junction',
        'relative_gene_expression', 'cohort_relative_gene_expression',
        'relative_gene_expression_motr', 'cohort_relative_gene_expression_motr',
        'gene_expression_zscore', 'gene_expression_outlier', 'n_cohort',
        'variant_ID', 'variant_GT_nanoTS', 'variant_GT_longcallR', 'variant_GT_clair3-RNA', 'variant_GT_deepvariant',
        'variant_gnomAD_AF',  'variant_CLNSIG', 'variant_CADD_PHRED', 'variant_SpliceAI', 'variant_consequence', 'variant_num_callers', 'variant_nsamples',
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

    if args.outfile:
        outdir = os.path.dirname(args.outfile)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir)

    hit_df.to_csv(args.outfile, sep='\t', index=False)
    print(f"Saved merged hits to {args.outfile}")

if __name__ == "__main__":
    main()