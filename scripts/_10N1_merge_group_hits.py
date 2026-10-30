# 10N1 - the per-sample gene hit table. pulls together variants, ASE, GTEx + cohort junction outliers,
# isoform outliers, OMIM, expression and AMELIE, then tiers and ranks the genes

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tier_rules
from transcript_outliers import QC_FAIL, QC_PASS
from gtf_utils import parse_gtf_junctions

warnings.filterwarnings("ignore", category=PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


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

_CALLER_PHASE_PRIORITY = ('wgs', 'nanoTS', 'clair3_rna', 'longcallR')


# genotype per caller, and whether two variants are in cis or trans
def _caller_from_name(name):
    for caller in _CALLER_PHASE_PRIORITY + ('deepvariant',):
        if name.endswith('_' + caller):
            return caller
    return None

_CALLER_GT_DISPLAY = {
    'wgs':        'WGS',
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


def _with_event_id(events, event_ids):
    return [f"{e} [{i}]" if isinstance(i, str) and i.strip() not in ('', '.', 'nan') and pd.notna(e) else e
            for e, i in zip(events, event_ids)]


def load_junction_df(path):
    cols = ['sample', 'gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event']
    df = pd.read_csv(path, sep='\t', usecols=lambda c: c in cols + ['event_id'])
    if 'event_id' in df.columns:
        df['event'] = _with_event_id(df['event'], df['event_id'].astype(str))
        df = df.drop(columns=['event_id'])
    return df[cols].drop_duplicates()


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

_ITEM_SEP = ';'

_ANNOT_SEP = ','


def _cell(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return '.'
    s = str(v).strip()
    return '.' if s in ('', 'nan', 'NaN', 'None') else s


def _join_items(values, sep=_ITEM_SEP):
    return sep.join(_cell(v) for v in values)


# one row per gene with bulk/hap1/hap2 columns side by side
def build_phased_junction_df(df, prefix, delta_cols=('delta_PSI',)):
    has_tissue = 'gtex_tissue' in df.columns
    tiers = {}
    df = df.copy()
    df['event'] = df['event'].apply(
        lambda e: _ANNOT_SEP.join(x.strip() for x in str(e).replace(';', ',').split(',') if x.strip())
        if pd.notna(e) else e
    )
    for phasing in ('bulk', 'hap1', 'hap2'):
        sep = _ITEM_SEP
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
                        row[col] = _ANNOT_SEP.join(f"{_cell(v)} ({t})" for v, t in zip(g[col], tissues))
                    row['event'] = _ANNOT_SEP.join(f"{_cell(v)} ({t})" for v, t in zip(g['event'], tissues))
                    row['sample_count'] = _ANNOT_SEP.join(f"{_cell(v)} ({t})" for v, t in zip(g['sample_count'], tissues))
                    collapsed_rows.append(row)
                sub = pd.DataFrame(collapsed_rows)
            else:
                sub = sub.drop(columns=['gtex_tissue'])

        agg_kwargs = {
            prefix + phasing + '_jxns':         ('junction', lambda x, sep=sep: _join_items(x, sep)),
            prefix + phasing + '_jxn_coverage': ('jxn_coverage', lambda x, sep=sep: _join_items(x, sep)),
        }
        for col in delta_cols:
            suffix = _DELTA_COL_TO_SUFFIX[col]
            agg_kwargs[prefix + phasing + '_' + suffix] = (col, lambda x, sep=sep: _join_items(x, sep))
        agg_kwargs[prefix + phasing + '_jxn_annotation'] = ('annotation', lambda x, sep=sep: _join_items(x, sep))
        agg_kwargs[prefix + phasing + '_jxn_event']      = ('event', lambda x, sep=sep: _join_items(x, sep))
        agg_kwargs[prefix + phasing + '_jxn_nsamples']   = ('sample_count', lambda x, sep=sep: _join_items(x, sep))
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


def build_hit_table(variant_df, ase_df, junction_df, cohort_junction_df, sample_name, omim_df=None,
                     gene_expression_df=None, gene_expression_motr_df=None, isoform_genes=None):

    gt_by_caller = _build_gt_by_caller(variant_df)
    gt_caller_cols = [f'variant_GT_{suffix}' for suffix in _CALLER_GT_DISPLAY.values()]
    # variants -> one row per gene. everything is joined in the same order so the columns still line up
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

    # outer join all the evidence on gene, then tack on OMIM
    hit_df = pd.merge(mod_variant_df, ase_df, on='gene', how='outer')
    hit_df = pd.merge(hit_df, mod_junction_df, on='gene', how='outer')
    if mod_cohort_junction_df is not None:
        hit_df = pd.merge(hit_df, mod_cohort_junction_df, on='gene', how='outer')
    else:
        for phasing in ('bulk', 'hap1', 'hap2'):
            for suffix in ('jxns', 'jxn_coverage', 'deltaPSI', 'deltaPSIapprox', 'delta5ssIR',
                           'delta3ssIR', 'deltaFullIR', 'deltaIPA', 'jxn_annotation', 'jxn_event'):
                hit_df['cohort_' + phasing + '_' + suffix] = '.'
    if isoform_genes:
        hit_df = pd.merge(hit_df, pd.DataFrame({'gene': sorted(isoform_genes), '_isoform': True}),
                          on='gene', how='outer')
    else:
        hit_df['_isoform'] = False
    hit_df['_isoform'] = hit_df['_isoform'].eq(True)
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

    # junction strength comes from the biggest delta-PSI in bulk or either haplotype
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
        return tier_rules.junction_strength(*max_deltas(row, prefix))
    hit_df['outlier_junction'] = hit_df.apply(inspect_row, axis=1)
    if mod_cohort_junction_df is not None:
        hit_df['cohort_outlier_junction'] = hit_df.apply(lambda row: inspect_row(row, prefix='cohort_'), axis=1)
    else:
        hit_df['cohort_outlier_junction'] = 'None'

    hit_df['RNA_dysregulation'] = (
        (hit_df['ASE'] == True)
        | (hit_df['outlier_junction'] != 'None')
        | (hit_df['cohort_outlier_junction'] != 'None')
        | hit_df['_isoform']
    )

    hit_df = hit_df.astype(object)
    hit_df.fillna(".", inplace=True)

    def _clean_inheritance(s):
        return str(s).strip().strip('"') if pd.notna(s) else ''

    # how many P/LP and VUS, homozygous, in trans, etc.
    def _classify_gene_variants(gene_variant_rows):
        empty = dict(
            n_pathogenic_variants=0, has_pathogenic=False, has_vus=False,
            has_homozygous_pathogenic=False, has_homozygous_vus=False,
            has_2_pathogenic_trans=False, has_2_pathogenic_unclear=False,
            has_1_pathogenic_1_vus_trans=False,
            has_2_vus_trans=False, has_2_vus_unclear=False, has_variants_in_trans=False,
        )
        if gene_variant_rows.empty:
            return empty

        variants = {}
        for variant_id, sub in gene_variant_rows.groupby('variant_ID'):
            clnsig_values = sub['CLNSIG'].astype(str)
            is_pathogenic = clnsig_values.str.contains(r'Pathogenic|Likely_pathogenic', regex=True).any()
            is_homozygous = sub['GT'].isin(['1/1', '1|1']).any()
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
            has_variants_in_trans=has_2_pathogenic_trans or has_1_pathogenic_1_vus_trans or has_2_vus_trans,
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
        has_2_vus_trans=False, has_2_vus_unclear=False, has_variants_in_trans=False,
    )

    # actual rules are in tier_rules.py
    def _assign_tier(row):
        vc = variant_class_by_gene.get(row['gene'], _empty_variant_class)
        return tier_rules.assign_tier(
            vc,
            ase=row['ASE'] is True,
            haploinsufficient=row['haploinsufficient'] is True,
            strength=tier_rules.strongest(row['outlier_junction'], row['cohort_outlier_junction']),
            dominant=tier_rules.is_dominant(_clean_inheritance(row['inheritance_patterns'])),
            isoform=bool(row['_isoform'] == True),
        )

    hit_df['tier'] = hit_df.apply(_assign_tier, axis=1)
    hit_df['variants_in_trans'] = hit_df['gene'].map(
        lambda g: variant_class_by_gene.get(g, _empty_variant_class)['has_variants_in_trans']
    )
    hit_df['_tb_n_pathogenic'] = hit_df['gene'].map(
        lambda g: variant_class_by_gene.get(g, _empty_variant_class)['n_pathogenic_variants']
    )
    hit_df['_tb_variants'] = hit_df['variant']

    # tie-breakers for ordering genes inside a tier
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
        by=['tier', '_tb_n_pathogenic', '_tb_variants', 'ASE', '_tb_max_bulk_delta', '_isoform', '_tb_max_cadd',
            'gene'],
        ascending=[True, False, False, False, False, False, False, True],
        na_position='last',
        inplace=True
    )
    hit_df.drop(columns=['_tb_n_pathogenic', '_tb_variants', '_tb_max_bulk_delta', '_tb_max_cadd'], inplace=True)
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
        hit_df, gene_expression_motr_df,
        'relative_gene_expression_motr', 'cohort_relative_gene_expression_motr',
        n_cohort_col='n_cohort',
    )

    hit_df['sample'] = sample_name
    hit_df = hit_df[[
        'sample', 'gene', 'phenotypes', 'inheritance_patterns', 'haploinsufficient', 'ranking', 'tier', 'variant', 'pathogenic_variant', 'ASE', 'outlier_junction', 'cohort_outlier_junction', 'RNA_dysregulation',
        'relative_gene_expression_motr', 'cohort_relative_gene_expression_motr',
        'n_cohort',
        'variant_ID', 'variant_GT_WGS', 'variant_GT_nanoTS', 'variant_GT_longcallR', 'variant_GT_clair3-RNA', 'variant_GT_deepvariant',
        'variant_gnomAD_AF',  'variant_CLNSIG', 'variant_CADD_PHRED', 'variant_SpliceAI', 'variant_consequence', 'variant_num_callers', 'variant_nsamples', 'variants_in_trans',
        'ASE_ratio', 'ASE_nsamples', 
        'bulk_jxns', 'bulk_jxn_coverage', 'bulk_deltaPSI', 'bulk_jxn_annotation', 'bulk_jxn_event', 'bulk_jxn_nsamples',
        'hap1_jxns', 'hap1_jxn_coverage', 'hap1_deltaPSI', 'hap1_jxn_annotation', 'hap1_jxn_event', 'hap1_jxn_nsamples',
        'hap2_jxns', 'hap2_jxn_coverage', 'hap2_deltaPSI', 'hap2_jxn_annotation', 'hap2_jxn_event', 'hap2_jxn_nsamples',
        'cohort_bulk_jxns', 'cohort_bulk_jxn_coverage',
        'cohort_bulk_deltaPSI', 'cohort_bulk_deltaPSIapprox', 'cohort_bulk_delta5ssIR', 'cohort_bulk_delta3ssIR', 'cohort_bulk_deltaFullIR', 'cohort_bulk_deltaIPA',
        'cohort_bulk_jxn_annotation', 'cohort_bulk_jxn_event',
        'cohort_hap1_jxns', 'cohort_hap1_jxn_coverage',
        'cohort_hap1_deltaPSI', 'cohort_hap1_deltaPSIapprox', 'cohort_hap1_delta5ssIR', 'cohort_hap1_delta3ssIR', 'cohort_hap1_deltaFullIR', 'cohort_hap1_deltaIPA',
        'cohort_hap1_jxn_annotation', 'cohort_hap1_jxn_event',
        'cohort_hap2_jxns', 'cohort_hap2_jxn_coverage',
        'cohort_hap2_deltaPSI', 'cohort_hap2_deltaPSIapprox', 'cohort_hap2_delta5ssIR', 'cohort_hap2_delta3ssIR', 'cohort_hap2_deltaFullIR', 'cohort_hap2_deltaIPA',
        'cohort_hap2_jxn_annotation', 'cohort_hap2_jxn_event'
    ]]
    return hit_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build every sample's ranked candidate-hits table for a group, directly from the "
                    "group's already-merged variant/ASE/junction tables, and write one concatenated "
                    "all_hits.tsv -- replaces split_group_hits_by_sample.py + per-sample merge_hits.py + "
                    "the old awk concat.")
    parser.add_argument("--outfile", required=True, help="Path to write the group-level all_hits.tsv to.")
    parser.add_argument("--variant-tsv", required=True,
        help="Group-level merged variant hits (output of _10K_merge_group_variants.py).")
    parser.add_argument("--ase-tsv", required=True,
        help="Group-level merged ASE hits (output of _10L_merge_group_ase.py).")
    parser.add_argument("--tissues", nargs="*", default=[],
        help="Tissues present in this group (same order as --junction-files).")
    parser.add_argument("--junction-files", nargs="*", default=[],
        help="Per-tissue merged junction hit files (output of _10M_merge_group_junctions.py, "
             "one per --tissues entry). Concatenated across tissues before filtering by sample.")
    parser.add_argument("--cohort-junction-tsv", required=False, default=None,
        help="Group-level cohort-comparison junction outliers (rules/7_cohort_junction_analysis.smk's "
             "*_outliers_filtered.tsv). Raw column names differ from --junction-files (and vary by "
             "whether this group used the beta_binomial or modified_zscore method) -- normalized to "
             "merge_hits.load_cohort_junction_df's schema here. If omitted/missing, every sample's "
             "cohort_* columns are filled with '.' (merge_hits.build_hit_table's existing fallback).")
    parser.add_argument("--drop-cohort-events", nargs="*", default=[],
        help="Cohort event labels to remove before merging (e.g. 5ss_IR 3ss_IR). A junction left with no event "
             "is dropped, so those calls add no evidence and play no part in tiering.")
    parser.add_argument("--samples", nargs="+", required=True, help="Every sample name in this group.")
    parser.add_argument("--omim", required=False, default=None,
        help="Path to OMIM data. If omitted, phenotypes/inheritance_patterns/haploinsufficient are filled with '.'/False.")
    parser.add_argument("--gene-expression-matrix", required=False, default=None,
        help="Group-level targeted-panel CPTM matrix (rule _10H's <outprefix>_matrix.tsv, from "
             "_10_merge_gene_expression_by_assignment.py) -- one row per gene, one column per sample in this "
             "(bed_id, sample_type) group. If omitted/missing, relative_gene_expression/n_cohort "
             "are filled with '.' (merge_hits.build_hit_table's existing fallback convention for "
             "optional inputs).")
    parser.add_argument("--gene-expression-matrix-motr", required=False, default=None,
        help="Group-level targeted-panel MOTR matrix (rule _10H's <outprefix>_matrix_motr.tsv, "
             "from _10_merge_gene_expression_by_assignment.py's compute_size_factors() -- DESeq2-style "
             "median-of-ratios normalization restricted to BED-panel genes). Same shape as "
             "--gene-expression-matrix. If omitted/missing, relative_gene_expression_motr/"
             "cohort_relative_gene_expression_motr are filled with '.'.")
    parser.add_argument("--debug-sample", required=False, default=None,
        help="If set, print diagnostic detail (to stderr) for this one sample: the gene list surviving "
             "the per-sample variant_df slice, and the gene list in build_hit_table's output -- useful "
             "for tracking down a gene that unexpectedly disappears for one specific sample. No effect "
             "on output files. Off by default.")
    parser.add_argument("--debug-gene", required=False, default=None,
        help="Used together with --debug-sample: also print the matching variant_df row(s) for this "
             "gene (if any) within that sample's slice, and whether the gene made it into "
             "build_hit_table's output. Ignored if --debug-sample isn't also set.")
    parser.add_argument("--bed", required=False, default=None,
        help="Panel BED (gene symbol in column 4). Every gene in it gets a row for every sample, "
             "with hit=FALSE and '.' in the hit columns when the gene has no hit in that sample. "
             "If omitted, only genes with hits are written (hit=TRUE for all rows).")
    parser.add_argument("--gtf", required=False, default=None,
        help="Reference annotation GTF, used to get each gene's canonical (Ensembl_canonical) "
             "junctions for n_canonical_jxns / n_canonical_jxns_tested.")
    parser.add_argument("--avg-coverage-files", nargs="*", default=[],
        help="Each sample's {sample}_avg_coverage.tsv (rule _6D), in the same order as --samples. "
             "Its average_coverage / breadth_coverage columns become the hit table's "
             "average_canonical_coverage / breadth_canonical_coverage.")
    parser.add_argument("--gtex-all-junction-files", nargs="*", default=[],
        help="Every sample's per-tissue GTEx {sample}_gtex_{tissue}_all_junctions.tsv (rule _5C). "
             "A canonical junction counts as tested in GTEx if its bulk row has a numeric padj.")
    parser.add_argument("--cohort-results-dir", required=False, default=None,
        help="This group's cohort scoring results directory ({group}_results_{method}/ from rule "
             "_8B), holding one {gene}.tsv per gene. A canonical junction counts as tested in the "
             "cohort if its bulk row has a numeric junction_PSI_approx test value.")
    parser.add_argument("--amelie-files", nargs="*", default=[],
        help="Per-sample {sample}_amelie.tsv files from rule _6E (samples with HPO terms). Adds AMELIE_rank, "
             "AMELIE_score, AMELIE_pmids and AMELIE_n_genes for every panel gene; '.' for samples without HPO "
             "terms, and AMELIE_rank 'error' where the AMELIE query failed.")
    parser.add_argument("--isoform-outlier-files", nargs="*", default=[],
        help="This group's transcript usage outlier files from rule _9G (individual isoform outliers and "
             "isoform mix outliers). Every sample/gene in either file is a hit and counts as RNA dysregulation, "
             "with the same weight as a weak outlier splicing event when tiering; isoform_outlier is TRUE for "
             "those genes and FALSE otherwise. Calls from samples that failed 9F's isoform sample QC (qc column "
             "'fail') are not used for tiering and are marked isoform_outlier 'fail'; isoform_qc gives each "
             "sample's QC result. Without these files both columns are '.'.")
    args = parser.parse_args()
    if len(args.tissues) != len(args.junction_files):
        parser.error("--tissues and --junction-files must have the same number of entries")
    if args.avg_coverage_files and len(args.avg_coverage_files) != len(args.samples):
        parser.error("--avg-coverage-files must have one entry per --samples entry")
    return args


_TESTABILITY_COLS = ['average_canonical_coverage', 'breadth_canonical_coverage', 'n_canonical_jxns', 'n_canonical_jxns_tested']
_LEAD_COLS = ['sample', 'gene', 'phenotypes', 'inheritance_patterns', 'haploinsufficient']


def _load_bed_genes(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 4 and fields[3] not in genes:
                genes.append(fields[3])
    return genes


def _load_coverage(samples, files):
    out = {}
    for sample, path in zip(samples, files):
        if not os.path.isfile(path):
            print(f"WARNING: coverage file missing for {sample}: {path}")
            continue
        df = pd.read_csv(path, sep='\t')
        for col in ('average_coverage', 'breadth_coverage'):
            if col not in df.columns:
                df[col] = float('nan')
        for gene, avg, breadth in zip(df['gene'], df['average_coverage'], df['breadth_coverage']):
            out[(str(sample), str(gene))] = (avg, breadth)
    return out


_AMELIE_COLS = ['AMELIE_rank', 'AMELIE_score', 'AMELIE_pmids', 'AMELIE_n_genes']


def _load_amelie(files):
    out = {}
    for path in files:
        if not os.path.isfile(path):
            print(f"WARNING: AMELIE file missing: {path}")
            continue
        df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        for row in df.itertuples(index=False):
            if row.AMELIE_status == 'error':
                values = ('error', '.', '.', row.AMELIE_n_genes)
            else:
                values = (row.AMELIE_rank, row.AMELIE_score, row.AMELIE_pmids, row.AMELIE_n_genes)
            out[(str(row.sample), str(row.gene))] = values
    return out


def _load_isoform_outliers(files):
    pairs, failed_pairs, failed = set(), set(), set()
    for path in files:
        if not os.path.isfile(path):
            print(f"WARNING: isoform outlier file missing: {path}")
            continue
        df = pd.read_csv(path, sep='\t', dtype=str)
        qc = df['qc'].astype(str).str.lower() if 'qc' in df.columns else pd.Series(QC_PASS, index=df.index)
        bad = qc.eq(QC_FAIL)
        pairs.update(zip(df.loc[~bad, 'sample'].astype(str), df.loc[~bad, 'gene'].astype(str)))
        failed_pairs.update(zip(df.loc[bad, 'sample'].astype(str), df.loc[bad, 'gene'].astype(str)))
        failed.update(df.loc[bad, 'sample'].astype(str))
    if failed:
        print(f"Isoform sample QC: outliers of {len(failed)} failed sample(s) are not used for tiering: "
              + ", ".join(sorted(failed)))
    return pairs, failed_pairs, failed


# respects the isoform sample QC, failed samples don't get isoform hits
def add_isoform_outlier(hit_df, sample, pairs, failed_pairs, failed, available):
    at = list(hit_df.columns).index('hit') + 1 if 'hit' in hit_df.columns else len(hit_df.columns)
    def value(g):
        if not available:
            return '.'
        if (str(sample), str(g)) in pairs:
            return 'TRUE'
        return QC_FAIL if (str(sample), str(g)) in failed_pairs else 'FALSE'
    hit_df.insert(at, 'isoform_outlier', [value(g) for g in hit_df['gene']])
    qc = (QC_FAIL if str(sample) in failed else QC_PASS) if available else '.'
    hit_df.insert(at + 1, 'isoform_qc', qc)
    return hit_df


def add_amelie(hit_df, sample, amelie):
    values = [amelie.get((str(sample), str(g)), ('.', '.', '.', '.')) for g in hit_df['gene']]
    at = list(hit_df.columns).index('tier') + 1 if 'tier' in hit_df.columns else len(hit_df.columns)
    for k, col in enumerate(_AMELIE_COLS):
        hit_df.insert(at + k, col, [v[k] for v in values])
    return hit_df


def _numeric_mask(series):
    return pd.to_numeric(series, errors='coerce').notna()


# canonical junctions GTEx or the cohort analysis actually managed to test
def _tested_canonical_junctions(canonical, samples, gtex_files, cohort_dir):
    canonical_all = set().union(*canonical.values()) if canonical else set()
    samples = {str(s) for s in samples}
    tested = {s: set() for s in samples}

    def _add(df, test_col):
        df = df[(df['phasing'] == 'bulk') & df['junction'].isin(canonical_all)
                & df['sample'].astype(str).isin(samples) & _numeric_mask(df[test_col])]
        for sample, jxn in zip(df['sample'].astype(str), df['junction']):
            tested[sample].add(jxn)

    for path in gtex_files:
        if not os.path.isfile(path):
            print(f"WARNING: GTEx all-junctions file missing: {path}")
            continue
        df = pd.read_csv(path, sep='\t', dtype=str, usecols=['sample', 'phasing', 'junction', 'padj'])
        _add(df, 'padj')

    if cohort_dir and os.path.isdir(cohort_dir):
        want = ('p_value_junction_PSI_approx', 'modz_junction_PSI_approx')
        for gene in canonical:
            path = os.path.join(cohort_dir, f"{gene}.tsv")
            if not os.path.isfile(path):
                continue
            df = pd.read_csv(path, sep='\t', dtype=str,
                             usecols=lambda c: c in ('sample', 'phasing', 'junction') + want)
            test_col = next((c for c in want if c in df.columns), None)
            if test_col is not None:
                _add(df, test_col)
    elif cohort_dir:
        print(f"WARNING: cohort results directory not found: {cohort_dir}")
    return tested


def _fmt(v):
    return '.' if v is None or (isinstance(v, float) and pd.isna(v)) else v


# panel genes without hits still go in, so you can see whether they were even testable
def add_testability_and_non_hits(hit_df, sample, panel_genes, omim_df, coverage, canonical, tested):
    hit_df = hit_df.copy()
    hit_df['hit'] = 'TRUE'
    hit_genes = set(hit_df['gene'].astype(str))
    missing = [g for g in (panel_genes or []) if g not in hit_genes]
    if missing:
        filler = pd.DataFrame('.', index=range(len(missing)), columns=hit_df.columns)
        filler['sample'] = sample
        filler['gene'] = missing
        filler['hit'] = 'FALSE'
        filler['haploinsufficient'] = False
        if omim_df is not None:
            omim = omim_df.drop_duplicates('gene').set_index('gene')
            for col in ('phenotypes', 'inheritance_patterns'):
                filler[col] = [_fmt(omim[col].get(g)) if g in omim.index else '.' for g in missing]
            filler['haploinsufficient'] = [bool(omim['haploinsufficient'].get(g)) if g in omim.index else False
                                           for g in missing]
        hit_df = pd.concat([hit_df, filler], ignore_index=True)

    def _cov(gene, i):
        v = coverage.get((str(sample), str(gene)))
        return _fmt(v[i]) if v is not None else '.'

    hit_df['average_canonical_coverage'] = [_cov(g, 0) for g in hit_df['gene']]
    hit_df['breadth_canonical_coverage'] = [_cov(g, 1) for g in hit_df['gene']]
    done = tested.get(str(sample), set())
    hit_df['n_canonical_jxns'] = [len(canonical[g]) if g in canonical else '.' for g in hit_df['gene']]
    hit_df['n_canonical_jxns_tested'] = [len(canonical[g] & done) if g in canonical else '.'
                                         for g in hit_df['gene']]

    rest = [c for c in hit_df.columns if c not in _LEAD_COLS + ['hit'] + _TESTABILITY_COLS]
    return hit_df[_LEAD_COLS + ['hit'] + _TESTABILITY_COLS + rest]


# reshape the 8C outliers so they look like the GTEx junction table
def normalize_cohort_junction_df(df):
    out = pd.DataFrame()
    out['sample'] = df['sample']
    out['phasing'] = df['phasing']
    out['gene'] = df['gene']
    out['junction'] = df['junction']
    out['jxn_coverage'] = df['junction_coverage'] if 'junction_coverage' in df.columns else ''

    _METRIC_EVENTS = {
        'junction_PSI_approx':    {'alt_5ss_approx', 'alt_3ss_approx', 'exon_skipping_approx', 'exon_inclusion_approx', 'complex_approx'},
        'junction_PSI':           {'alt_5ss', 'alt_3ss', 'exon_skipping', 'exon_inclusion', 'complex'},
        '5ss_IR_ratio':           {'5ss_IR'},
        '3ss_IR_ratio':           {'3ss_IR'},
        'junction_full_IR_ratio': {'full_IR'},
        'junction_IPA_ratio':     {'IPA'},
    }
    _METRIC_TO_OUTCOL = {
        'junction_PSI':           'delta_PSI',
        'junction_PSI_approx':    'delta_PSI_approx',
        '5ss_IR_ratio':           'delta_5ss_IR',
        '3ss_IR_ratio':           'delta_3ss_IR',
        'junction_full_IR_ratio': 'delta_full_IR',
        'junction_IPA_ratio':     'delta_IPA',
    }

    event_type = df['event_type'] if 'event_type' in df.columns else pd.Series([''] * len(df), index=df.index)
    fired_events = event_type.fillna('').apply(lambda s: set(str(s).split(',')) if s else set())

    for metric, outcol in _METRIC_TO_OUTCOL.items():
        if f'delta_{metric}' in df.columns:
            raw = df[f'delta_{metric}']
        elif f'modz_{metric}' in df.columns:
            raw = df[f'modz_{metric}']
        else:
            out[outcol] = ''
            continue
        is_outlier_for_metric = fired_events.apply(lambda evs, m=metric: bool(evs & _METRIC_EVENTS[m]))
        out[outcol] = raw.where(is_outlier_for_metric, other='')

    out['annotation'] = df['junction_type'] if 'junction_type' in df.columns else '.'
    out['event'] = df['event_type'] if 'event_type' in df.columns else '.'
    if 'event_id' in df.columns and 'event_type' in df.columns:
        out['event'] = _with_event_id(df['event_type'], df['event_id'].astype(str))
    out['sample_count'] = df.groupby(['gene', 'junction'])['sample'].transform('nunique')
    return out


# e.g. the default filter set throws out the 5'/3' IR events
def drop_cohort_events(df, drop):
    def keep(text):
        labels = [e.strip() for e in str(text).split(',') if e.strip() and e.strip() not in drop]
        return ','.join(labels) if labels else 'none'
    n_before = len(df)
    df = df.assign(event_type=df['event_type'].fillna('').apply(keep))
    df = df[df['event_type'] != 'none'].reset_index(drop=True)
    print(f"Removed cohort event label(s) {', '.join(sorted(drop))}: {n_before - len(df)} of {n_before} cohort "
          f"junction row(s) had no other event and were dropped")
    return df


_EMPTY_JUNCTION_COLS = ['sample', 'gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event']


def main():
    args = parse_args()

    print(f"\nBuilding group-level hits for {len(args.samples)} sample(s)...")

    variant_df = load_variant_df(args.variant_tsv)
    ase_df = load_ase_df(args.ase_tsv)

    print("Reading junction hits (concatenated across tissues)...")
    per_tissue_dfs = []
    for tissue, jxn_file in zip(args.tissues, args.junction_files):
        if os.path.isfile(jxn_file):
            tissue_df = load_junction_df(jxn_file)
            tissue_df['gtex_tissue'] = tissue
            per_tissue_dfs.append(tissue_df)
        else:
            print(f"WARNING: No junction hits file found for tissue {tissue} at {jxn_file}")
    junction_df = (pd.concat(per_tissue_dfs, ignore_index=True) if per_tissue_dfs
                   else pd.DataFrame(columns=_EMPTY_JUNCTION_COLS))

    print("Reading cohort-comparison junction hits...")
    cohort_junction_df = None
    if args.cohort_junction_tsv and os.path.isfile(args.cohort_junction_tsv):
        _required_cols = {'sample', 'phasing', 'gene', 'junction'}
        raw_cohort_df = pd.read_csv(args.cohort_junction_tsv, sep='\t', dtype=str)
        if args.drop_cohort_events and 'event_type' in raw_cohort_df.columns:
            raw_cohort_df = drop_cohort_events(raw_cohort_df, set(args.drop_cohort_events))
        if _required_cols.issubset(raw_cohort_df.columns) and not raw_cohort_df.empty:
            cohort_junction_df = normalize_cohort_junction_df(raw_cohort_df)

    omim_df = load_omim_df(args.omim) if args.omim else None

    gene_expression_df = None
    if args.gene_expression_matrix and os.path.isfile(args.gene_expression_matrix):
        gene_expression_df = load_gene_expression_df(args.gene_expression_matrix)

    gene_expression_motr_df = None
    if args.gene_expression_matrix_motr and os.path.isfile(args.gene_expression_matrix_motr):
        gene_expression_motr_df = load_gene_expression_df(args.gene_expression_matrix_motr)

    panel_genes = _load_bed_genes(args.bed) if args.bed else None
    coverage = _load_coverage(args.samples, args.avg_coverage_files) if args.avg_coverage_files else {}
    isoform_pairs, isoform_failed_pairs, isoform_failed = _load_isoform_outliers(args.isoform_outlier_files)
    # one table per sample
    hit_dfs = []
    for sample in args.samples:
        sample_variant_df = variant_df[variant_df['sample'].astype(str) == str(sample)]
        sample_ase_df = ase_df[ase_df['sample'].astype(str) == str(sample)]
        sample_junction_df = junction_df[junction_df['sample'].astype(str) == str(sample)]
        sample_cohort_junction_df = (
            cohort_junction_df[cohort_junction_df['sample'].astype(str) == str(sample)]
            if cohort_junction_df is not None else None
        )

        if args.debug_sample and str(sample) == args.debug_sample:
            genes_in_slice = sorted(sample_variant_df['gene'].dropna().unique().tolist())
            print(f"  [DEBUG] {sample}: {len(sample_variant_df)} variant_df row(s) after per-sample slice, "
                  f"genes: {genes_in_slice}", file=sys.stderr)
            if args.debug_gene:
                match = sample_variant_df[sample_variant_df['gene'] == args.debug_gene]
                print(f"  [DEBUG] {sample}: {len(match)} row(s) for gene=={args.debug_gene!r} in the slice:",
                      file=sys.stderr)
                if not match.empty:
                    print(match.to_string(), file=sys.stderr)

        hit_df = build_hit_table(
            sample_variant_df, sample_ase_df, sample_junction_df, sample_cohort_junction_df,
            sample, omim_df, gene_expression_df,
            gene_expression_motr_df=gene_expression_motr_df,
            isoform_genes={g for s, g in isoform_pairs if s == str(sample)},
        )

        if args.debug_sample and str(sample) == args.debug_sample:
            genes_in_output = sorted(hit_df['gene'].dropna().unique().tolist())
            print(f"  [DEBUG] {sample}: {len(hit_df)} row(s) in build_hit_table's output, "
                  f"genes: {genes_in_output}", file=sys.stderr)
            if args.debug_gene:
                present = args.debug_gene in genes_in_output
                print(f"  [DEBUG] {sample}: gene=={args.debug_gene!r} present in output: {present}", file=sys.stderr)

        hit_dfs.append((sample, hit_df))

    # testability, the non-hit genes and AMELIE ranks
    canonical = {}
    if args.gtf:
        genes_for_gtf = list(dict.fromkeys(
            list(panel_genes or []) + [str(g) for _, d in hit_dfs for g in d['gene']]))
        print(f"Parsing canonical junctions for {len(genes_for_gtf)} gene(s) from {args.gtf} ...")
        canonical = {g: d['canonical_junctions']
                     for g, d in parse_gtf_junctions(args.gtf, genes_for_gtf).items()}
    print("Collecting tested canonical junctions (GTEx + cohort)...")
    tested = _tested_canonical_junctions(canonical, args.samples, args.gtex_all_junction_files,
                                         args.cohort_results_dir)

    amelie = _load_amelie(args.amelie_files)
    out_dfs = []
    for sample, hit_df in hit_dfs:
        n_hits = len(hit_df)
        hit_df = add_testability_and_non_hits(hit_df, sample, panel_genes, omim_df,
                                              coverage, canonical, tested)
        hit_df = add_amelie(hit_df, sample, amelie)
        hit_df = add_isoform_outlier(hit_df, sample, isoform_pairs, isoform_failed_pairs, isoform_failed,
                                     bool(args.isoform_outlier_files))
        out_dfs.append(hit_df)
        print(f"  {sample}: {n_hits} gene(s) with hits, {len(hit_df)} row(s) total")
    hit_dfs = out_dfs

    all_hits = pd.concat(hit_dfs, ignore_index=True) if hit_dfs else pd.DataFrame()

    outdir = os.path.dirname(args.outfile)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    all_hits.to_csv(args.outfile, sep='\t', index=False)
    print(f"\nSaved group-level hits ({len(all_hits)} total row(s)) to {args.outfile}")


if __name__ == "__main__":
    main()