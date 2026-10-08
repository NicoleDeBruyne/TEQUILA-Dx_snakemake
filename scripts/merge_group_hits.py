

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_hits
from identify_cohort_junction_outliers import parse_gtf_junctions as _parse_gtf_junctions


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build every sample's ranked candidate-hits table for a group, directly from the "
                    "group's already-merged variant/ASE/junction tables, and write one concatenated "
                    "all_hits.tsv -- replaces split_group_hits_by_sample.py + per-sample merge_hits.py + "
                    "the old awk concat.")
    parser.add_argument("--outfile", required=True, help="Path to write the group-level all_hits.tsv to.")
    parser.add_argument("--variant-tsv", required=True,
        help="Group-level merged variant hits (output of merge_and_filter_variants.py).")
    parser.add_argument("--ase-tsv", required=True,
        help="Group-level merged ASE hits (output of merge_and_filter_ase_results.py).")
    parser.add_argument("--tissues", nargs="*", default=[],
        help="Tissues present in this group (same order as --junction-files).")
    parser.add_argument("--junction-files", nargs="*", default=[],
        help="Per-tissue merged junction hit files (output of merge_and_filter_junction_results.py, "
             "one per --tissues entry). Concatenated across tissues before filtering by sample.")
    parser.add_argument("--cohort-junction-tsv", required=False, default=None,
        help="Group-level cohort-comparison junction outliers (rules/7_cohort_junction_analysis.smk's "
             "*_outliers_filtered.tsv). Raw column names differ from --junction-files (and vary by "
             "whether this group used the beta_binomial or modified_zscore method) -- normalized to "
             "merge_hits.load_cohort_junction_df's schema here. If omitted/missing, every sample's "
             "cohort_* columns are filled with '.' (merge_hits.build_hit_table's existing fallback).")
    parser.add_argument("--samples", nargs="+", required=True, help="Every sample name in this group.")
    parser.add_argument("--omim", required=False, default=None,
        help="Path to OMIM data. If omitted, phenotypes/inheritance_patterns/haploinsufficient are filled with '.'/False.")
    parser.add_argument("--gene-expression-matrix", required=False, default=None,
        help="Group-level targeted-panel CPTM matrix (rule _9H1's <outprefix>_matrix.tsv, from "
             "quantify_gene_by_assignment.py) -- one row per gene, one column per sample in this "
             "(bed_id, sample_type) group. If omitted/missing, relative_gene_expression/n_cohort "
             "are filled with '.' (merge_hits.build_hit_table's existing fallback convention for "
             "optional inputs).")
    parser.add_argument("--gene-expression-matrix-motr", required=False, default=None,
        help="Group-level targeted-panel MOTR matrix (rule _9H1's <outprefix>_matrix_motr.tsv, "
             "from quantify_gene_by_assignment.py's compute_size_factors() -- DESeq2-style "
             "median-of-ratios normalization restricted to BED-panel genes). Same shape as "
             "--gene-expression-matrix. If omitted/missing, relative_gene_expression_motr/"
             "cohort_relative_gene_expression_motr are filled with '.'.")
    parser.add_argument("--gene-expression-zscores", required=False, default=None,
        help="Group-level MOTR-normalized low-expression z-score matrix (rule _9H1's "
             "<outprefix>_zscores_motr.tsv, from quantify_gene_by_assignment.py -- see "
             "scripts/expression_outliers.py's module docstring for the algorithm). Same shape "
             "as --gene-expression-matrix. If omitted/missing, gene_expression_zscore_motr is "
             "filled with '.'. Annotation only -- does NOT currently factor into a gene's tier.")
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
    args = parser.parse_args()
    if len(args.tissues) != len(args.junction_files):
        parser.error("--tissues and --junction-files must have the same number of entries")
    if args.avg_coverage_files and len(args.avg_coverage_files) != len(args.samples):
        parser.error("--avg-coverage-files must have one entry per --samples entry")
    return args


# --------------------------------------------------------------------------------------------
# Gene x sample testability: coverage + canonical-junction testing, reported for every gene
# --------------------------------------------------------------------------------------------
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
    """{(sample, gene): (average_coverage, breadth_coverage)}"""
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


def _numeric_mask(series):
    return pd.to_numeric(series, errors='coerce').notna()


def _tested_canonical_junctions(canonical, samples, gtex_files, cohort_dir):
    """{sample: set(junction)} of canonical junctions whose bulk PSI_approx was tested in GTEx or the cohort."""
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


def add_testability_and_non_hits(hit_df, sample, panel_genes, omim_df, coverage, canonical, tested):
    """Add hit + the four testability columns to one sample's hit table, and (if panel_genes is
    given) one '.'-filled row per panel gene that has no hit. Hits keep their ranking order and
    come first; non-hit rows follow, sorted by gene."""
    hit_df = hit_df.copy()
    hit_df['hit'] = 'TRUE'
    hit_genes = set(hit_df['gene'].astype(str))
    missing = [g for g in (panel_genes or []) if g not in hit_genes]
    if missing:
        filler = pd.DataFrame('.', index=range(len(missing)), columns=hit_df.columns)
        filler['sample'] = sample
        filler['gene'] = missing
        filler['hit'] = 'FALSE'
        # Gene-level OMIM annotations, filled the same way as for hit rows (merge_hits.build_hit_table):
        # phenotypes/inheritance_patterns '.' when the gene isn't in OMIM; haploinsufficient False.
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
        out['event'] = merge_hits._with_event_id(df['event_type'], df['event_id'].astype(str))
    out['sample_count'] = df['n_sample_outlier_junction_PSI'] if 'n_sample_outlier_junction_PSI' in df.columns else ''
    return out


_EMPTY_JUNCTION_COLS = ['sample', 'gene', 'phasing', 'junction', 'jxn_coverage', 'delta_PSI', 'sample_count', 'annotation', 'event']


def main():
    args = parse_args()

    print(f"\nBuilding group-level hits for {len(args.samples)} sample(s)...")

    variant_df = merge_hits.load_variant_df(args.variant_tsv)
    ase_df = merge_hits.load_ase_df(args.ase_tsv)

    print("Reading junction hits (concatenated across tissues)...")
    per_tissue_dfs = []
    for tissue, jxn_file in zip(args.tissues, args.junction_files):
        if os.path.isfile(jxn_file):
            tissue_df = merge_hits.load_junction_df(jxn_file)
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
        if _required_cols.issubset(raw_cohort_df.columns) and not raw_cohort_df.empty:
            cohort_junction_df = normalize_cohort_junction_df(raw_cohort_df)

    omim_df = merge_hits.load_omim_df(args.omim) if args.omim else None

    gene_expression_df = None
    if args.gene_expression_matrix and os.path.isfile(args.gene_expression_matrix):
        gene_expression_df = merge_hits.load_gene_expression_df(args.gene_expression_matrix)

    gene_expression_motr_df = None
    if args.gene_expression_matrix_motr and os.path.isfile(args.gene_expression_matrix_motr):
        gene_expression_motr_df = merge_hits.load_gene_expression_df(args.gene_expression_matrix_motr)

    gene_expression_zscore_df = None
    if args.gene_expression_zscores and os.path.isfile(args.gene_expression_zscores):
        gene_expression_zscore_df = merge_hits.load_gene_expression_zscore_df(args.gene_expression_zscores)

    panel_genes = _load_bed_genes(args.bed) if args.bed else None
    coverage = _load_coverage(args.samples, args.avg_coverage_files) if args.avg_coverage_files else {}
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

        hit_df = merge_hits.build_hit_table(
            sample_variant_df, sample_ase_df, sample_junction_df, sample_cohort_junction_df,
            sample, omim_df, gene_expression_df,
            gene_expression_motr_df=gene_expression_motr_df,
            gene_expression_zscore_df=gene_expression_zscore_df,
        )

        if args.debug_sample and str(sample) == args.debug_sample:
            genes_in_output = sorted(hit_df['gene'].dropna().unique().tolist())
            print(f"  [DEBUG] {sample}: {len(hit_df)} row(s) in build_hit_table's output, "
                  f"genes: {genes_in_output}", file=sys.stderr)
            if args.debug_gene:
                present = args.debug_gene in genes_in_output
                print(f"  [DEBUG] {sample}: gene=={args.debug_gene!r} present in output: {present}", file=sys.stderr)

        hit_dfs.append((sample, hit_df))

    # Canonical junctions for every panel gene plus any hit gene outside the panel
    canonical = {}
    if args.gtf:
        genes_for_gtf = list(dict.fromkeys(
            list(panel_genes or []) + [str(g) for _, d in hit_dfs for g in d['gene']]))
        print(f"Parsing canonical junctions for {len(genes_for_gtf)} gene(s) from {args.gtf} ...")
        canonical = {g: d['canonical_junctions']
                     for g, d in _parse_gtf_junctions(args.gtf, genes_for_gtf).items()}
    print("Collecting tested canonical junctions (GTEx + cohort)...")
    tested = _tested_canonical_junctions(canonical, args.samples, args.gtex_all_junction_files,
                                         args.cohort_results_dir)

    out_dfs = []
    for sample, hit_df in hit_dfs:
        n_hits = len(hit_df)
        hit_df = add_testability_and_non_hits(hit_df, sample, panel_genes, omim_df,
                                              coverage, canonical, tested)
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