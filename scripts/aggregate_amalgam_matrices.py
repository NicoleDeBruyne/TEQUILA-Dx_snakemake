
import argparse

import pandas as pd

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-sample AMALGAM transcript quantification files into "
                    "cohort-wide transcript_matrix.tsv and gene_matrix.tsv.")
    parser.add_argument('--infiles', required=True, nargs='+',
        help="Every sample's <sample>_transcript_quantification.tsv in this group.")
    parser.add_argument('--samples', required=True, nargs='+',
        help="Sample name for each --infiles entry, same order/length.")
    parser.add_argument('--outprefix', required=True,
        help="Writes <outprefix>_transcript_matrix.tsv and <outprefix>_gene_matrix.tsv "
             "(raw, genome-wide, gene_id-keyed) -- see module docstring.")
    add_alias_map_arg(parser)
    args = parser.parse_args()
    if len(args.infiles) != len(args.samples):
        parser.error("--infiles and --samples must have the same number of entries")
    return args


def main():
    args = parse_args()

    tx_index = None
    for f in args.infiles:
        ids = pd.read_csv(f, sep='\t', usecols=['transcript_id', 'gene_id'])
        idx = pd.MultiIndex.from_frame(ids)
        tx_index = idx.unique() if tx_index is None else tx_index.union(idx, sort=False)

    transcript_matrix = pd.DataFrame(index=tx_index)
    gene_tables = []
    for sample, f in zip(args.samples, args.infiles):
        df = pd.read_csv(f, sep='\t', usecols=['transcript_id', 'gene_id', 'count'])
        transcript_matrix[sample] = (
            df.set_index(['transcript_id', 'gene_id'])['count'].reindex(tx_index)
        )
        gene = df.groupby('gene_id')['count'].sum().to_frame(name=sample)
        gene_tables.append(gene)
        print(f"Processed {sample}", flush=True)

    transcript_matrix = transcript_matrix.fillna(0).reset_index()

    gene_matrix = pd.concat(gene_tables, axis=1).fillna(0)

    out_transcript = args.outprefix + '_transcript_matrix.tsv'
    out_gene = args.outprefix + '_gene_matrix.tsv'
    transcript_matrix.to_csv(out_transcript, sep='\t', index=False)
    gene_matrix.to_csv(out_gene, sep='\t')
    print(f"Done: {out_transcript} and {out_gene} written.", flush=True)

    alias_map = parse_alias_map(args.alias_map)
    alias_gene_matrix = gene_matrix.copy()
    alias_gene_matrix.columns = resolve_all(alias_gene_matrix.columns, alias_map)
    out_gene_alias = args.outprefix + '_gene_matrix_alias.tsv'
    alias_gene_matrix.to_csv(out_gene_alias, sep='\t')
    print(f"Done: {out_gene_alias} written.", flush=True)


if __name__ == '__main__':
    main()