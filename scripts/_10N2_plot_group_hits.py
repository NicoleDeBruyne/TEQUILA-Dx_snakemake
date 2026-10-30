# 10N2 - stacked bars of hit genes per sample, coloured by evidence type

import argparse
import os, pandas as pd, matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator

from palette import RED_FAMILY, TEAL_FAMILY, PURPLE_FAMILY, GOLD_FAMILY, use_fonts
from common import read_hit_rows
use_fonts()

rcParams['pdf.fonttype'] = 42


def parse_args():

    parser = argparse.ArgumentParser(description="Plot clustermaps for feature counts")
    parser.add_argument("--infile", required=True, help="Input file containing candidate genes. Expected columns: \
                                    sample, gene, variant, pathogenic_variant, ASE, outlier_junction, and optionally max_AMELIE_score and AMELIE_rank")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--amelie_score_threshold", type=float, help="AMELIE score threshold for plots.")
    parser.add_argument("--amelie_rank_threshold", type=float, help="AMELIE rank threshold for plots.")

    return parser.parse_args()

def plot(df, samples, catcol, categories, colors, outfile, title, legend_title):
    counts = df.groupby(["sample", catcol])["gene"].nunique().unstack(fill_value=0)
    counts = counts.reindex(samples, columns=categories, fill_value=0)

    fig, ax = plt.subplots(figsize=(12,6))
    counts.plot(kind="bar", stacked=True, ax=ax, color=colors)
    ax.set_xticklabels(samples, rotation=90)
    ax.set_xlabel("")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel("Gene Count")
    ax.set_title(title)
    ax.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    barplot_file = f"{os.path.splitext(outfile)[0]}_barplot{os.path.splitext(outfile)[1]}"
    plt.savefig(barplot_file)
    plt.close(fig)
    print(f"Bar plot saved to {barplot_file}")


def main():

    args = parse_args()

    df = read_hit_rows(args.infile, sep='\t', keep_default_na=False)
    df = df.astype(object)
    df.fillna('.', inplace=True)

    isoform = (df['isoform_outlier'].astype(str).str.upper() == 'TRUE') if 'isoform_outlier' in df.columns else False
    # lump the RNA evidence together so each gene lands in one category
    df['RNA_dysregulation'] = (df['ASE'] | (df['outlier_junction'] != 'None') | isoform)
    df['final_grouping'] = df.apply(lambda row: 
        "RNA dysregulation with candidate variant" if row['RNA_dysregulation'] and row['variant'] else
        ("RNA dysregulation" if row['RNA_dysregulation'] else
        ("Candidate variant" if row['variant'] else "No findings")), axis=1
    )

    samples = sorted(df['sample'].unique())

    plot(df[df['pathogenic_variant']].copy(), samples, "pathogenic_variant", [True, False], [RED_FAMILY[0], RED_FAMILY[2]], f'{args.outdir}/genes_with_pathogenic_variant.pdf', "Candidate Genes with Pathogenic Variant", "Pathogenic Variant")
    plot(df, samples, "ASE", [True], [PURPLE_FAMILY[0]], f'{args.outdir}/genes_with_ASE.pdf', "Candidate Genes", "Allelic Imbalance")
    plot(df, samples, "outlier_junction", ["Strong", "Moderate", "Weak"], TEAL_FAMILY[:3], f'{args.outdir}/genes_with_outlier_junction.pdf', "Candidate Genes", "Alternative Splicing")
    plot(df, samples, "final_grouping", ["RNA dysregulation", "RNA dysregulation with candidate variant"], [GOLD_FAMILY[0], GOLD_FAMILY[2]], f'{args.outdir}/genes_with_RNA_dysregulation.pdf', "Candidate Genes", "Group")

    # and again, but only for genes that pass the AMELIE cutoff
    if args.amelie_score_threshold:
        filtered_df = df[pd.to_numeric(df['max_AMELIE_score'], errors='coerce') >= args.amelie_score_threshold]
        plot(filtered_df[filtered_df['pathogenic_variant']].copy(), samples, "pathogenic_variant", [True, False], [RED_FAMILY[0], RED_FAMILY[2]], f'{args.outdir}/genes_with_pathogenic_variant_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes with Pathogenic Variant, AMELIE Score >= {args.amelie_score_threshold}", "Pathogenic Variant")
        plot(filtered_df, samples, "ASE", [True], [PURPLE_FAMILY[0]], f'{args.outdir}/genes_with_ASE_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes, AMELIE Score >= {args.amelie_score_threshold}", "Allelic Imbalance")
        plot(filtered_df, samples, "outlier_junction", ["Strong", "Moderate", "Weak"], TEAL_FAMILY[:3], f'{args.outdir}/genes_with_outlier_junction_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes, AMELIE Score >= {args.amelie_score_threshold}", "Outlier Junction")
        plot(filtered_df, samples, "final_grouping", ["RNA dysregulation", "RNA dysregulation with candidate variant"], [GOLD_FAMILY[0], GOLD_FAMILY[2]], f'{args.outdir}/genes_with_RNA_dysregulation_amelie_score_{args.amelie_score_threshold}.pdf', f"Candidate Genes, AMELIE Score >= {args.amelie_score_threshold}", "Group")

    if args.amelie_rank_threshold:
        filtered_df = df[pd.to_numeric(df['AMELIE_rank'], errors='coerce') <= args.amelie_rank_threshold]
        plot(filtered_df[filtered_df['pathogenic_variant']].copy(), samples, "pathogenic_variant", [True, False], [RED_FAMILY[0], RED_FAMILY[2]], f'{args.outdir}/genes_with_pathogenic_variant_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes with Pathogenic Variant, AMELIE Rank <= {args.amelie_rank_threshold}", "Pathogenic Variant")
        plot(filtered_df, samples, "ASE", [True], [PURPLE_FAMILY[0]], f'{args.outdir}/genes_with_ASE_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes, AMELIE Rank <= {args.amelie_rank_threshold}", "Allelic Imbalance")
        plot(filtered_df, samples, "outlier_junction", ["Strong", "Moderate", "Weak"], TEAL_FAMILY[:3], f'{args.outdir}/genes_with_outlier_junction_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes, AMELIE Rank <= {args.amelie_rank_threshold}", "Outlier Junction")
        plot(filtered_df, samples, "final_grouping", ["RNA dysregulation", "RNA dysregulation with candidate variant"], [GOLD_FAMILY[0], GOLD_FAMILY[2]], f'{args.outdir}/genes_with_RNA_dysregulation_amelie_rank_{args.amelie_rank_threshold}.pdf', f"Candidate Genes, AMELIE Rank <= {args.amelie_rank_threshold}", "Group")

if __name__ == "__main__":
    main()
