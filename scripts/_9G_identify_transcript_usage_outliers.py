# 9G - filter the 9F scores down to isoform / landscape outliers, with sample QC and recurrence

import argparse
import os

import pandas as pd

from sample_alias import add_alias_map_arg, parse_alias_map, resolve_all
from transcript_outliers import QC_PASS, call_outliers, load_sample_qc


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter the transcript usage scores to outliers. An isoform is an outlier when it is "
                    "significant (padj, or |modified z|) and its usage lies at least --delta outside the range of "
                    "the other samples. A gene's isoform mix is an outlier when it is significant (padj, or "
                    "modified z of its TVD from the cohort median) and its TVD to the nearest other sample "
                    "(nearest_tvd) is at least --tvd. n_outlier_samples counts how many samples share the "
                    "outlier (the same isoform, or the same gene for isoform-mix outliers); samples that failed 9F's "
                    "sample QC are left out of that count. qc is pass/fail from 9F's sample QC. Alias-labelled copies "
                    "of both outputs are written too.")
    p.add_argument("--isoform-scores", required=True, help="<prefix>_transcript_usage_scores.tsv.gz from 9F.")
    p.add_argument("--landscape-scores", required=True, help="<prefix>_gene_isoform_landscape_scores.tsv.gz from 9F.")
    p.add_argument("--method", required=True, choices=["beta_binomial", "modified_zscore"])
    p.add_argument("--padj", type=float, default=0.05, help="beta_binomial: maximum padj.")
    p.add_argument("--z", type=float, default=3.5, help="modified_zscore: minimum |modified z|.")
    p.add_argument("--delta", type=float, default=0.1, help="Minimum |delta_usage| for an isoform outlier.")
    p.add_argument("--tvd", type=float, default=0.1,
                   help="Minimum TVD to the nearest other sample for a gene isoform-mix outlier (default 0.1).")
    p.add_argument("--sample-qc", default=None, help="<prefix>_transcript_usage_sample_qc.tsv from 9F.")
    p.add_argument("--isoform-outfile", required=True)
    p.add_argument("--landscape-outfile", required=True)
    p.add_argument("--isoform-alias-outfile", required=True)
    p.add_argument("--landscape-alias-outfile", required=True)
    add_alias_map_arg(p)
    return p.parse_args()


def main():
    args = parse_args()
    iso = pd.read_csv(args.isoform_scores, sep="\t")
    land = pd.read_csv(args.landscape_scores, sep="\t")

    iso_out, land_out = call_outliers(iso, land, args.method, args.padj, args.z, args.delta, args.tvd)
    # recurrence only counts QC-passing samples, otherwise one bad sample makes everything look recurrent
    qc = load_sample_qc(args.sample_qc) if args.sample_qc and os.path.isfile(args.sample_qc) else {}
    iso_out = iso_out.assign(qc=iso_out["sample"].astype(str).map(qc).fillna(QC_PASS))
    land_out = land_out.assign(qc=land_out["sample"].astype(str).map(qc).fillna(QC_PASS))
    ok_iso = iso_out[iso_out["qc"] == QC_PASS]
    ok_land = land_out[land_out["qc"] == QC_PASS]
    iso_out = iso_out.assign(n_outlier_samples=iso_out.set_index(["gene", "transcript_id"]).index.map(
        ok_iso.groupby(["gene", "transcript_id"])["sample"].nunique()).fillna(0).astype(int))
    land_out = land_out.assign(n_outlier_samples=land_out["gene"].map(
        ok_land.groupby("gene")["sample"].nunique()).fillna(0).astype(int))
    failed = sorted({s for s, q in qc.items() if q != QC_PASS})
    if failed:
        print(f"{len(failed)} sample(s) failed isoform sample QC; their calls are kept with qc=fail: "
              + ", ".join(failed))
    alias_map = parse_alias_map(args.alias_map)
    for df, path, alias_path in ((iso_out, args.isoform_outfile, args.isoform_alias_outfile),
                                 (land_out, args.landscape_outfile, args.landscape_alias_outfile)):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        df.to_csv(path, sep="\t", index=False)
        df.assign(sample=resolve_all(df["sample"].astype(str), alias_map)).to_csv(alias_path, sep="\t", index=False)
    print(f"{len(iso_out)} isoform outlier(s) in {iso_out['sample'].nunique()} sample(s) -> {args.isoform_outfile}")
    print(f"{len(land_out)} gene isoform-mix outlier(s) in {land_out['sample'].nunique()} sample(s) "
          f"-> {args.landscape_outfile}")


if __name__ == "__main__":
    main()
