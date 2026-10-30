# 9F - isoform usage outliers across a group. scored per isoform, and per gene for the overall isoform mix

import argparse
import concurrent.futures
import os
import time
import zlib

import numpy as np
import pandas as pd
from scipy.special import digamma, gammaln, polygamma
from scipy.stats import beta, betabinom
from statsmodels.stats.multitest import multipletests

from transcript_outliers import QC_FAIL, QC_PASS, call_outliers, outlier_gene_counts, qc_threshold

MODZ_CONST = 0.6745
MEANAD_CONST = 0.7979
OTHER = "other"
TRIM_MIN_SAMPLES = 30


def parse_args():
    p = argparse.ArgumentParser(
        description="Score transcript (isoform) usage outliers across a cohort from the AMALGAM transcript count "
                    "matrix. For each gene, samples with enough gene reads are compared with each other: every "
                    "well-used isoform is tested on its own, and the gene's whole isoform mix is tested as one. "
                    "beta_binomial (large cohorts): a beta distribution is fitted to each isoform's usage and the "
                    "sample's isoform/gene read counts are tested with a beta-binomial; a Dirichlet is fitted to the "
                    "isoform mix and the sample's counts are tested with a Dirichlet-multinomial (simulated p-value). "
                    "modified_zscore (small cohorts): modified z-scores of each isoform's usage and of the total "
                    "variation distance (TVD) of the sample's isoform mix from the cohort median. Every sample also "
                    "gets nearest_tvd: the TVD to the closest other sample (fewer than 30 samples), or to the "
                    "(k+1)-th closest after skipping ceil(1% of the others) (30 or more), the effect size that "
                    "isoform mix outliers are gated on. Sample QC: outliers are called with the 9G thresholds and "
                    "any sample whose number of genes with an isoform or isoform-mix outlier exceeds the group "
                    "median + --qc-mad-multiplier x MAD (MAD floored at 1) fails QC; the cohort is then refitted "
                    "with failing samples left out of every reference distribution, and they are still scored "
                    "against it.")
    p.add_argument("--tx-matrix", required=True,
                   help="transcript_amalgam_matrix_raw.tsv: transcript_id, gene_id (gene symbol), one count column "
                        "per sample.")
    p.add_argument("--outprefix", required=True,
                   help="Writes <outprefix>_transcript_usage_scores.tsv.gz, "
                        "<outprefix>_gene_isoform_landscape_scores.tsv.gz, <outprefix>_transcript_usage_matrix.tsv and "
                        "<outprefix>_transcript_usage_gene_status.tsv (every gene: tested or why not, and each "
                        "sample's gene read count).")
    p.add_argument("--method", required=True, choices=["beta_binomial", "modified_zscore"])
    p.add_argument("--n-threshold", type=int, required=True,
                   help="Minimum number of samples passing --min-gene-reads for a gene to be tested.")
    p.add_argument("--min-gene-reads", type=float, default=50,
                   help="A sample is used for a gene only if it has at least this many reads on the gene.")
    p.add_argument("--min-usage", type=float, default=0.05,
                   help="An isoform is tested on its own if, in at least one used sample, it has at least this "
                        "usage and --min-reads reads; all other isoforms of the gene are pooled into 'other'.")
    p.add_argument("--min-reads", type=float, default=10)
    p.add_argument("--rescale-factor", type=float, default=1e-3,
                   help="Usage is rescaled away from 0 and 1 before fitting: u*(1-K*eps)+eps for K categories.")
    p.add_argument("--dm-simulations", type=int, default=10000,
                   help="Maximum Dirichlet-multinomial simulations per sample and gene (beta_binomial method).")
    p.add_argument("--padj", type=float, default=0.05, help="Sample QC outlier calling (as 9G): maximum padj.")
    p.add_argument("--z", type=float, default=3.5, help="Sample QC outlier calling (as 9G): minimum |modified z|.")
    p.add_argument("--delta", type=float, default=0.1, help="Sample QC outlier calling (as 9G): minimum |delta_usage|.")
    p.add_argument("--tvd", type=float, default=0.1, help="Sample QC outlier calling (as 9G): minimum nearest_tvd.")
    p.add_argument("--qc-mad-multiplier", type=float, default=5.0,
                   help="A sample fails QC when its outlier gene count exceeds median + this x MAD. 0 turns QC off.")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threads", type=int, default=1)
    return p.parse_args()


# 1st-99th percentile of everyone except the sample itself, otherwise an outlier widens its own range
def leave_one_out_range(x, ref):
    lo, hi = np.full(len(x), np.nan), np.full(len(x), np.nan)
    vals = np.sort(x[ref], kind="stable")
    n = len(vals)
    for j in range(len(x)):
        if ref[j]:
            others = np.delete(vals, np.searchsorted(vals, x[j]))
        else:
            others = vals
        m = len(others)
        if m < 1:
            continue
        k = 0 if n < TRIM_MIN_SAMPLES else int(np.ceil((n - 1) * 0.01))
        k = min(k, (m - 1) // 2)
        lo[j], hi[j] = others[k], others[m - 1 - k]
    return lo, hi


# how far each sample's isoform mix is from its closest neighbour
def nearest_tvd(usage, ref):
    n_ref = int(ref.sum())
    dist = 0.5 * np.abs(usage[:, :, None] - usage[:, None, :]).sum(axis=0)
    np.fill_diagonal(dist, np.inf)
    dist = dist[:, ref]
    k = 0 if n_ref < TRIM_MIN_SAMPLES else int(np.ceil((n_ref - 1) * 0.01))
    return np.sort(dist, axis=1)[:, min(k, max(n_ref - 2, 0))]


def fit_beta(x, tol):
    var = float(np.var(x))
    if var < tol:
        m = float(np.clip(np.median(x), tol, 1.0 - tol))
        k = max(m * (1.0 - m) / tol - 1.0, tol)
        return m * k, (1.0 - m) * k
    try:
        a, b, *_ = beta.fit(x, floc=0, fscale=1)
        if np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0:
            return float(a), float(b)
    except Exception:
        pass
    mean = float(np.mean(x))
    k = max(mean * (1.0 - mean) / var - 1.0, tol)
    return mean * k, (1.0 - mean) * k


def betabinom_two_sided(k, n, a, b):
    lte = betabinom.cdf(k, n, a, b)
    gte = betabinom.cdf(n - k, n, b, a)
    return float(np.clip(2.0 * min(lte, gte), 0.0, 1.0))


def modz_fit(x, tol):
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    if mad < tol:
        meanad = float(np.mean(np.abs(x - med)))
        if meanad < tol:
            return med, None
        mad = meanad * (MODZ_CONST / MEANAD_CONST)
    return med, mad


# Dirichlet-multinomial fit for the whole isoform mix + a Monte Carlo p-value
def inv_digamma(y):
    x = np.where(y >= -2.22, np.exp(y) + 0.5, -1.0 / (y - digamma(1.0)))
    for _ in range(5):
        x = x - (digamma(x) - y) / polygamma(1, x)
    return x


def fit_dirichlet(P):
    logp = np.log(P).mean(axis=0)
    mean = P.mean(axis=0)
    var = P[:, 0].var()
    s = mean[0] * (1 - mean[0]) / var - 1 if var > 0 else 100.0
    alpha = mean * (s if np.isfinite(s) and s > 0 else 1.0)
    for _ in range(1000):
        new = inv_digamma(digamma(alpha.sum()) + logp)
        if np.max(np.abs(new - alpha)) < 1e-9 * np.max(alpha):
            alpha = new
            break
        alpha = new
    return alpha


def dm_logpmf(x, alpha):
    n = x.sum(axis=-1)
    a0 = alpha.sum()
    return (gammaln(n + 1) - gammaln(x + 1).sum(axis=-1) + gammaln(a0) - gammaln(n + a0)
            + (gammaln(x + alpha) - gammaln(alpha)).sum(axis=-1))


def dm_pvalue(x, alpha, max_sims, rng):
    n = int(x.sum())
    if n == 0:
        return 1.0
    obs = dm_logpmf(x.astype(float), alpha) + 1e-9
    done, extreme, batch = 0, 0, 1000
    while done < max_sims:
        b = min(batch, max_sims - done)
        probs = rng.dirichlet(alpha, size=b)
        probs = np.clip(probs, 0, None)
        probs /= probs.sum(axis=1, keepdims=True)
        sims = rng.multinomial(n, probs)
        extreme += int((dm_logpmf(sims.astype(float), alpha) <= obs).sum())
        done += b
        if (extreme + 1) / (done + 1) > 0.1:
            break
        batch = max_sims
    return (extreme + 1) / (done + 1)


def score_gene(task):
    gene, tx_ids, counts, samples, ref_all, a = task
    eps = a["rescale_factor"]
    gene_reads = counts.sum(axis=0)
    use = gene_reads >= a["min_gene_reads"]
    n_used = int((use & ref_all).sum())
    status = dict(gene=gene, n_samples_used=n_used, n_isoforms_tested=0,
                  **{s: r for s, r in zip(samples, np.round(gene_reads, 2))})
    if n_used < a["n_threshold"]:
        return gene, None, None, None, dict(status, status="low_n")
    # low usage / low count isoforms get pooled into 'other', they're mostly noise
    C = counts[:, use]
    N = gene_reads[use]
    S = [s for s, u in zip(samples, use) if u]
    ref = ref_all[use]
    n_scored = len(S)
    U = C / N
    keep = ((U[:, ref] >= a["min_usage"]) & (C[:, ref] >= a["min_reads"])).any(axis=1)
    if not keep.any():
        return gene, None, None, None, dict(status, status="no_isoform")
    names = [t for t, k in zip(tx_ids, keep) if k]
    cat_counts = C[keep]
    other = C[~keep].sum(axis=0)
    if (~keep).any() and other.sum() > 0:
        names = names + [OTHER]
        cat_counts = np.vstack([cat_counts, other])
    K = len(names)
    if K < 2:
        return gene, None, None, None, dict(status, status="single_isoform", n_isoforms_tested=1)
    cat_usage = cat_counts / N
    n_tested = int(keep.sum())
    rng = np.random.default_rng([a["seed"], zlib.crc32(gene.encode())])

    iso_rows = []
    for i in range(n_tested):
        u = cat_usage[i]
        r = u * (1 - 2 * eps) + eps
        p1, p99 = leave_one_out_range(u, ref)
        delta = np.where(u > p99, u - p99, np.where(u < p1, u - p1, 0.0))
        base = dict(gene=gene, transcript_id=names[i], n=n_used)
        if a["method"] == "beta_binomial":
            al, be = fit_beta(r[ref], eps)
            for j in range(n_scored):
                k = int(round(cat_counts[i, j]))
                nn = int(round(N[j]))
                iso_rows.append(dict(sample=S[j], **base, gene_reads=N[j], transcript_reads=cat_counts[i, j],
                                     usage=u[j], p1=p1[j], p99=p99[j], delta_usage=delta[j], alpha=al, beta=be,
                                     expected_usage=al / (al + be),
                                     p_value=betabinom_two_sided(min(k, nn), nn, al, be)))
        else:
            med, mad = modz_fit(r[ref], eps)
            for j in range(n_scored):
                z = MODZ_CONST * (r[j] - med) / mad if mad else np.nan
                iso_rows.append(dict(sample=S[j], **base, gene_reads=N[j], transcript_reads=cat_counts[i, j],
                                     usage=u[j], p1=p1[j], p99=p99[j], delta_usage=delta[j], median=med,
                                     mad=mad if mad else "no_variance", modz=z if mad else "no_variance"))

    # per gene: TVD between each sample's mix and the median of the others
    tvd = np.empty(n_scored)
    shift_tx, shift_val = [], []
    for j in range(n_scored):
        others = ref.copy()
        others[j] = False
        center = np.median(cat_usage[:, others], axis=1)
        center = center / center.sum() if center.sum() > 0 else np.full(K, 1.0 / K)
        d = cat_usage[:, j] - center
        tvd[j] = 0.5 * np.abs(d).sum()
        top = int(np.argmax(np.abs(d)))
        shift_tx.append(names[top])
        shift_val.append(d[top])

    near = nearest_tvd(cat_usage, ref)
    land_rows = []
    other_usage = cat_usage[-1] if names[-1] == OTHER else np.zeros(n_scored)
    base = dict(gene=gene, n_isoforms_tested=n_tested, n=n_used)
    if a["method"] == "beta_binomial":
        P = (cat_usage * (1 - K * eps) + eps).T
        alpha = fit_dirichlet(P[ref])
        X = np.rint(cat_counts).astype(np.int64).T
        for j in range(n_scored):
            land_rows.append(dict(sample=S[j], **base, gene_reads=N[j], other_usage=other_usage[j], tvd=tvd[j],
                                  nearest_tvd=near[j], max_shift_transcript=shift_tx[j], max_shift=shift_val[j],
                                  dirichlet_precision=float(alpha.sum()),
                                  p_value=dm_pvalue(X[j], alpha, a["dm_simulations"], rng)))
    else:
        med, mad = modz_fit(tvd[ref], eps)
        for j in range(n_scored):
            z = MODZ_CONST * (tvd[j] - med) / mad if mad else np.nan
            land_rows.append(dict(sample=S[j], **base, gene_reads=N[j], other_usage=other_usage[j], tvd=tvd[j],
                                  nearest_tvd=near[j], max_shift_transcript=shift_tx[j], max_shift=shift_val[j],
                                  tvd_median=med,
                                  tvd_mad=mad if mad else "no_variance", modz=z if mad else "no_variance"))

    mat = pd.DataFrame(np.full((K, len(samples)), np.nan), columns=samples)
    mat.loc[:, S] = cat_usage
    mat.insert(0, "transcript_id", names)
    mat.insert(0, "gene", gene)
    return gene, iso_rows, land_rows, mat, dict(status, status="tested", n_isoforms_tested=n_tested)


def add_padj(df):
    if df.empty or "p_value" not in df.columns:
        return df
    df["padj"] = np.nan
    for _, idx in df.groupby("sample").groups.items():
        p = df.loc[idx, "p_value"].to_numpy(dtype=float)
        ok = np.isfinite(p)
        if ok.any():
            adj = np.full(len(p), np.nan)
            adj[ok] = multipletests(p[ok], method="fdr_bh")[1]
            df.loc[idx, "padj"] = adj
    return df


def score_cohort(tx, samples, ref_all, settings, threads):
    tasks = [(str(g), sub["transcript_id"].astype(str).tolist(), sub[samples].to_numpy(dtype=float), samples,
              ref_all, settings) for g, sub in tx.groupby("gene_id", sort=True)]
    iso, land, mats, status, status_rows = [], [], [], {}, []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, threads)) as pool:
        for gene, iso_rows, land_rows, mat, info in pool.map(score_gene, tasks, chunksize=8):
            st = info["status"]
            status[st] = status.get(st, 0) + 1
            status_rows.append(info)
            if st == "tested":
                iso.extend(iso_rows)
                land.extend(land_rows)
                mats.append(mat)
    print("Genes: " + ", ".join(f"{v} {k}" for k, v in sorted(status.items())))
    if settings["method"] == "beta_binomial":
        iso_cols = ["alpha", "beta", "expected_usage", "p_value"]
        land_cols = ["dirichlet_precision", "p_value"]
    else:
        iso_cols = ["median", "mad", "modz"]
        land_cols = ["tvd_median", "tvd_mad", "modz"]
    iso_df = pd.DataFrame(iso, columns=["sample", "gene", "transcript_id", "gene_reads", "transcript_reads", "usage",
                                        "n", "p1", "p99", "delta_usage"] + iso_cols)
    land_df = pd.DataFrame(land, columns=["sample", "gene", "gene_reads", "n_isoforms_tested", "other_usage", "n",
                                          "tvd", "nearest_tvd", "max_shift_transcript", "max_shift"] + land_cols)
    return add_padj(iso_df), add_padj(land_df), mats, status_rows


def main():
    print("\n\n\n******************************************************************************************")
    print("Scoring transcript usage outliers...")
    print("******************************************************************************************\n")
    args = parse_args()
    t0 = time.time()
    tx = pd.read_csv(args.tx_matrix, sep="\t")
    samples = [c for c in tx.columns if c not in ("transcript_id", "gene_id")]
    print(f"{len(tx)} transcripts, {tx['gene_id'].nunique()} genes, {len(samples)} samples; method {args.method}")
    settings = dict(method=args.method, n_threshold=args.n_threshold, min_gene_reads=args.min_gene_reads,
                    min_usage=args.min_usage, min_reads=args.min_reads, rescale_factor=args.rescale_factor,
                    dm_simulations=args.dm_simulations, seed=args.seed)
    ref_all = np.ones(len(samples), dtype=bool)
    iso_df, land_df, mats, status_rows = score_cohort(tx, samples, ref_all, settings, args.threads)

    # sample QC: samples with way too many outlier genes get dropped and everything is refit
    iso_out, land_out = call_outliers(iso_df, land_df, args.method, args.padj, args.z, args.delta, args.tvd)
    counts = outlier_gene_counts(iso_out, land_out, samples)
    flagged = []
    thr, med, mad = (np.inf, np.nan, np.nan)
    if args.qc_mad_multiplier > 0:
        thr, med, mad = qc_threshold(counts.values, args.qc_mad_multiplier)
        flagged = [s for s, n in counts.items() if n > thr]
    print(f"Sample QC: genes with an isoform outlier per sample: median {med:g}, MAD {mad:g}, "
          f"threshold {thr:g}; {len(flagged)} sample(s) fail: {', '.join(flagged) if flagged else 'none'}")
    if flagged:
        ref_all = np.array([s not in set(flagged) for s in samples])
        print(f"Refitting with {int(ref_all.sum())} reference sample(s) (failing samples are scored but not used "
              f"as reference)...")
        iso_df, land_df, mats, status_rows = score_cohort(tx, samples, ref_all, settings, args.threads)
    qc = pd.DataFrame(dict(sample=samples, n_outlier_genes=[int(counts.get(str(s), 0)) for s in samples],
                           median=med, mad=mad, threshold=thr,
                           qc=[QC_FAIL if s in set(flagged) else QC_PASS for s in samples],
                           reference=ref_all))

    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)), exist_ok=True)
    iso_df.to_csv(args.outprefix + "_transcript_usage_scores.tsv.gz", sep="\t", index=False)
    land_df.to_csv(args.outprefix + "_gene_isoform_landscape_scores.tsv.gz", sep="\t", index=False)
    mat_df = pd.concat(mats, ignore_index=True) if mats else pd.DataFrame(columns=["gene", "transcript_id"] + samples)
    mat_df.to_csv(args.outprefix + "_transcript_usage_matrix.tsv", sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(status_rows, columns=["gene", "status", "n_samples_used", "n_isoforms_tested"] + samples).to_csv(
        args.outprefix + "_transcript_usage_gene_status.tsv", sep="\t", index=False)
    qc.to_csv(args.outprefix + "_transcript_usage_sample_qc.tsv", sep="\t", index=False)
    print(f"{len(iso_df)} isoform rows, {len(land_df)} gene rows written to {args.outprefix}_* "
          f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
