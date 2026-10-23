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
                    "isoform mix outliers are gated on.")
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
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threads", type=int, default=1)
    return p.parse_args()


def leave_one_out_range(x):
    n = len(x)
    order = np.argsort(x, kind="stable")
    vals = x[order]
    rank = np.empty(n, dtype=np.int64)
    rank[order] = np.arange(n)
    m = n - 1
    if m < 1:
        return np.full(n, np.nan), np.full(n, np.nan)
    k = 0 if n < TRIM_MIN_SAMPLES else int(np.ceil(m * 0.01))
    lo = k + (k >= rank)
    hi = (m - 1 - k) + ((m - 1 - k) >= rank)
    return vals[lo], vals[hi]


def nearest_tvd(usage):
    n = usage.shape[1]
    dist = 0.5 * np.abs(usage[:, :, None] - usage[:, None, :]).sum(axis=0)
    np.fill_diagonal(dist, np.inf)
    k = 0 if n < TRIM_MIN_SAMPLES else int(np.ceil((n - 1) * 0.01))
    return np.sort(dist, axis=1)[:, min(k, n - 2)]


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
    gene, tx_ids, counts, samples, a = task
    eps = a["rescale_factor"]
    gene_reads = counts.sum(axis=0)
    use = gene_reads >= a["min_gene_reads"]
    n_used = int(use.sum())
    status = dict(gene=gene, n_samples_used=n_used, n_isoforms_tested=0,
                  **{s: r for s, r in zip(samples, np.round(gene_reads, 2))})
    if n_used < a["n_threshold"]:
        return gene, None, None, None, dict(status, status="low_n")
    C = counts[:, use]
    N = gene_reads[use]
    S = [s for s, u in zip(samples, use) if u]
    U = C / N
    keep = ((U >= a["min_usage"]) & (C >= a["min_reads"])).any(axis=1)
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
        p1, p99 = leave_one_out_range(u)
        delta = np.where(u > p99, u - p99, np.where(u < p1, u - p1, 0.0))
        base = dict(gene=gene, transcript_id=names[i], n=n_used)
        if a["method"] == "beta_binomial":
            al, be = fit_beta(r, eps)
            for j in range(n_used):
                k = int(round(cat_counts[i, j]))
                nn = int(round(N[j]))
                iso_rows.append(dict(sample=S[j], **base, gene_reads=N[j], transcript_reads=cat_counts[i, j],
                                     usage=u[j], p1=p1[j], p99=p99[j], delta_usage=delta[j], alpha=al, beta=be,
                                     expected_usage=al / (al + be),
                                     p_value=betabinom_two_sided(min(k, nn), nn, al, be)))
        else:
            med, mad = modz_fit(r, eps)
            for j in range(n_used):
                z = MODZ_CONST * (r[j] - med) / mad if mad else np.nan
                iso_rows.append(dict(sample=S[j], **base, gene_reads=N[j], transcript_reads=cat_counts[i, j],
                                     usage=u[j], p1=p1[j], p99=p99[j], delta_usage=delta[j], median=med,
                                     mad=mad if mad else "no_variance", modz=z if mad else "no_variance"))

    tvd = np.empty(n_used)
    shift_tx, shift_val = [], []
    for j in range(n_used):
        ref = np.median(np.delete(cat_usage, j, axis=1), axis=1)
        ref = ref / ref.sum() if ref.sum() > 0 else np.full(K, 1.0 / K)
        d = cat_usage[:, j] - ref
        tvd[j] = 0.5 * np.abs(d).sum()
        top = int(np.argmax(np.abs(d)))
        shift_tx.append(names[top])
        shift_val.append(d[top])

    near = nearest_tvd(cat_usage)
    land_rows = []
    other_usage = cat_usage[-1] if names[-1] == OTHER else np.zeros(n_used)
    base = dict(gene=gene, n_isoforms_tested=n_tested, n=n_used)
    if a["method"] == "beta_binomial":
        P = (cat_usage * (1 - K * eps) + eps).T
        alpha = fit_dirichlet(P)
        X = np.rint(cat_counts).astype(np.int64).T
        for j in range(n_used):
            land_rows.append(dict(sample=S[j], **base, gene_reads=N[j], other_usage=other_usage[j], tvd=tvd[j],
                                  nearest_tvd=near[j], max_shift_transcript=shift_tx[j], max_shift=shift_val[j],
                                  dirichlet_precision=float(alpha.sum()),
                                  p_value=dm_pvalue(X[j], alpha, a["dm_simulations"], rng)))
    else:
        med, mad = modz_fit(tvd, eps)
        for j in range(n_used):
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
    tasks = [(str(g), sub["transcript_id"].astype(str).tolist(), sub[samples].to_numpy(dtype=float), samples,
              settings) for g, sub in tx.groupby("gene_id", sort=True)]

    iso, land, mats, status, status_rows = [], [], [], {}, []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
        for gene, iso_rows, land_rows, mat, info in pool.map(score_gene, tasks, chunksize=8):
            st = info["status"]
            status[st] = status.get(st, 0) + 1
            status_rows.append(info)
            if st == "tested":
                iso.extend(iso_rows)
                land.extend(land_rows)
                mats.append(mat)
    print("Genes: " + ", ".join(f"{v} {k}" for k, v in sorted(status.items())))

    if args.method == "beta_binomial":
        iso_cols = ["alpha", "beta", "expected_usage", "p_value"]
        land_cols = ["dirichlet_precision", "p_value"]
    else:
        iso_cols = ["median", "mad", "modz"]
        land_cols = ["tvd_median", "tvd_mad", "modz"]
    iso_df = pd.DataFrame(iso, columns=["sample", "gene", "transcript_id", "gene_reads", "transcript_reads", "usage",
                                        "n", "p1", "p99", "delta_usage"] + iso_cols)
    land_df = pd.DataFrame(land, columns=["sample", "gene", "gene_reads", "n_isoforms_tested", "other_usage", "n",
                                          "tvd", "nearest_tvd", "max_shift_transcript", "max_shift"] + land_cols)
    iso_df = add_padj(iso_df)
    land_df = add_padj(land_df)

    os.makedirs(os.path.dirname(os.path.abspath(args.outprefix)), exist_ok=True)
    iso_df.to_csv(args.outprefix + "_transcript_usage_scores.tsv.gz", sep="\t", index=False)
    land_df.to_csv(args.outprefix + "_gene_isoform_landscape_scores.tsv.gz", sep="\t", index=False)
    mat_df = pd.concat(mats, ignore_index=True) if mats else pd.DataFrame(columns=["gene", "transcript_id"] + samples)
    mat_df.to_csv(args.outprefix + "_transcript_usage_matrix.tsv", sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(status_rows, columns=["gene", "status", "n_samples_used", "n_isoforms_tested"] + samples).to_csv(
        args.outprefix + "_transcript_usage_gene_status.tsv", sep="\t", index=False)
    print(f"{len(iso_df)} isoform rows, {len(land_df)} gene rows written to {args.outprefix}_* "
          f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
