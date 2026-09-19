"""
rules/5_junction_analysis.smk
Per-sample splice junction outlier analysis vs. GTEx reference tissues.
See docs/rules/5_junction_analysis.md for details.
"""

# ---------------------------------------------------------------------------
# Get splice junction counts  [always 1 thread]
# ---------------------------------------------------------------------------
rule _5A_get_junction_counts:
    input:
        mapping = "{outdir}/phased_reads/{sample}_gene_bam_mapping_file.tsv",
        bam     = lambda wc: SAMPLES[wc.sample]["bam"],
    output:
        jxn_counts = "{outdir}/junction_analysis/junction_counts/{sample}_splice_junction_counts.tsv",
        jxn_matrix = "{outdir}/junction_analysis/junction_counts/{sample}_junction_count_matrix.tsv",
    params:
        script_jxn   = workflow.basedir + "/scripts/get_splice_junction_counts_by_region.py",
        script_matrix= workflow.basedir + "/scripts/make_junction_count_matrix.py",
    threads: 1
    resources:
        # make_junction_count_matrix.py scans the whole BAM genome-wide into
        # memory (unlike get_splice_junction_counts_by_region.py, which is
        # restricted to the BED panel), so this can OOM on large BAMs. Scale
        # with BAM size, same pattern as the variant callers.
        mem_mb     = lambda wc, attempt: max(4096, attempt * _bam_size_gb(wc) * 2 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_junction_counts.log"
    shell:
        """
        mkdir -p $(dirname {output.jxn_counts})
        python -u {params.script_jxn} \\
            --mapping-file {input.mapping} \\
            --outfile      {output.jxn_counts} \\
        2>&1 | tee {log}
        python -u {params.script_matrix} \\
            --bam     {input.bam} \\
            --outfile {output.jxn_matrix} \\
        2>&1 | tee -a {log}
        """


# ---------------------------------------------------------------------------
# Fit beta distributions on GTEx reference junctions -- genome-wide, once
# per (tissue, thresholds) for the whole run, NOT per sample. See
# scripts/fit_gtex_beta_distributions.py's module docstring for why this is
# safe to share: the fit depends only on GTEx data + these three thresholds,
# never on any individual sample's own read counts. Every _5C job for a
# given tissue depends on this same output, so Snakemake computes it once
# and every sample's job waits on/reuses it rather than each one refitting
# from scratch (which is what this pipeline used to do).
# ---------------------------------------------------------------------------
def _gtex_file(tissue):
    return (str(config['gtex_data_dir']) + '/gtex_' + str(tissue) + '_jxn_counts.txt')

# A plain string TEMPLATE containing the literal wildcard token "{tissue}"
# (config values are baked in now, since they're the same for the whole
# run; "tissue" stays a real Snakemake wildcard, filled in by whichever
# rule instance references it). Snakemake's `output:` must be a static
# pattern it can match filenames against -- not a Python callable -- so
# this is used directly as a string, the same way every other
# "{outdir}/..." path in this pipeline is, rather than as
# `lambda wc: _gtex_beta_fits_file(wc.tissue)`.
#
# Lives under this RUN's own output_dir (not the shared gtex_data_dir
# resource folder) -- still shared across every sample in this run (one
# file per tissue+thresholds, not per sample), just scoped to this run's
# output rather than the reference-data directory.
_GTEX_BETA_FITS_TEMPLATE = (
    str(config['output_dir']) + '/gtex_beta_distributions/gtex_{tissue}_beta_fits'
    '_cov' + str(config['gtex_coverage_threshold']) +
    '_n' + str(config['gtex_n_threshold']) +
    '_tol' + str(config['PSI_rescale_factor']) + '.tsv'
)

def _gtex_fit_mem(wc, threads, attempt):
    return max(4096, attempt * threads * (6 if wc.tissue == "brain" else 2) * 1024)


rule _5B1_fit_gtex_beta_distributions:
    input:
        gtex_file = lambda wc: _gtex_file(wc.tissue),
        # Union of every sample's BED panel in this run -- restricts fitting
        # to junctions any sample here could actually query (see
        # all_bed_files()'s docstring in the Snakefile). Only narrows WHICH
        # junctions get the expensive beta.fit() call; coverage/PSI are
        # still computed genome-wide first so the kept junctions' numbers
        # are unchanged from a full genome-wide fit.
        bed_files = all_bed_files(),
    output:
        beta_fits = _GTEX_BETA_FITS_TEMPLATE,
    params:
        gtex_cov  = config["gtex_coverage_threshold"],
        gtex_n    = config["gtex_n_threshold"],
        psi_tol   = config["PSI_rescale_factor"],
        script    = workflow.basedir + "/scripts/fit_gtex_beta_distributions.py",
    # No {sample} wildcard on this rule (it runs once per tissue for the
    # whole run, not per sample) -- _rule_threads() indexes SAMPLES[wc.sample]
    # and would raise on a wildcards object with no .sample attribute here,
    # so this uses a plain config lookup instead, same fallback-to-default
    # pattern as _rule_threads(), just without the per-sample override.
    threads: lambda wc: int(config.get("fit_gtex_beta_distributions_threads", config["threads"]))
    resources:
        mem_mb    = _gtex_fit_mem,
        runtime   = config["time"],
    log:
        str(config['output_dir']) + "/gtex_beta_distributions/logs/gtex_{tissue}_fit_beta_distributions.log"
    shell:
        """
        mkdir -p $(dirname {output.beta_fits})
        python -u {params.script} \\
            --gtexfile              {input.gtex_file} \\
            --bed-files               {input.bed_files} \\
            --outfile                {output.beta_fits} \\
            --gtex-coverage-threshold {params.gtex_cov} \\
            --gtex-n-threshold        {params.gtex_n} \\
            --PSI-rescale-factor      {params.psi_tol} \\
            --threads                 {threads} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Fit beta distributions for this sample's own junctions that are entirely
# absent from the GTEx reference matrix ("novel" junctions) -- can't be
# precomputed ahead of time the way _5B1 is, since which junctions are
# novel depends on this specific sample's own splicing. See
# scripts/fit_novel_junction_beta_distributions.py's module docstring.
# [per-sample configurable via fit_novel_junction_beta_distributions_threads]
# ---------------------------------------------------------------------------
rule _5B2_fit_novel_junction_beta_distributions:
    input:
        jxn_counts = "{outdir}/junction_analysis/junction_counts/{sample}_splice_junction_counts.tsv",
        gtex_file  = lambda wc: _gtex_file(wc.tissue),
        gtex_fits  = _GTEX_BETA_FITS_TEMPLATE,
    output:
        novel_fits = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_novel_junction_beta_fits.tsv",
    params:
        gtex_cov  = config["gtex_coverage_threshold"],
        gtex_n    = config["gtex_n_threshold"],
        psi_tol   = config["PSI_rescale_factor"],
        script    = workflow.basedir + "/scripts/fit_novel_junction_beta_distributions.py",
    threads: lambda wc: _rule_threads(wc, "fit_novel_junction_beta_distributions")
    resources:
        mem_mb    = _gtex_fit_mem,
        runtime   = config["time"],
    log:
        "{outdir}/../logs/{sample}_{tissue}_fit_novel_junction_beta_distributions.log"
    shell:
        """
        mkdir -p $(dirname {output.novel_fits})
        python -u {params.script} \\
            --jxn-info-file           {input.jxn_counts} \\
            --gtexfile                {input.gtex_file} \\
            --gtex-beta-fits          {input.gtex_fits} \\
            --outfile                 {output.novel_fits} \\
            --gtex-coverage-threshold {params.gtex_cov} \\
            --gtex-n-threshold        {params.gtex_n} \\
            --PSI-rescale-factor      {params.psi_tol} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Beta-binomial tests vs GTEx (one job per sample × tissue). No longer fits
# anything itself -- merges the precomputed _5B1/_5B2 fits against this
# sample's junction data and runs the significance test.
# [per-sample configurable via perform_binomial_tests_threads in run config]
# ---------------------------------------------------------------------------
rule _5C_perform_binomial_tests:
    input:
        jxn_counts = "{outdir}/junction_analysis/junction_counts/{sample}_splice_junction_counts.tsv",
        gtex_fits  = _GTEX_BETA_FITS_TEMPLATE,
        novel_fits = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_novel_junction_beta_fits.tsv",
    output:
        all_jxns = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_all_junctions.tsv",
    params:
        annotation = config["annotation"],
        samp_cov   = config["sample_coverage_threshold"],
        gtex_n     = config["gtex_n_threshold"],
        psi_tol    = config["PSI_rescale_factor"],
        phasing_thr= config["jxn_phasing_threshold"],
        script     = workflow.basedir + "/scripts/perform_splice_junction_beta_binomial_tests.py",
    threads: lambda wc: _rule_threads(wc, "perform_binomial_tests")
    resources:
        # No longer reads the raw (large) per-GTEx-sample count matrix --
        # only the small precomputed fit tables -- so this no longer needs
        # to scale with tissue size the way _5B1/_5B2 do.
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * threads * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_{tissue}_betabinom.log"
    shell:
        """
        mkdir -p $(dirname {output.all_jxns})
        python -u {params.script} \\
            --jxn-info-file           {input.jxn_counts} \\
            --gtex-beta-fits          {input.gtex_fits} \\
            --novel-beta-fits         {input.novel_fits} \\
            --outfile                 {output.all_jxns} \\
            --sample-coverage-threshold {params.samp_cov} \\
            --gtex-n-threshold          {params.gtex_n} \\
            --PSI-rescale-factor        {params.psi_tol} \\
            --phasing-threshold         {params.phasing_thr} \\
            --annotation-file           {params.annotation} \\
            --threads                   {threads} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Filter to significant outlier junctions
# [per-sample configurable via identify_junction_outliers_threads in run config]
# ---------------------------------------------------------------------------
rule _5D_identify_junction_outliers:
    input:
        all_jxns = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_all_junctions.tsv",
    output:
        outliers = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_outlier_junctions.tsv",
    params:
        padj_thr   = config["padj_threshold"],
        dpsi_thr   = config["delta_psi_threshold"],
        script     = workflow.basedir + "/scripts/identify_splice_junction_outliers.py",
    threads: 1   # single vectorized pandas filter, no per-item work to parallelize -- see the script's own header comment
    resources:
        mem_mb     = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_{tissue}_jxn_outliers.log"
    shell:
        """
        python -u {params.script} \\
            --infile             {input.all_jxns} \\
            --outfile            {output.outliers} \\
            --padj-threshold     {params.padj_thr} \\
            --delta-PSI-threshold {params.dpsi_thr} \\
        2>&1 | tee {log}
        """
