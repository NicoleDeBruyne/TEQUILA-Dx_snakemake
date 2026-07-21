"""
rules/junction_analysis.smk
1. get_junction_counts          – per-sample junction counts + count matrix       [1 thread, always]
2. perform_binomial_tests       – beta-binomial tests against each GTEx tissue    [per-sample configurable]
3. identify_junction_outliers   – filter to significant outlier junctions         [per-sample configurable]

(the per-gene BAM mapping file these depend on is built once, in rules/phase_reads.smk)
"""

# ---------------------------------------------------------------------------
# Get splice junction counts  [always 1 thread]
# ---------------------------------------------------------------------------
rule get_junction_counts:
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
        # make_junction_count_matrix.py scans the ENTIRE BAM genome-wide into
        # in-memory dicts/dataframes (unlike get_splice_junction_counts_by_region.py,
        # which is restricted to the BED panel's regions), so this can get
        # silently OOM-killed on larger BAMs right at the tail end of the
        # job (after useful output has already been logged, which is why
        # the log can look like it finished successfully even though the
        # job failed). Scale with BAM size, same pattern as the variant callers.
        mem_mb     = lambda wc, attempt: max(4096, attempt * _bam_size_gb(wc) * 2 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/logs/{sample}_junction_counts.log"
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
# Beta-binomial tests vs GTEx (one job per sample × tissue)
# [per-sample configurable via perform_binomial_tests_threads in run config]
# ---------------------------------------------------------------------------
def _gtex_file(tissue):
    return f"{config['gtex_data_dir']}/gtex_{tissue}_jxn_counts.txt"

def _gtex_mem(wc, threads, attempt):
    return max(4096, attempt * threads * (5 if wc.tissue == "brain" else 2) * 1024)


rule perform_binomial_tests:
    input:
        jxn_counts = "{outdir}/junction_analysis/junction_counts/{sample}_splice_junction_counts.tsv",
        gtex_file  = lambda wc: _gtex_file(wc.tissue),
    output:
        all_jxns = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_all_junctions.tsv",
    params:
        annotation = config["annotation"],
        samp_cov   = config["sample_coverage_threshold"],
        gtex_cov   = config["gtex_coverage_threshold"],
        gtex_n     = config["gtex_n_threshold"],
        phasing_thr= config["jxn_phasing_threshold"],
        script     = workflow.basedir + "/scripts/perform_splice_junction_beta_binomial_tests.py",
    threads: lambda wc: _rule_threads(wc, "perform_binomial_tests")
    resources:
        mem_mb     = _gtex_mem,
        runtime    = config["time"],
    log:
        "{outdir}/logs/{sample}_{tissue}_betabinom.log"
    shell:
        """
        mkdir -p $(dirname {output.all_jxns})
        python -u {params.script} \\
            --jxn-info-file           {input.jxn_counts} \\
            --gtexfile                {input.gtex_file} \\
            --outfile                 {output.all_jxns} \\
            --sample-coverage-threshold {params.samp_cov} \\
            --gtex-coverage-threshold   {params.gtex_cov} \\
            --gtex-n-threshold          {params.gtex_n} \\
            --phasing-threshold         {params.phasing_thr} \\
            --annotation-file           {params.annotation} \\
            --threads                   {threads} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Filter to significant outlier junctions
# [per-sample configurable via identify_junction_outliers_threads in run config]
# ---------------------------------------------------------------------------
rule identify_junction_outliers:
    input:
        all_jxns = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_all_junctions.tsv",
    output:
        outliers = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_outlier_junctions.tsv",
    params:
        padj_thr   = config["padj_threshold"],
        dpsi_thr   = config["delta_psi_threshold"],
        script     = workflow.basedir + "/scripts/identify_splice_junction_outliers.py",
    threads: lambda wc: _rule_threads(wc, "identify_junction_outliers")
    resources:
        mem_mb     = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/logs/{sample}_{tissue}_jxn_outliers.log"
    shell:
        """
        python -u {params.script} \\
            --infile             {input.all_jxns} \\
            --outfile            {output.outliers} \\
            --padj-threshold     {params.padj_thr} \\
            --delta-PSI-threshold {params.dpsi_thr} \\
            --threads            {threads} \\
        2>&1 | tee {log}
        """
