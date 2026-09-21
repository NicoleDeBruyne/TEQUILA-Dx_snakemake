
# Per-sample splice junction usage vs. GTEx: count junctions, fit GTEx as a reference distribution,
# then test each sample's junctions against it and flag outliers
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


# GTEx junction-count file for a tissue, used as the reference cohort for outlier calling
def _gtex_file(tissue):
    return (str(config['gtex_data_dir']) + '/gtex_' + str(tissue) + '_jxn_counts.txt')

_GTEX_BETA_FITS_TEMPLATE = (
    str(config['output_dir']) + '/gtex_beta_distributions/gtex_{tissue}_beta_fits'
    '_cov' + str(config['gtex_coverage_threshold']) +
    '_n' + str(config['gtex_n_threshold']) +
    '_tol' + str(config['psi_rescale_factor']) + '.tsv'
)

def _gtex_fit_mem(wc, threads, attempt):
    return max(4096, attempt * threads * (6 if wc.tissue == "brain" else 2) * 1024)


# Fit a beta distribution to each junction's PSI-approx values across GTEx samples (per tissue)
rule _5B1_fit_gtex_beta_distributions:
    input:
        gtex_file = lambda wc: _gtex_file(wc.tissue),
        bed_files = all_bed_files(),
    output:
        beta_fits = _GTEX_BETA_FITS_TEMPLATE,
    params:
        gtex_cov  = config["gtex_coverage_threshold"],
        gtex_n    = config["gtex_n_threshold"],
        psi_tol   = config["psi_rescale_factor"],
        script    = workflow.basedir + "/scripts/fit_gtex_beta_distributions.py",
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


# Same as _5B1, but for junctions this sample has that aren't already in the GTEx reference
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
        psi_tol   = config["psi_rescale_factor"],
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


# Test this sample's junction usage against the fitted GTEx beta distribution
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
        psi_tol    = config["psi_rescale_factor"],
        phasing_thr= config["jxn_phasing_threshold"],
        script     = workflow.basedir + "/scripts/perform_splice_junction_beta_binomial_tests.py",
    threads: lambda wc: _rule_threads(wc, "perform_binomial_tests")
    resources:
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


# Apply padj/delta-PSI thresholds to flag this sample's outlier junctions
rule _5D_identify_junction_outliers:
    input:
        all_jxns = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_all_junctions.tsv",
    output:
        outliers = "{outdir}/junction_analysis/gtex_{tissue}/{sample}_gtex_{tissue}_outlier_junctions.tsv",
    params:
        padj_thr   = config["psi_padj_threshold"],
        dpsi_thr   = config["delta_psi_threshold"],
        script     = workflow.basedir + "/scripts/identify_splice_junction_outliers.py",
    threads: 1
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
