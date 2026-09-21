
# Per-sample allele-specific expression: binomial test on hap1 vs hap2 read counts per gene
rule _4A_detect_ase_outliers:
    input:
        infile = "{outdir}/phased_reads/{sample}_phasing_summary.tsv",
    output:
        tsv      = "{outdir}/ase_analysis/{sample}_binomial_ase_results.tsv",
        outliers = "{outdir}/ase_analysis/{sample}_binomial_ase_outliers.tsv",
    params:
        sample_cov_thr = config["sample_coverage_threshold"],
        padj_thr       = config["ase_padj_threshold"],
        minor_hap_freq_thr = config["minor_haplotype_frequency_threshold"],
        phasing_thr    = config["ase_phasing_threshold"],
        outprefix      = "{outdir}/ase_analysis/{sample}_binomial",
        script         = workflow.basedir + "/scripts/detect_ase_outliers.py",
    threads: 1
    resources:
        mem_mb     = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_ase_outliers.log"
    shell:
        """
        python -u {params.script} \\
            --infile                    {input.infile} \\
            --sample-coverage-threshold {params.sample_cov_thr} \\
            --padj-threshold            {params.padj_thr} \\
            --minor-haplotype-frequency-threshold {params.minor_hap_freq_thr} \\
            --phasing-threshold         {params.phasing_thr} \\
            --outprefix                 {params.outprefix} \\
        2>&1 | tee {log}
        """
