"""
rules/4_ase_analysis.smk
Runs binomial ASE outlier detection on the per-gene haplotype table written
by phase_reads.py. See docs/rules/4_ase_analysis.md for details.
"""

rule _4A_detect_ase_outliers:
    input:
        infile = "{outdir}/phased_reads/{sample}_phasing_summary.tsv",
    output:
        tsv = "{outdir}/ase_analysis/{sample}_binomial_ase_results.tsv",
    params:
        sample_cov_thr = config["sample_coverage_threshold"],
        padj_thr       = config["padj_threshold"],
        hap_ratio_thr  = config["haplotype_ratio_threshold"],
        phasing_thr    = config["ase_phasing_threshold"],
        outprefix      = "{outdir}/ase_analysis/{sample}_binomial",
        script         = workflow.basedir + "/scripts/detect_ase_outliers.py",
    threads: 1
    resources:
        mem_mb     = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/logs/{sample}_ase_outliers.log"
    shell:
        """
        python -u {params.script} \\
            --infile                    {input.infile} \\
            --sample-coverage-threshold {params.sample_cov_thr} \\
            --padj-threshold            {params.padj_thr} \\
            --haplotype-ratio-threshold {params.hap_ratio_thr} \\
            --phasing-threshold         {params.phasing_thr} \\
            --outprefix                 {params.outprefix} \\
        2>&1 | tee {log}
        """