"""
rules/phase_reads.smk
Builds a merged VCF for phasing then runs whatshap-based read phasing per gene.
Depends on all four variant callers completing first.

phase_reads.py also writes the gene/BAM mapping file (indexing its own
per-gene BAM outputs) directly, as a natural continuation of phasing itself,
for get_junction_counts in junction_analysis.smk.
"""

rule build_vcf_for_phasing:
    input:
        longcallr = "{outdir}/variant_calling/longcallR/{sample}_longcallR_norm.vcf.gz",
        nanots    = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz",
        clair3    = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz",
        deepvar   = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz",
    output:
        vcf_gz  = "{outdir}/phased_reads/{sample}_for_phasing.vcf.gz",
        vcf_tbi = "{outdir}/phased_reads/{sample}_for_phasing.vcf.gz.tbi",
    params:
        script    = workflow.basedir + "/scripts/build_vcf_for_phasing.py",
    threads: 1
    resources:
        mem_mb     = lambda wc, attempt: max(4096, attempt * 8 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/logs/{sample}_build_vcf_for_phasing.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf_gz})
        python -u {params.script} \\
            --longcallR-vcf  {input.longcallr} \\
            --nanoTS-vcf     {input.nanots} \\
            --clair3-vcf     {input.clair3} \\
            --deepvariant-vcf {input.deepvar} \\
            --outfile        {output.vcf_gz} \\
            --sample-name    {wildcards.sample} \\
        2>&1 | tee {log}
        tabix -f -p vcf {output.vcf_gz}
        """


rule phase_reads:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        vcf = "{outdir}/phased_reads/{sample}_for_phasing.vcf.gz",
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        ase_infile = "{outdir}/phased_reads/{sample}_phasing_summary.tsv",
        mapping    = "{outdir}/phased_reads/{sample}_gene_bam_mapping_file.tsv",
    params:
        genome         = config["genome"],
        phased_dir     = "{outdir}/phased_reads",
        phasing_thr    = config["phasing_threshold"],
        terminal_prop  = config["terminal_variant_proportion"],
        min_dist       = config["min_dist_from_read_end_variant_phasing"],
        script         = workflow.basedir + "/scripts/phase_reads.py",
    threads: lambda wc: _rule_threads(wc, "phase_reads")
    resources:
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * threads * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/logs/{sample}_phase_reads.log"
    shell:
        """
        mkdir -p $(dirname {output.ase_infile})
        python -u {params.script} \\
            --bam        {input.bam} \\
            --bed        {input.bed} \\
            --vcf        {input.vcf} \\
            --genome     {params.genome} \\
            --outdir     {params.phased_dir} \\
            --name       {wildcards.sample} \\
            --threads    {threads} \\
            --phasing-threshold             {params.phasing_thr} \\
            --terminal-variant-proportion   {params.terminal_prop} \\
            --min-distance-from-read-end    {params.min_dist} \\
        2>&1 | tee {log}
        """