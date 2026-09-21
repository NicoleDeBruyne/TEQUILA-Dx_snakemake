
# Assign reads to haplotypes per gene using trusted het variants, for use in ASE/junction analysis
rule _3A_phase_reads:
    input:
        bam     = lambda wc: SAMPLES[wc.sample]["bam"],
        nanots  = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz",
        clair3  = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz",
        deepvar = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz",
        bed     = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        ase_infile = "{outdir}/phased_reads/{sample}_phasing_summary.tsv",
        mapping    = "{outdir}/phased_reads/{sample}_gene_bam_mapping_file.tsv",
    params:
        genome              = config["genome"],
        phased_dir          = "{outdir}/phased_reads",
        min_dp              = config["phasing_min_dp"],
        min_af              = config["phasing_min_af"],
        phasing_thr         = config["phasing_threshold"],
        terminal_prop       = config["terminal_variant_proportion"],
        min_dist            = config["min_dist_from_read_end_variant_phasing"],
        script              = workflow.basedir + "/scripts/phase_reads.py",
    threads: lambda wc: _rule_threads(wc, "phase_reads")
    resources:
        mem_mb = lambda wc, threads, attempt: max(4096, int(attempt * threads * 1.5 * 1024)),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_phase_reads.log"
    shell:
        """
        mkdir -p $(dirname {output.ase_infile})
        python -u {params.script} \\
            --bam        {input.bam} \\
            --bed        {input.bed} \\
            --nanoTS-vcf       {input.nanots} \\
            --clair3-vcf       {input.clair3} \\
            --deepvariant-vcf  {input.deepvar} \\
            --min-dp     {params.min_dp} \\
            --min-af     {params.min_af} \\
            --genome     {params.genome} \\
            --outdir     {params.phased_dir} \\
            --name       {wildcards.sample} \\
            --threads    {threads} \\
            --phasing-threshold             {params.phasing_thr} \\
            --terminal-variant-proportion   {params.terminal_prop} \\
            --min-distance-from-read-end    {params.min_dist} \\
        2>&1 | tee {log}
        """
