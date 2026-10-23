def _phase_vcf_inputs(wc):
    if SAMPLES[wc.sample].get("wgs_vcf"):
        return {"wgs": str(wc.outdir) + '/variant_calling/wgs/' + str(wc.sample) + '_wgs_phased.vcf.gz'}
    return {
        "nanots":  str(wc.outdir) + '/variant_calling/nanoTS/' + str(wc.sample) + '_nanoTS_norm.vcf.gz',
        "clair3":  str(wc.outdir) + '/variant_calling/clair3_rna/' + str(wc.sample) + '_clair3_rna_norm.vcf.gz',
        "deepvar": str(wc.outdir) + '/variant_calling/deepvariant/' + str(wc.sample) + '_deepvariant_norm.vcf.gz',
    }

def _phase_vcf_args(wc, input):
    if "wgs" in input.keys():
        return "--wgs-vcf " + str(input.wgs)
    return ("--nanoTS-vcf " + str(input.nanots) + " --clair3-vcf " + str(input.clair3)
            + " --deepvariant-vcf " + str(input.deepvar))


rule _3A_phase_reads:
    input:
        unpack(_phase_vcf_inputs),
        bam     = lambda wc: SAMPLES[wc.sample]["bam"],
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
        vcf_args            = _phase_vcf_args,
        script              = workflow.basedir + "/scripts/_3A_phase_reads.py",
    threads: lambda wc: _rule_threads(wc, "phase_reads")
    resources:
        mem_mb = lambda wc, threads, attempt: max(4096, int(attempt * threads * 2 * 1024)),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_3A_phase_reads.log"
    shell:
        """
        mkdir -p {params.phased_dir}
        find {params.phased_dir} -mindepth 1 -maxdepth 1 -type d -exec rm -rf {{}} +
        python -u {params.script} \\
            --bam        {input.bam} \\
            --bed        {input.bed} \\
            {params.vcf_args} \\
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
