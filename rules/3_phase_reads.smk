# step 3 - phase reads per gene and write haplotype BAMs

# WGS variants if we have them, otherwise the RNA caller VCFs
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


# 3B/3C only run for samples with a WGS VCF. 3A stays the primary (WGS-aided) phasing that everything
# downstream reads; 3B redoes the phasing from the RNA calls alone so we can check it against WGS.
# Kept in a sibling folder on purpose: 3A wipes every subfolder of phased_reads/ when it runs.
def _rna_phase_vcfs(wc):
    od, s = str(wc.outdir), str(wc.sample)
    return {
        "nanots":  od + '/variant_calling/nanoTS/' + s + '_nanoTS_norm.vcf.gz',
        "clair3":  od + '/variant_calling/clair3_rna/' + s + '_clair3_rna_norm.vcf.gz',
        "deepvar": od + '/variant_calling/deepvariant/' + s + '_deepvariant_norm.vcf.gz',
    }


rule _3B_phase_reads_wgs_blinded:
    input:
        unpack(_rna_phase_vcfs),
        bam     = lambda wc: SAMPLES[wc.sample]["bam"],
        bed     = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        ase_infile = "{outdir}/phased_reads_wgs_blinded/{sample}_phasing_summary.tsv",
        mapping    = "{outdir}/phased_reads_wgs_blinded/{sample}_gene_bam_mapping_file.tsv",
    params:
        genome              = config["genome"],
        phased_dir          = "{outdir}/phased_reads_wgs_blinded",
        min_dp              = config["phasing_min_dp"],
        min_af              = config["phasing_min_af"],
        phasing_thr         = config["phasing_threshold"],
        terminal_prop       = config["terminal_variant_proportion"],
        min_dist            = config["min_dist_from_read_end_variant_phasing"],
        script              = workflow.basedir + "/scripts/_3A_phase_reads.py",
    threads: lambda wc: _rule_threads(wc, "phase_reads")
    resources:
        mem_mb = lambda wc, threads, attempt: max(4096, int(attempt * threads * 2 * 1024)),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_3B_phase_reads_wgs_blinded.log"
    shell:
        """
        mkdir -p {params.phased_dir}
        find {params.phased_dir} -mindepth 1 -maxdepth 1 -type d -exec rm -rf {{}} +
        python -u {params.script} \\
            --bam        {input.bam} \\
            --bed        {input.bed} \\
            --nanoTS-vcf {input.nanots} \\
            --clair3-vcf {input.clair3} \\
            --deepvariant-vcf {input.deepvar} \\
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


# truth is the 1E panel VCF (plain WGS phase), not 1F, which already has the NanoTS blocks merged in.
# the WGS-aided run is checked too, mostly as a sanity check
rule _3C_check_phasing_vs_wgs:
    input:
        wgs_vcf = "{outdir}/variant_calling/wgs/{sample}_wgs_panel.vcf.gz",
        blinded = "{outdir}/phased_reads_wgs_blinded/{sample}_gene_bam_mapping_file.tsv",
        aided   = "{outdir}/phased_reads/{sample}_gene_bam_mapping_file.tsv",
    output:
        sites = "{outdir}/phasing_vs_wgs/{sample}_phasing_vs_wgs_sites.tsv",
        genes = "{outdir}/phasing_vs_wgs/{sample}_phasing_vs_wgs_genes.tsv",
    params:
        min_depth = config.get("phasing_check_min_depth", 5),
        min_frac  = config.get("phasing_check_min_frac", 0.8),
        outprefix = "{outdir}/phasing_vs_wgs/{sample}_phasing_vs_wgs",
        script    = workflow.basedir + "/scripts/_3C_check_phasing_vs_wgs.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4 * 1024,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_3C_check_phasing_vs_wgs.log"
    shell:
        """
        python -u {params.script} \\
            --sample    {wildcards.sample} \\
            --wgs-vcf   {input.wgs_vcf} \\
            --mode      wgs_blinded={input.blinded} wgs_aided={input.aided} \\
            --min-depth {params.min_depth} \\
            --min-frac  {params.min_frac} \\
            --outprefix {params.outprefix} \\
        2>&1 | tee {log}
        """
