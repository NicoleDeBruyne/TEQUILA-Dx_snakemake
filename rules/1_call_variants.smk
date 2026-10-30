# step 1 - small variant calling on the RNA BAM, plus pulling in / phasing a WGS VCF if there is one

def _bam(wc):
    return SAMPLES[wc.sample]["bam"]

def _outdir(wc):
    return SAMPLES[wc.sample]["outdir"]

def _bam_size_gb(wc):
    import os
    bam = SAMPLES[wc.sample]["bam"]
    try:
        return max(1, (os.path.getsize(bam) + 1073741823) // 1073741824)
    except FileNotFoundError:
        return 10

# singularity needs every directory we touch bind-mounted
def _bind_dirs(files=(), dirs=()):
    import os
    out = [os.path.dirname(os.path.abspath(f)) for f in files]
    out += [os.path.abspath(d) for d in dirs]
    return ",".join(dict.fromkeys(d for d in out if d))


# 1A NanoTS
rule _1A_nanots:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz.tbi",
    params:
        genome        = config["genome"],
        image         = config["nanots_image"],
        binds         = lambda wc: _bind_dirs(
                            files=[_bam(wc), config["genome"],
                                   config["nanots_model_unphased"], config["nanots_model_phased"]],
                            dirs=[_outdir(wc)]),
        model_unphased= config["nanots_model_unphased"],
        model_phased  = config["nanots_model_phased"],
        work_dir      = "{outdir}/variant_calling/nanoTS/work",
    threads: lambda wc: _rule_threads(wc, "nanots")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(4096, attempt * threads * _bam_size_gb(wc) * 1 * 1024),
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_1A_nanots.log"
    shell:
        """
        mkdir -p {params.work_dir}
        singularity exec -B {params.binds} \\
            {params.image} nanoTS full_pipeline \\
            --bam {input.bam} \\
            --ref {params.genome} \\
            --threads {threads} \\
            --model_unphased {params.model_unphased} \\
            --model_phased {params.model_phased} \\
            --outdir {params.work_dir} \\
        2>&1 | tee {log}
        bcftools sort {params.work_dir}/phased_predict.pass.vcf \
            | bcftools norm -f {params.genome} -m -both -O z -o {output.vcf_gz} \
        2>&1 | tee -a {log}
        tabix -f -p vcf {output.vcf_gz}
        rm -f {params.work_dir}/suffix_qname.pgbam
        rm -f {params.work_dir}/suffix_qname.pgbai
        rm -f {params.work_dir}/suffix_qname.bam
        rm -f {params.work_dir}/suffix_qname.bam.bai
        """


# 1B longcallR
rule _1B_longcallr:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/longcallR/{sample}_longcallR_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/longcallR/{sample}_longcallR_norm.vcf.gz.tbi",
    params:
        genome    = config["genome"],
        longcallr = config["longcallr_bin"],
        work_dir  = "{outdir}/variant_calling/longcallR/work",
    threads: lambda wc: _rule_threads(wc, "longcallr")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(4096, int(attempt * threads * _bam_size_gb(wc) * 2.5 * 1024)),
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_1B_longcallr.log"
    shell:
        """
        mkdir -p {params.work_dir}
        {params.longcallr} \\
            --bam-path {input.bam} \\
            --ref-path {params.genome} \\
            --threads {threads} \\
            --preset ont-cdna \\
            --no-bam-output \\
            --output {params.work_dir}/{wildcards.sample}_longcallR \\
        2>&1 | tee {log}
        bcftools sort {params.work_dir}/{wildcards.sample}_longcallR.vcf \
            | bcftools norm -f {params.genome} -m -both -O z -o {output.vcf_gz} \
        2>&1 | tee -a {log}
        tabix -f -p vcf {output.vcf_gz}
        """


# 1C Clair3-RNA
rule _1C_clair3_rna:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz.tbi",
    params:
        genome   = config["genome"],
        image    = config["clair3_rna_image"],
        binds    = lambda wc: _bind_dirs(
                       files=[_bam(wc), config["genome"]],
                       dirs=[_outdir(wc), config["conda_env"]]),
        whatshap = lambda wc: (str(config['conda_env']) + '/bin/whatshap'),
        work_dir = "{outdir}/variant_calling/clair3_rna/work",
    threads: lambda wc: _rule_threads(wc, "clair3_rna")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(4096, attempt * threads * _bam_size_gb(wc) * 1 * 1024),
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_1C_clair3_rna.log"
    shell:
        """
        mkdir -p {params.work_dir}
        singularity exec -B {params.binds} \\
            {params.image} /opt/bin/run_clair3_rna \\
            --bam_fn {input.bam} \\
            --ref_fn {params.genome} \\
            --threads {threads} \\
            --platform ont_r10_dorado_cdna \\
            --output_dir {params.work_dir} \\
            --enable_phasing_model \\
            --whatshap {params.whatshap} \\
            --conda_prefix /opt/conda/envs/clair3_rna \\
            --sample_name {wildcards.sample}_clair3_rna \\
            --output_prefix {wildcards.sample}_clair3_rna \\
        2>&1 | tee {log}
        bcftools sort {params.work_dir}/{wildcards.sample}_clair3_rna_enable_phasing.vcf.gz \
            | bcftools norm -f {params.genome} -m -both -O z -o {output.vcf_gz} \
        2>&1 | tee -a {log}
        tabix -f -p vcf {output.vcf_gz}
        rm -rf {params.work_dir}/tmp
        """


# 1D DeepVariant (long-read RNA model)
rule _1D_deepvariant:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz.tbi",
    params:
        genome   = config["genome"],
        image    = config["deepvariant_image"],
        binds    = lambda wc: _bind_dirs(
                       files=[_bam(wc), config["genome"]],
                       dirs=[_outdir(wc)]),
        work_dir = "{outdir}/variant_calling/deepvariant/work",
    threads: lambda wc: _rule_threads(wc, "deepvariant")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(8192, attempt * threads * _bam_size_gb(wc) * 1 * 1024),
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_1D_deepvariant.log"
    shell:
        """
        mkdir -p {params.work_dir}
        singularity exec -B {params.binds} \\
            {params.image} /opt/deepvariant/bin/run_deepvariant \\
            --reads {input.bam} \\
            --ref {params.genome} \\
            --num_shards {threads} \\
            --model_type ONT_R104 \\
            --intermediate_results_dir {params.work_dir} \\
            --output_vcf {params.work_dir}/{wildcards.sample}_deepvariant.vcf.gz \\
            --sample_name {wildcards.sample}_deepvariant \\
        2>&1 | tee {log}
        bcftools sort {params.work_dir}/{wildcards.sample}_deepvariant.vcf.gz \
            | bcftools norm -f {params.genome} -m -both -O z -o {output.vcf_gz} \
        2>&1 | tee -a {log}
        tabix -f -p vcf {output.vcf_gz}
        rm -f {params.work_dir}/make_examples.tfrecord-*.gz
        rm -f {params.work_dir}/make_examples.tfrecord-*.gz.example_info.json
        rm -f {params.work_dir}/call_variants_output-*.tfrecord.gz
        """


# 1E/1F only run for samples that come with a WGS VCF
rule _1E_prepare_wgs_vcf:
    input:
        vcf = lambda wc: SAMPLES[wc.sample]["wgs_vcf"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        norm_vcf  = "{outdir}/variant_calling/wgs/{sample}_wgs_norm.vcf.gz",
        norm_tbi  = "{outdir}/variant_calling/wgs/{sample}_wgs_norm.vcf.gz.tbi",
        panel_vcf = "{outdir}/variant_calling/wgs/{sample}_wgs_panel.vcf.gz",
        panel_tbi = "{outdir}/variant_calling/wgs/{sample}_wgs_panel.vcf.gz.tbi",
    params:
        genome  = config["genome"],
        padding = config.get("wgs_panel_padding", 10000),
        script  = workflow.basedir + "/scripts/_1E_prepare_wgs_vcf.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4 * 1024,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_1E_prepare_wgs_vcf.log"
    shell:
        """
        python -u {params.script} \\
            --vcf       {input.vcf} \\
            --genome    {params.genome} \\
            --bed       {input.bed} \\
            --padding   {params.padding} \\
            --norm-vcf  {output.norm_vcf} \\
            --panel-vcf {output.panel_vcf} \\
        2>&1 | tee {log}
        """


rule _1F_phase_wgs_variants:
    input:
        wgs_vcf    = "{outdir}/variant_calling/wgs/{sample}_wgs_panel.vcf.gz",
        nanots_vcf = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz",
    output:
        vcf_gz  = "{outdir}/variant_calling/wgs/{sample}_wgs_phased.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/wgs/{sample}_wgs_phased.vcf.gz.tbi",
    params:
        script = workflow.basedir + "/scripts/_1F_phase_wgs_variants.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4 * 1024,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_1F_phase_wgs_variants.log"
    shell:
        """
        python -u {params.script} \\
            --wgs-vcf    {input.wgs_vcf} \\
            --nanoTS-vcf {input.nanots_vcf} \\
            --outfile    {output.vcf_gz} \\
        2>&1 | tee {log}
        """
