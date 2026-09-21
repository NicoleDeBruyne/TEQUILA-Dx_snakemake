
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

def _bind_dirs(*paths):
    import os
    dirs = []
    for p in paths:
        d = p if os.path.isdir(p) else os.path.dirname(os.path.abspath(p))
        if d and d not in dirs:
            dirs.append(d)
    return ",".join(dirs)


rule _1A_nanots:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz.tbi",
    params:
        genome        = config["genome"],
        image         = config["nanots_image"],
        binds         = lambda wc: _bind_dirs(_bam(wc), config["genome"], _outdir(wc),
                                               config["nanots_model_unphased"], config["nanots_model_phased"]),
        model_unphased= config["nanots_model_unphased"],
        model_phased  = config["nanots_model_phased"],
        work_dir      = "{outdir}/variant_calling/nanoTS/work",
    threads: lambda wc: _rule_threads(wc, "nanots")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(4096, attempt * threads * _bam_size_gb(wc) * 1 * 1024),
        runtime = config["time"],
    log:
    shell:


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
    shell:


rule _1C_clair3_rna:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz.tbi",
    params:
        genome   = config["genome"],
        image    = config["clair3_rna_image"],
        binds    = lambda wc: _bind_dirs(_bam(wc), config["genome"], _outdir(wc), config["conda_env"]),
        whatshap = lambda wc: (str(config['conda_env']) + '/bin/whatshap'),
        work_dir = "{outdir}/variant_calling/clair3_rna/work",
    threads: lambda wc: _rule_threads(wc, "clair3_rna")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(4096, attempt * threads * _bam_size_gb(wc) * 1 * 1024),
        runtime = config["time"],
    log:
    shell:


rule _1D_deepvariant:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz.tbi",
    params:
        genome   = config["genome"],
        image    = config["deepvariant_image"],
        binds    = lambda wc: _bind_dirs(_bam(wc), config["genome"], _outdir(wc)),
        work_dir = "{outdir}/variant_calling/deepvariant/work",
    threads: lambda wc: _rule_threads(wc, "deepvariant")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(8192, attempt * threads * _bam_size_gb(wc) * 1 * 1024),
        runtime = config["time"],
    log:
    shell:
