
def _gnomad_vcf_list(base, mito_vcf):
    chroms = config["gnomad_chroms"]
    vcfs   = [(str(base) + '/gnomad.genomes.v4.1.sites.' + str(c) + '.vcf.bgz') for c in chroms]
    vcfs.append(mito_vcf)
    return ",".join(vcfs)

def _clnsig_args(wc):
    return " ".join(('"' + str(s) + '"') for s in config["clnsig_filter"])

def _cadd_clnsig_args(wc):
    return " ".join(('"' + str(s) + '"') for s in config["cadd_clnsig_filter"])

def _spliceai_clnsig_args(wc):
    return " ".join(('"' + str(s) + '"') for s in config["spliceai_clnsig_filter"])

def _final_dp_flag(wc):
    v = config.get("final_dp_threshold", "")
    return ('--final-DP-threshold ' + str(v)) if str(v).strip() != "" else ""

def _final_af_flag(wc):
    v = config.get("final_af_threshold", "")
    return ('--final-AF-threshold ' + str(v)) if str(v).strip() != "" else ""


rule _2A_compile_variants:
    input:
        longcallr  = "{outdir}/variant_calling/longcallR/{sample}_longcallR_norm.vcf.gz",
        nanots     = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz",
        clair3     = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz",
        deepvar    = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz",
        bed        = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv          = "{outdir}/variant_calling/compiled_variants/{sample}_compiled_variants.tsv",
    params:
        genome      = config["genome"],
        annotation  = config["annotation"],
        gnomad_vcf           = _gnomad_vcf_list(_resolved_gnomad_base(), _resolved_gnomad_mito_vcf()),
        gnomad_vcf_fallback  = _gnomad_vcf_list(_REMOTE_GNOMAD_BASE, _REMOTE_GNOMAD_MITO_VCF),
        clinvar_vcf          = _resolved_clinvar_vcf(),
        clinvar_vcf_fallback = _REMOTE_CLINVAR_VCF,
        annovar_dir = config["annovar_dir"],
        cadd_script         = lambda wc: config["cadd_script"] if _cadd_use_local() else "",
        cadd_data_dir       = lambda wc: config["cadd_data_dir"] if _cadd_use_local() else "",
        cadd_local_prescored_snv   = config["cadd_local_prescored_snv"],
        cadd_local_prescored_indel = config["cadd_local_prescored_indel"],
        cadd_prescored_url  = config["cadd_prescored_url"],
        spliceai_annotation = config["spliceai_annotation"],
        spliceai_prescored_snv   = config["spliceai_prescored_snv_vcf"],
        spliceai_prescored_indel = config["spliceai_prescored_indel_vcf"],
        spliceai_force_prescored = lambda wc: "--SpliceAI-force-prescored-lookup" if config["spliceai_force_prescored_lookup"] else "",
        cadd_gnomad_af = config["cadd_gnomad_af_threshold"],
        cadd_clnsig    = _cadd_clnsig_args,
        spliceai_gnomad_af = config["spliceai_gnomad_af_threshold"],
        spliceai_clnsig    = _spliceai_clnsig_args,
        outprefix   = "{outdir}/variant_calling/compiled_variants/{sample}",
        script      = workflow.basedir + "/scripts/compile_variants.py",
        conda_env_compile_variants = config["conda_env_compile_variants"],
    threads: lambda wc: _rule_threads(wc, "compile_variants")
    resources:
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * threads * 16 * 1024),
        runtime    = config["time"],
    log:
    shell:


rule _2B_filter_variants:
    input:
        tsv = "{outdir}/variant_calling/compiled_variants/{sample}_compiled_variants.tsv",
    output:
        filtered_tsv = "{outdir}/variant_calling/compiled_variants/{sample}_filtered_variants.tsv",
    params:
        clnsig       = _clnsig_args,
        gnomad_af    = config["gnomad_af_threshold"],
        cadd_thr     = config["cadd_threshold"],
        spliceai_thr = config["spliceai_threshold"],
        final_dp_flag = _final_dp_flag,
        final_af_flag = _final_af_flag,
        script       = workflow.basedir + "/scripts/filter_variants.py",
        conda_env_compile_variants = config["conda_env_compile_variants"],
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(2048, attempt * 4 * 1024),
        runtime = config["time"],
    log:
    shell:
