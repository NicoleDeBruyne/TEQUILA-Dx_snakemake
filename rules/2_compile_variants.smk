# step 2 - merge the caller VCFs, annotate, filter down to candidate variants

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

# WGS VCF replaces the RNA callers when we have one
def _has_wgs(wc):
    return bool(SAMPLES[wc.sample].get("wgs_vcf"))

def _compile_vcfs(wc):
    if _has_wgs(wc):
        return [str(wc.outdir) + '/variant_calling/wgs/' + str(wc.sample) + '_wgs_phased.vcf.gz']
    return [
        str(wc.outdir) + '/variant_calling/longcallR/' + str(wc.sample) + '_longcallR_norm.vcf.gz',
        str(wc.outdir) + '/variant_calling/nanoTS/' + str(wc.sample) + '_nanoTS_norm.vcf.gz',
        str(wc.outdir) + '/variant_calling/clair3_rna/' + str(wc.sample) + '_clair3_rna_norm.vcf.gz',
        str(wc.outdir) + '/variant_calling/deepvariant/' + str(wc.sample) + '_deepvariant_norm.vcf.gz',
    ]

def _compile_vcf_names(wc):
    callers = ["wgs"] if _has_wgs(wc) else ["longcallR", "nanoTS", "clair3_rna", "deepvariant"]
    return " ".join(str(wc.sample) + '_' + c for c in callers)

def _compile_bed_padding(wc):
    return config.get("wgs_panel_padding", 10000) if _has_wgs(wc) else 0

def _compile_vaf_flag(wc):
    return "--VAF-from-AD" if _has_wgs(wc) else ""

# blank in the config = don't pass the flag at all
def _final_dp_flag(wc):
    v = config.get("final_dp_threshold", "")
    return ('--final-DP-threshold ' + str(v)) if str(v).strip() != "" else ""

def _final_af_flag(wc):
    v = config.get("final_af_threshold", "")
    return ('--final-AF-threshold ' + str(v)) if str(v).strip() != "" else ""


# 2A merge + annotate (gnomAD, ClinVar, CADD, SpliceAI)
rule _2A_compile_variants:
    input:
        vcfs       = _compile_vcfs,
        bed        = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv          = "{outdir}/variant_calling/compiled_variants/{sample}_compiled_variants.tsv",
    params:
        genome      = config["genome"],
        annotation  = config["annotation"],
        vcf_names   = _compile_vcf_names,
        bed_padding = _compile_bed_padding,
        vaf_flag    = _compile_vaf_flag,
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
        script      = workflow.basedir + "/scripts/_2A_compile_variants.py",
        conda_env_compile_variants = config["conda_env_compile_variants"],
    threads: lambda wc: _rule_threads(wc, "compile_variants")
    resources:
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * threads * 16 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/../logs/{sample}_2A_compile_variants.log"
    shell:
        """
        mkdir -p $(dirname {params.outprefix})
        export PATH="{params.conda_env_compile_variants}/bin:$PATH"
        python -u {params.script} \\
            --vcf-files {input.vcfs} \\
            --vcf-file-names {params.vcf_names} \\
            --outprefix {params.outprefix} \\
            --sample-name {wildcards.sample} \\
            --bed {input.bed} \\
            --bed-padding {params.bed_padding} \\
            {params.vaf_flag} \\
            --genome {params.genome} \\
            --gtf {params.annotation} \\
            --gnomad-vcf {params.gnomad_vcf} \\
            --gnomad-vcf-fallback {params.gnomad_vcf_fallback} \\
            --clinvar-vcf {params.clinvar_vcf} \\
            --clinvar-vcf-fallback {params.clinvar_vcf_fallback} \\
            --ANNOVAR "{params.annovar_dir}" \\
            --CADD-script "{params.cadd_script}" \\
            --CADD-data-dir "{params.cadd_data_dir}" \\
            --CADD-local-prescored-snv "{params.cadd_local_prescored_snv}" \\
            --CADD-local-prescored-indel "{params.cadd_local_prescored_indel}" \\
            --CADD-prescored-url {params.cadd_prescored_url} \\
            --CADD-gnomadAF-threshold {params.cadd_gnomad_af} \\
            --CADD-CLNSIG-filter {params.cadd_clnsig} \\
            --include-SpliceAI-scores \\
            --SpliceAI-annotation {params.spliceai_annotation} \\
            --SpliceAI-prescored-snv-vcf "{params.spliceai_prescored_snv}" \\
            --SpliceAI-prescored-indel-vcf "{params.spliceai_prescored_indel}" \\
            {params.spliceai_force_prescored} \\
            --SpliceAI-gnomadAF-threshold {params.spliceai_gnomad_af} \\
            --SpliceAI-CLNSIG-filter {params.spliceai_clnsig} \\
            --threads {threads} \\
        2>&1 | tee {log}
        """


# 2B keep rare, potentially damaging stuff
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
        script       = workflow.basedir + "/scripts/_2B_filter_variants.py",
        conda_env_compile_variants = config["conda_env_compile_variants"],
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4 * 1024,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_2B_filter_variants.log"
    shell:
        """
        export PATH="{params.conda_env_compile_variants}/bin:$PATH"
        python -u {params.script} \\
            --compiled-tsv {input.tsv} \\
            --outfile {output.filtered_tsv} \\
            --final-gnomadAF-threshold {params.gnomad_af} \\
            --final-CLNSIG-filter {params.clnsig} \\
            --final-CADD-phred-threshold {params.cadd_thr} \\
            --final-SpliceAI-threshold {params.spliceai_thr} \\
            {params.final_dp_flag} \\
            {params.final_af_flag} \\
        2>&1 | tee {log}
        """
