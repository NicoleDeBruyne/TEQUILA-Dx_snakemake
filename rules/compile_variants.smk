"""
rules/compile_variants.smk
Merges longcallR / NanoTS / Clair3-RNA / DeepVariant VCFs, annotates with
ANNOVAR, gnomAD, ClinVar, CADD, and SpliceAI.
"""

def _gnomad_vcf_list(base, mito_vcf):
    """Build the comma-separated gnomAD VCF list from a given base URL/path
    and mito VCF -- used for both the primary (resolved) and fallback
    (always-canonical-remote) lists below."""
    chroms = config["gnomad_chroms"]
    vcfs   = [f"{base}/gnomad.genomes.v4.1.sites.{c}.vcf.bgz" for c in chroms]
    vcfs.append(mito_vcf)
    return ",".join(vcfs)

def _clnsig_args(wc):
    return " ".join(f'"{s}"' for s in config["clnsig_filter"])

def _spliceai_clnsig_args(wc):
    return " ".join(f'"{s}"' for s in config["spliceai_clnsig_filter"])


rule compile_variants:
    input:
        longcallr  = "{outdir}/variant_calling/longcallR/{sample}_longcallR_norm.vcf.gz",
        nanots     = "{outdir}/variant_calling/nanoTS/{sample}_nanoTS_norm.vcf.gz",
        clair3     = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz",
        deepvar    = "{outdir}/variant_calling/deepvariant/{sample}_deepvariant_norm.vcf.gz",
        bed        = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv          = "{outdir}/variant_calling/compiled_variants/{sample}_compiled_variants.tsv",
        filtered_tsv = "{outdir}/variant_calling/compiled_variants/{sample}_filtered_variants.tsv",
    params:
        genome      = config["genome"],
        annotation  = config["annotation"],
        # Primary: whatever's configured (local path, explicit URL, or the
        # "remote" sentinel resolved to the canonical public URL -- see
        # _resolved_gnomad_base() etc. in the Snakefile). Fallback: always
        # the canonical public URL, so compile_variants.py can retry there
        # if the primary (e.g. a configured local copy) fails at runtime.
        gnomad_vcf           = _gnomad_vcf_list(_resolved_gnomad_base(), _resolved_gnomad_mito_vcf()),
        gnomad_vcf_fallback  = _gnomad_vcf_list(_REMOTE_GNOMAD_BASE, _REMOTE_GNOMAD_MITO_VCF),
        clinvar_vcf          = _resolved_clinvar_vcf(),
        clinvar_vcf_fallback = _REMOTE_CLINVAR_VCF,
        annovar_dir = config["annovar_dir"],
        # cadd_script/cadd_data_dir are only passed through if a local CADD
        # install should be attempted at all (_cadd_use_local() is False when
        # cadd_script is the "remote" sentinel). cadd_local_prescored_snv/
        # cadd_local_prescored_indel and cadd_prescored_url are always
        # passed, so compile_variants.py can fall back to the local
        # pre-scored lookup, then the remote pre-scored SNV lookup, if
        # local scoring isn't configured/requested or fails at runtime
        # (indels are only covered by the local pre-scored file -- the
        # remote one is SNV-only).
        cadd_script         = lambda wc: config["cadd_script"] if _cadd_use_local() else "",
        cadd_data_dir       = lambda wc: config["cadd_data_dir"] if _cadd_use_local() else "",
        cadd_local_prescored_snv   = config["cadd_local_prescored_snv"],
        cadd_local_prescored_indel = config["cadd_local_prescored_indel"],
        cadd_prescored_url  = config["cadd_prescored_url"],
        spliceai_annotation = config["spliceai_annotation"],
        spliceai_prescored_snv   = config["spliceai_prescored_snv_vcf"],
        spliceai_prescored_indel = config["spliceai_prescored_indel_vcf"],
        # store_true flag in compile_variants.py -- either the flag itself or
        # nothing, not a value, so this is a lambda rather than a plain
        # config[...] lookup like the params around it.
        spliceai_force_prescored = lambda wc: "--SpliceAI-force-prescored-lookup" if config["spliceai_force_prescored_lookup"] else "",
        clnsig      = _clnsig_args,
        gnomad_af   = config["gnomad_af_threshold"],
        spliceai_gnomad_af = config["spliceai_gnomad_af_threshold"],
        spliceai_clnsig    = _spliceai_clnsig_args,
        cadd_thr    = config["cadd_threshold"],
        spliceai_thr= config["spliceai_threshold"],
        outprefix   = "{outdir}/variant_calling/compiled_variants/{sample}",
        script      = workflow.basedir + "/scripts/compile_variants.py",
        conda_env_compile_variants = config["conda_env_compile_variants"],
    threads: lambda wc: _rule_threads(wc, "compile_variants")
    resources:
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * threads * 16 * 1024),
        runtime    = config["time"],
    log:
        "{outdir}/logs/{sample}_compile_variants.log"
    shell:
        """
        mkdir -p $(dirname {params.outprefix})
        # This rule uses its own dedicated environment (not the main
        # conda_env every other rule's SLURM job activates via
        # profile/slurm-jobscript.sh) -- see conda_env_compile_variants's
        # comment in config.yaml for why: compile_variants.py needs to run
        # under this env's own Python (pandas/pysam/spliceai/etc.), and
        # directly shells out to several of this env's other tools
        # (ANNOVAR's perl, bcftools/tabix/bgzip, spliceai). Prepending here,
        # same PATH-prepend style as the jobscript itself, so only this
        # rule is affected -- every other rule keeps using conda_env
        # untouched. NOTE: CADD itself does NOT rely on this PATH export --
        # cadd_script points at CADD_wrapper.sh (see config.yaml), which is
        # self-contained and builds its own PATH rather than inheriting
        # this one, specifically so conda_env_compile_variants's own `perl`
        # can't leak into CADD.sh's per-rule conda environments.
        export PATH="{params.conda_env_compile_variants}/bin:$PATH"
        python -u {params.script} \\
            --vcf-files \\
                {input.longcallr} \\
                {input.nanots} \\
                {input.clair3} \\
                {input.deepvar} \\
            --vcf-file-names \\
                {wildcards.sample}_longcallR \\
                {wildcards.sample}_nanoTS \\
                {wildcards.sample}_clair3_rna \\
                {wildcards.sample}_deepvariant \\
            --outprefix {params.outprefix} \\
            --sample-name {wildcards.sample} \\
            --bed {input.bed} \\
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
            --CADD-gnomadAF-threshold {params.gnomad_af} \\
            --CADD-CLNSIG-filter {params.clnsig} \\
            --include-SpliceAI-scores \\
            --SpliceAI-annotation {params.spliceai_annotation} \\
            --SpliceAI-prescored-snv-vcf "{params.spliceai_prescored_snv}" \\
            --SpliceAI-prescored-indel-vcf "{params.spliceai_prescored_indel}" \\
            {params.spliceai_force_prescored} \\
            --SpliceAI-gnomadAF-threshold {params.spliceai_gnomad_af} \\
            --SpliceAI-CLNSIG-filter {params.spliceai_clnsig} \\
            --final-gnomadAF-threshold {params.gnomad_af} \\
            --final-CLNSIG-filter {params.clnsig} \\
            --final-CADD-phred-threshold {params.cadd_thr} \\
            --final-SpliceAI-threshold {params.spliceai_thr} \\
            --threads {threads} \\
        2>&1 | tee {log}
        """
