"""
rules/2_compile_variants.smk
Merges longcallR / NanoTS / Clair3-RNA / DeepVariant VCFs, annotates with
ANNOVAR, gnomAD, ClinVar, CADD, and SpliceAI.
See docs/rules/2_compile_variants.md for the CADD/SpliceAI fallback tiers,
the dedicated conda env, and the CADD_wrapper.sh PATH-isolation details.
"""

def _gnomad_vcf_list(base, mito_vcf):
    """Build the comma-separated gnomAD VCF list from a given base URL/path
    and mito VCF -- used for both the primary and fallback (always-remote) lists."""
    chroms = config["gnomad_chroms"]
    vcfs   = [f"{base}/gnomad.genomes.v4.1.sites.{c}.vcf.bgz" for c in chroms]
    vcfs.append(mito_vcf)
    return ",".join(vcfs)

def _clnsig_args(wc):
    return " ".join(f'"{s}"' for s in config["clnsig_filter"])

def _cadd_clnsig_args(wc):
    return " ".join(f'"{s}"' for s in config["cadd_clnsig_filter"])

def _spliceai_clnsig_args(wc):
    return " ".join(f'"{s}"' for s in config["spliceai_clnsig_filter"])

def _final_dp_flag(wc):
    v = config.get("final_dp_threshold", "")
    return f"--final-DP-threshold {v}" if str(v).strip() != "" else ""

def _final_af_flag(wc):
    v = config.get("final_af_threshold", "")
    return f"--final-AF-threshold {v}" if str(v).strip() != "" else ""


rule _2A_compile_variants:
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
        # Primary: whatever's configured (local, explicit URL, or "remote"
        # resolved to the canonical public URL). Fallback: always the
        # canonical public URL, for compile_variants.py to retry against.
        gnomad_vcf           = _gnomad_vcf_list(_resolved_gnomad_base(), _resolved_gnomad_mito_vcf()),
        gnomad_vcf_fallback  = _gnomad_vcf_list(_REMOTE_GNOMAD_BASE, _REMOTE_GNOMAD_MITO_VCF),
        clinvar_vcf          = _resolved_clinvar_vcf(),
        clinvar_vcf_fallback = _REMOTE_CLINVAR_VCF,
        annovar_dir = config["annovar_dir"],
        # Only passed through if a local CADD install should be attempted
        # (_cadd_use_local() is False when cadd_script is "remote"). The
        # pre-scored paths are always passed as further fallback tiers --
        # see docs/rules/2_compile_variants.md.
        cadd_script         = lambda wc: config["cadd_script"] if _cadd_use_local() else "",
        cadd_data_dir       = lambda wc: config["cadd_data_dir"] if _cadd_use_local() else "",
        cadd_local_prescored_snv   = config["cadd_local_prescored_snv"],
        cadd_local_prescored_indel = config["cadd_local_prescored_indel"],
        cadd_prescored_url  = config["cadd_prescored_url"],
        spliceai_annotation = config["spliceai_annotation"],
        spliceai_prescored_snv   = config["spliceai_prescored_snv_vcf"],
        spliceai_prescored_indel = config["spliceai_prescored_indel_vcf"],
        # store_true flag in compile_variants.py -- either the flag itself
        # or nothing, hence the lambda instead of a plain config[...] lookup.
        spliceai_force_prescored = lambda wc: "--SpliceAI-force-prescored-lookup" if config["spliceai_force_prescored_lookup"] else "",
        clnsig      = _clnsig_args,
        gnomad_af   = config["gnomad_af_threshold"],
        cadd_gnomad_af = config["cadd_gnomad_af_threshold"],
        cadd_clnsig    = _cadd_clnsig_args,
        spliceai_gnomad_af = config["spliceai_gnomad_af_threshold"],
        spliceai_clnsig    = _spliceai_clnsig_args,
        cadd_thr    = config["cadd_threshold"],
        spliceai_thr= config["spliceai_threshold"],
        final_dp_flag = _final_dp_flag,
        final_af_flag = _final_af_flag,
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
        # This rule uses its own dedicated environment (conda_env_compile_variants,
        # not the main conda_env every other rule uses) -- see
        # docs/rules/2_compile_variants.md for why. CADD itself does NOT rely
        # on this PATH export -- cadd_script points at CADD_wrapper.sh,
        # which builds its own PATH so this env's `perl` can't leak into
        # CADD.sh's per-rule conda environments.
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
            --CADD-gnomadAF-threshold {params.cadd_gnomad_af} \\
            --CADD-CLNSIG-filter {params.cadd_clnsig} \\
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
            {params.final_dp_flag} \\
            {params.final_af_flag} \\
            --threads {threads} \\
        2>&1 | tee {log}
        """
