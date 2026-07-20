"""
rules/clair3_rna.smk
Runs Clair3-RNA per sample.

Output structure:
    variant_calling/clair3_rna/work/   — raw tool output (intermediate files)
    variant_calling/clair3_rna/        — final normalized/compressed/indexed VCF

Rule:
  1. Runs Clair3-RNA into work/
  2. Sorts with bcftools sort
  3. Normalizes with bcftools norm -f <genome> -m -both (piped from sort)
  4. Compresses (-Oz) and writes to clair3_rna/{sample}_clair3_rna_norm.vcf.gz
  5. Indexes with tabix -f

Snakemake only declares the final _norm.vcf.gz and .tbi as outputs.
"""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
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
        return 10  # fallback

def _bind_dirs(*paths):
    """Return a comma-separated, de-duplicated list of directories to bind into a
    singularity container (singularity exec -B), derived from actual file/directory
    paths this specific rule needs -- rather than a fixed, environment-specific list
    of paths in config.yaml, which breaks the moment the pipeline runs somewhere the
    data doesn't happen to live at those exact locations."""
    import os
    dirs = []
    for p in paths:
        d = p if os.path.isdir(p) else os.path.dirname(os.path.abspath(p))
        if d and d not in dirs:
            dirs.append(d)
    return ",".join(dirs)


# ---------------------------------------------------------------------------
# Clair3-RNA
# ---------------------------------------------------------------------------
rule clair3_rna:
    input:
        bam = _bam,
    output:
        vcf_gz  = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz",
        vcf_tbi = "{outdir}/variant_calling/clair3_rna/{sample}_clair3_rna_norm.vcf.gz.tbi",
    params:
        genome   = config["genome"],
        image    = config["clair3_rna_image"],
        binds    = lambda wc: _bind_dirs(_bam(wc), config["genome"], _outdir(wc), config["conda_env"]),
        whatshap = lambda wc: f"{config['conda_env']}/bin/whatshap",
        work_dir = "{outdir}/variant_calling/clair3_rna/work",
    threads: lambda wc: _rule_threads(wc, "clair3_rna")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(4096, attempt * threads * _bam_size_gb(wc) * 1 * 1024),
        runtime = config["time"],
    log:
        "{outdir}/logs/{sample}_clair3_rna.log"
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
