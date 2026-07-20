"""
rules/deepvariant.smk
Runs DeepVariant per sample.

Output structure:
    variant_calling/deepvariant/work/   — raw tool output (intermediate files)
    variant_calling/deepvariant/        — final normalized/compressed/indexed VCF

Rule:
  1. Runs DeepVariant into work/
  2. Sorts with bcftools sort
  3. Normalizes with bcftools norm -f <genome> -m -both (piped from sort)
  4. Compresses (-Oz) and writes to deepvariant/{sample}_deepvariant_norm.vcf.gz
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
# DeepVariant
# ---------------------------------------------------------------------------
rule deepvariant:
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
        "{outdir}/logs/{sample}_deepvariant.log"
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
