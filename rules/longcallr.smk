"""
rules/longcallr.smk
Runs longcallR per sample.

Output structure:
    variant_calling/longcallR/work/   — raw tool output (intermediate files)
    variant_calling/longcallR/        — final normalized/compressed/indexed VCF

Rule:
  1. Runs longcallR into work/
  2. Sorts with bcftools sort
  3. Normalizes with bcftools norm -f <genome> -m -both (piped from sort)
  4. Compresses (-Oz) and writes to longcallR/{sample}_longcallR_norm.vcf.gz
  5. Indexes with tabix -f

Snakemake only declares the final _norm.vcf.gz and .tbi as outputs.
"""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _bam(wc):
    return SAMPLES[wc.sample]["bam"]

def _bam_size_gb(wc):
    import os
    bam = SAMPLES[wc.sample]["bam"]
    try:
        return max(1, (os.path.getsize(bam) + 1073741823) // 1073741824)
    except FileNotFoundError:
        return 10  # fallback


# ---------------------------------------------------------------------------
# longcallR
# ---------------------------------------------------------------------------
rule longcallr:
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
        mem_mb  = lambda wc, threads, attempt: max(4096, attempt * threads * _bam_size_gb(wc) * 2 * 1024),
        runtime = config["time"],
    log:
        "{outdir}/logs/{sample}_longcallR.log"
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
