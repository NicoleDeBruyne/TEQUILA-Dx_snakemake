
rule _6A_get_on_target_rate:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/qc/{sample}_on_target.tsv",
    params:
        script = workflow.basedir + "/scripts/get_on_target_rate.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
    shell:


rule _6B_get_read_attributes:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/qc/{sample}_read_attributes.tsv",
    params:
        script = workflow.basedir + "/scripts/get_read_attributes.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
    shell:


rule _6C_get_full_length_ratio:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
        gtf = config["annotation"],
    output:
        tsv = "{outdir}/qc/{sample}_full_length_ratio.tsv",
    params:
        script = workflow.basedir + "/scripts/get_full_length_ratio_sample.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
    shell:
