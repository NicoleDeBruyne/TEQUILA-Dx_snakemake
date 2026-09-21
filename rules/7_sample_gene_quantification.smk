
rule _7A_quantify_gene_count:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/gene_quantification/{sample}_gene_count.tsv",
    params:
        script = workflow.basedir + "/scripts/quantify_gene_expression_sample.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 4,
        runtime = config["time"],
    log:
    shell:


rule _7B_quantify_gene_coverage:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/gene_quantification/{sample}_gene_coverage.tsv",
    params:
        script = workflow.basedir + "/scripts/quantify_gene_expression_sample.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 4,
        runtime = config["time"],
    log:
    shell:


rule _7C_quantify_gene_by_assignment:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
        gtf = config["annotation"],
    output:
        tsv   = "{outdir}/gene_quantification/{sample}_gene_assignment.tsv",
        stats = "{outdir}/gene_quantification/{sample}_read_outcomes.tsv",
    params:
        script = workflow.basedir + "/scripts/quantify_gene_by_assignment_sample.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
    shell:


rule _7D_amalgam_stringtie:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
    output:
        gtf = "{outdir}/gene_quantification/amalgam/{sample}_stringtie.gtf",
    params:
        amalgam_env = config["amalgam_dir"] + "/conda_env",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
    shell:
