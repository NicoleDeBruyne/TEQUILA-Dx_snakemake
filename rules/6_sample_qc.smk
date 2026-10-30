# step 6 - per-sample QC, and ranking genes by phenotype

rule _6A_get_on_target_rate:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/qc/{sample}_on_target.tsv",
    params:
        script = workflow.basedir + "/scripts/_6A_get_on_target_rate.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_6A_get_on_target_rate.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample  {wildcards.sample} \\
            --bam     {input.bam} \\
            --bed     {input.bed} \\
            --outfile {output.tsv} \\
        2>&1 | tee {log}
        """


rule _6B_get_read_attributes:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/qc/{sample}_read_attributes.tsv",
    params:
        script = workflow.basedir + "/scripts/_6B_get_read_attributes.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_6B_get_read_attributes.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample  {wildcards.sample} \\
            --bam     {input.bam} \\
            --bed     {input.bed} \\
            --outfile {output.tsv} \\
        2>&1 | tee {log}
        """


# 6C how full-length the reads are per gene
rule _6C_get_full_length_ratio:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
        gtf = config["annotation"],
    output:
        tsv = "{outdir}/qc/{sample}_full_length_ratio.tsv",
    params:
        script = workflow.basedir + "/scripts/_6C_get_full_length_ratio.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_6C_get_full_length_ratio.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample  {wildcards.sample} \\
            --bam     {input.bam} \\
            --bed     {input.bed} \\
            --gtf     {input.gtf} \\
            --outfile {output.tsv} \\
        2>&1 | tee {log}
        """


rule _6D_get_avg_coverage:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
        gtf = config["annotation"],
    output:
        tsv = "{outdir}/qc/{sample}_avg_coverage.tsv",
    params:
        script = workflow.basedir + "/scripts/_6D_get_avg_coverage.py",
        min_depth = config["merge_min_dp_snv"],
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_6D_get_avg_coverage.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample  {wildcards.sample} \\
            --bam     {input.bam} \\
            --bed     {input.bed} \\
            --gtf     {input.gtf} \\
            --min-depth {params.min_depth} \\
            --outfile {output.tsv} \\
        2>&1 | tee {log}
        """


# 6E AMELIE ranking from the HPO terms
rule _6E_rank_genes_amelie:
    input:
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv   = "{outdir}/phenotype/{sample}_amelie.tsv",
        terms = "{outdir}/phenotype/{sample}_hpo_terms.tsv",
    params:
        hpo     = lambda wc: ";".join(sample_hpo(wc.sample)),
        hpo_obo = hpo_obo_path(),
        timeout = config.get("amelie_timeout", 600),
        retries = config.get("amelie_retries", 3),
        script  = workflow.basedir + "/scripts/_6E_rank_genes_amelie.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 2,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_6E_rank_genes_amelie.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample  {wildcards.sample} \\
            --hpo     "{params.hpo}" \\
            --hpo-obo {params.hpo_obo} \\
            --terms-outfile {output.terms} \\
            --bed     {input.bed} \\
            --timeout {params.timeout} \\
            --retries {params.retries} \\
            --outfile {output.tsv} \\
        2>&1 | tee {log}
        """
