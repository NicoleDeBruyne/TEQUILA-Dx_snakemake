# step 7 - gene quantification per sample (counts, coverage, splice-site assignment, StringTie for AMALGAM)

rule _7A_quantify_gene_count:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/gene_quantification/{sample}_gene_count.tsv",
    params:
        script = workflow.basedir + "/scripts/_7_quantify_gene_count_coverage.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 4,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_7A_quantify_gene_count.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample  {wildcards.sample} \\
            --bam     {input.bam} \\
            --bed     {input.bed} \\
            --metric  count \\
            --outfile {output.tsv} \\
        2>&1 | tee {log}
        """


rule _7B_quantify_gene_coverage:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
    output:
        tsv = "{outdir}/gene_quantification/{sample}_gene_coverage.tsv",
    params:
        script = workflow.basedir + "/scripts/_7_quantify_gene_count_coverage.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 4,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_7B_quantify_gene_coverage.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample  {wildcards.sample} \\
            --bam     {input.bam} \\
            --bed     {input.bed} \\
            --metric  coverage \\
            --outfile {output.tsv} \\
        2>&1 | tee {log}
        """


# 7C reads assigned to genes by their splice sites
rule _7C_quantify_gene_by_assignment:
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
        bed = lambda wc: SAMPLES[wc.sample]["bed"],
        gtf = config["annotation"],
    output:
        tsv   = "{outdir}/gene_quantification/{sample}_gene_assignment.tsv",
        stats = "{outdir}/gene_quantification/{sample}_read_outcomes.tsv",
    params:
        script = workflow.basedir + "/scripts/_7C_quantify_gene_by_assignment.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_7C_quantify_gene_by_assignment.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --sample        {wildcards.sample} \\
            --bam           {input.bam} \\
            --bed           {input.bed} \\
            --gtf           {input.gtf} \\
            --outfile       {output.tsv} \\
            --stats-outfile {output.stats} \\
        2>&1 | tee {log}
        """


# 7D StringTie, feeds the AMALGAM transcriptome in step 9
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
        "{outdir}/../logs/{sample}_7D_amalgam_stringtie.log"
    shell:
        """
        mkdir -p $(dirname {output.gtf}) $(dirname {log})
        (
            export PATH="{params.amalgam_env}/bin:$PATH"
            stringtie -o {output.gtf} {input.bam}
            echo -e "\\nFinished running StringTie on {wildcards.sample}."
        ) 2>&1 | tee {log}
        """
