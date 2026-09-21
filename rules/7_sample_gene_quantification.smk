
# Per-sample gene expression, via 4 independent metrics (count, coverage, splice-site assignment, AMALGAM)
# Number of distinct reads overlapping each panel gene's region
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
        "{outdir}/../logs/{sample}_gene_count.log"
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


# Max per-base pileup depth over each panel gene's region
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
        "{outdir}/../logs/{sample}_gene_coverage.log"
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


# Assign each read to a gene by splice-site/exon overlap with the annotation (more precise than raw overlap)
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
        "{outdir}/../logs/{sample}_gene_assignment.log"
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


# Assemble transcripts with StringTie for AMALGAM-based gene quantification
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
        "{outdir}/../logs/{sample}_amalgam_stringtie.log"
    shell:
        """
        mkdir -p $(dirname {output.gtf}) $(dirname {log})
        (
            export PATH="{params.amalgam_env}/bin:$PATH"
            stringtie -o {output.gtf} {input.bam}
            echo -e "\\nFinished running StringTie on {wildcards.sample}."
        ) 2>&1 | tee {log}
        """
