"""
rules/7_sample_gene_quantification.smk
Per-sample relative gene expression, computed once per sample directly from
that sample's own BAM (+ BED, + GTF for _7C/_7D) -- no dependency on cohort
membership, so adding/removing a sample from a cohort, or changing an
alias, never reruns these. The cohort-level merge step
(rules/9_merge_results.smk's _9K/_9L/_9M/_9N1-_9N5) combines every sample's
own TSV (or, for AMALGAM, GTF) here into the cohort matrices and plots.

Four rules:
  1. _7A_quantify_gene_count    -- read-count-based proxy (number of
                                    distinct reads overlapping each gene's
                                    BED region), via
                                    scripts/quantify_gene_expression_sample.py
                                    --metric count.
  2. _7B_quantify_gene_coverage -- max-coverage-based proxy (peak per-base
                                    pileup depth in each gene's BED region),
                                    via the same script with --metric coverage.
  3. _7C_quantify_gene_by_assignment -- splice-site-sharing-based proxy:
                                    assigns each alignment to whichever GTF
                                    gene (genome-wide, not just genes on the
                                    BED panel) it shares the most annotated
                                    splice sites with, via
                                    scripts/quantify_gene_by_assignment_sample.py.
                                    See that script's module docstring for
                                    the assignment rules.
  4. _7D_amalgam_stringtie       -- step 1 of AMALGAM's own pipeline (see
                                    its README): de-novo transcriptome
                                    assembly per sample, short-read mode/
                                    default settings. The only step of
                                    AMALGAM's sub-pipeline that's genuinely
                                    per-sample with no cohort dependency --
                                    everything downstream (merging every
                                    sample's assembly into one shared
                                    transcriptome, then quantifying against
                                    it) is inherently cohort-level and stays
                                    in rules/9_merge_results.smk's _9N1-_9N5.

_7A and _7B report two things per gene: a raw value (read count, or max
depth) and a CPTM ("counts per target million") value = raw_value / (sum of
every gene's raw value for this one sample) * 1e6. _7C reports the
analogous raw assigned-read count and CPTM, but ALSO writes a raw count for
every GTF gene that received >=1 assigned read in this sample (assignment
isn't restricted to the BED panel, even though CPTM is). CPTM is what the
cohort-level matrices and boxplots primarily report -- see
scripts/quantify_gene_expression_sample.py's module docstring for why it's
normalized against the panel's own total signal, not overall sequencing
depth.
"""

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
