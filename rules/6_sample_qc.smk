"""
rules/6_sample_qc.smk
Per-sample QC metrics, computed once per sample directly from that sample's
own BAM (+ BED, + GTF for _6C) -- no dependency on cohort membership, so
adding/removing a sample from a cohort, or changing an alias, never
reruns these. The cohort-level merge step (rules/9_merge_results.smk's
_9G/_9H/_9J) combines every sample's own TSV here into the cohort matrices
and plots.

Three rules, one per QC analysis:
  1. _6A_get_on_target_rate  -- total/mapped/on-target read counts, via
                                scripts/get_on_target_rate.py. Mapping rate
                                and on-target rate are cheap ratios derived
                                at merge time, not computed here.
  2. _6B_get_read_attributes -- read-length five-number summary (min, Q1,
                                median, Q3, max) per target type
                                (on_target/off_target/mapped/unmapped), via
                                scripts/get_read_attributes.py. Only the
                                five-number summary is kept, not every
                                individual read length -- see that script's
                                module docstring.
  3. _6C_get_full_length_ratio -- per-gene "full-length ratio" (FLR): the
                                average, across every read overlapping that
                                gene's canonical transcript in this sample's
                                own BAM, of the fraction of the transcript's
                                exonic positions the read's alignment spans
                                via a non-N CIGAR operation. Flags samples/
                                genes with unusually truncated read coverage
                                (degraded RNA, partial-length library prep,
                                etc.), via scripts/get_full_length_ratio_sample.py.
"""

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
        "{outdir}/../logs/{sample}_on_target_rate.log"
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
        script = workflow.basedir + "/scripts/get_read_attributes.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_read_attributes.log"
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


rule _6C_get_full_length_ratio:
    input:
        # Reads this sample's own original BAM directly (fetch()-ing each
        # gene's canonical-transcript region) rather than phase_reads.py's
        # per-gene bulk BAM files -- see scripts/get_full_length_ratio_sample.py's
        # module docstring for why. No dependency on phase_reads.py.
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
        "{outdir}/../logs/{sample}_full_length_ratio.log"
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
