"""
rules/9_plot_cohort_info.smk
Cohort-wide sequencing QC plots, run once per BED panel -- pooling every
sample on that panel across all its sample_types (unlike rules/6_merge_hits.smk
and rules/7_cohort_junction_analysis.smk, which both operate per (bed,
sample_type) group). See docs/rules/9_plot_cohort_info.md for details.

Split into two independent rules:
  1. _9A_plot_on_target_rates -- per-sample mapping rate and on-target rate
                                 (fraction of mapped reads falling in that
                                 sample's own BED panel), via
                                 scripts/plot_on_target_rates.py.
  2. _9B_plot_read_attributes -- per-sample read length and quality
                                 distributions, split by on-target /
                                 off-target / unmapped status, via
                                 scripts/plot_read_attributes.py. Computed
                                 directly from each BAM's primary alignments
                                 (no FASTQ input required -- see the script's
                                 docstring for the approximation this makes).

Both rules build their own "name, bam, bed[, group]" mapping file from
bed_samples(wc.bed_id) inline in the shell command (mirroring
rules/7_cohort_junction_analysis.smk's cohort_mapping construction) rather
than depending on a Snakemake input file, since the mapping file's content is
fully determined by the run config already in memory.
"""

# NOTE: output:/log: path templates below use string concatenation, not
# f-strings, to combine a config value with a literal Snakemake wildcard
# placeholder like "{bed_id}" -- an f-string's "{{bed_id}}" escape (to
# produce a literal "{bed_id}") does not survive Snakemake's own rule
# parsing and raises a NameError at load time.
_cohort_outdir = config["output_dir"] + "/cohort"


rule _9A_plot_on_target_rates:
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in bed_samples(wc.bed_id)],
        beds = lambda wc: [SAMPLES[s]["bed"] for s in bed_samples(wc.bed_id)],
    output:
        mapping_file = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates_mapping_file.tsv",
        tsv          = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates.tsv",
        mapping_pdf  = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates_mapping.pdf",
        ontarget_pdf = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates_ontarget.pdf",
    params:
        names     = lambda wc: _quoted(bed_samples(wc.bed_id)),
        bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in bed_samples(wc.bed_id)]),
        beds      = lambda wc: _quoted([SAMPLES[s]["bed"] for s in bed_samples(wc.bed_id)]),
        groups    = lambda wc: _quoted([SAMPLES[s]["sample_type"] for s in bed_samples(wc.bed_id)]),
        outprefix = lambda wc: (str(_cohort_outdir) + '/' + str(wc.bed_id) + '/output/cohort_qc/on_target_rates/' + str(wc.bed_id) + '_on_target_rates'),
        title     = lambda wc: config.get("cohort_qc_title", (str(wc.bed_id) + ' on-target rates')),
        script    = workflow.basedir + "/scripts/plot_on_target_rates.py",
    threads: lambda wc: _group_threads(wc.bed_id, "plot_on_target_rates", config["threads"])
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(bed_samples(wc.bed_id))),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_on_target_rates.log"
    shell:
        """
        mkdir -p $(dirname {output.mapping_file})
        mkdir -p $(dirname {log})

        paste <(printf '%s\\n' {params.names}) \\
              <(printf '%s\\n' {params.bams})  \\
              <(printf '%s\\n' {params.beds})  \\
              <(printf '%s\\n' {params.groups}) \\
        > {output.mapping_file}

        python -u {params.script} \\
            --mapping-file {output.mapping_file} \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """


rule _9B_plot_read_attributes:
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in bed_samples(wc.bed_id)],
        beds = lambda wc: [SAMPLES[s]["bed"] for s in bed_samples(wc.bed_id)],
    output:
        mapping_file       = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_mapping_file.tsv",
        read_attributes    = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_attributes.tsv",
        summary            = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_summary.tsv",
        length_boxplot     = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_lengths_boxplot.pdf",
        quality_boxplot    = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_qualities_boxplot.pdf",
        length_violin      = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_lengths_violin.pdf",
        quality_violin     = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_qualities_violin.pdf",
    params:
        names     = lambda wc: _quoted(bed_samples(wc.bed_id)),
        bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in bed_samples(wc.bed_id)]),
        beds      = lambda wc: _quoted([SAMPLES[s]["bed"] for s in bed_samples(wc.bed_id)]),
        outprefix = lambda wc: (str(_cohort_outdir) + '/' + str(wc.bed_id) + '/output/cohort_qc/read_attributes/' + str(wc.bed_id) + '_read_attributes'),
        title     = lambda wc: config.get("cohort_qc_title", (str(wc.bed_id) + ' read attributes')),
        script    = workflow.basedir + "/scripts/plot_read_attributes.py",
    threads: lambda wc: _group_threads(wc.bed_id, "plot_read_attributes", config["threads"])
    resources:
        # Reads through every alignment in every BAM (twice, when a bed is
        # given: once for on-target IDs, once for length/quality) -- scale
        # with cohort size like _9A, with a higher floor since this rule is
        # the more memory-hungry of the two (holds every read's length +
        # quality in memory per worker before down-sampling to 100k/group).
        mem_mb  = lambda wc, attempt: attempt * 2048 * max(8, len(bed_samples(wc.bed_id))),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_read_attributes.log"
    shell:
        """
        mkdir -p $(dirname {output.mapping_file})
        mkdir -p $(dirname {log})

        paste <(printf '%s\\n' {params.names}) \\
              <(printf '%s\\n' {params.bams})  \\
              <(printf '%s\\n' {params.beds})  \\
        > {output.mapping_file}

        python -u {params.script} \\
            --mapping-file {output.mapping_file} \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """
