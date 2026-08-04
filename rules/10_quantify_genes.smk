"""
rules/10_quantify_genes.smk
Approximate relative gene expression across a cohort, run once per (bed,
sample_type) group -- one gene expression matrix + one boxplot-per-gene PDF
set per quantification method. See docs/rules/10_quantify_genes.md.

Three rules:
  1. _10A_plot_relative_gene_count    -- read-count-based proxy (number of
                                         distinct reads overlapping each
                                         gene's BED region), via
                                         scripts/quantify_gene_expression.py
                                         --metric count.
  2. _10B_plot_relative_gene_coverage -- max-coverage-based proxy (peak
                                         per-base pileup depth in each
                                         gene's BED region), via the same
                                         script with --metric coverage.
  3. _10C_run_amalgam                 -- placeholder for a future
                                         independent quantification tool
                                         (e.g. featureCounts/Salmon-based).
                                         Currently a no-op: writes a marker
                                         file explaining nothing has run yet,
                                         so downstream rules/all_outputs()
                                         have a real, trackable target.

Both 10A and 10B report two things per gene per sample: a raw value (read
count, or max depth) and a CPTM ("counts per target million") value =
raw_value / (sum of every gene's raw value for that sample, i.e. every gene
on this BED panel) * 1e6. CPTM is what the matrix and boxplots primarily
report -- it's the "relative expression" number, comparable across samples
regardless of depth, and normalized against the panel's own total signal
(not overall sequencing depth) since this is targeted, not whole-transcriptome,
sequencing.

Only each matrix TSV is tracked as a Snakemake output: -- the per-gene
boxplot PDFs (one per gene on the panel, an a-priori-unknown count) are
written as a side effect of the same script run, the same convention
already used for _7B_identify_cohort_junction_outliers' per-metric
heatmaps/boxplots in rules/7_cohort_junction_analysis.smk.
"""

_cohort_outdir = config["output_dir"] + "/cohort"


rule _10A_plot_relative_gene_count:
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]],
        bed  = lambda wc: bed_path(wc.bed_id),
    output:
        mapping_file = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_mapping_file.tsv",
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_raw.tsv",
    params:
        names     = lambda wc: _quoted(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]),
        bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]]),
        outprefix = lambda wc: (_cohort_outdir + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_count/gene_count"),
        title     = lambda wc: config.get("gene_quant_title") or (str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (read count)"),
        script    = workflow.basedir + "/scripts/quantify_gene_expression.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "plot_relative_gene_count", config["threads"])
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_count_quantification.log"
    shell:
        """
        mkdir -p $(dirname {output.mapping_file}) $(dirname {log})

        paste <(printf '%s\\n' {params.names}) \\
              <(printf '%s\\n' {params.bams})  \\
        > {output.mapping_file}

        python -u {params.script} \\
            --mapping-file {output.mapping_file} \\
            --bed          {input.bed} \\
            --metric       count \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """


rule _10B_plot_relative_gene_coverage:
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]],
        bed  = lambda wc: bed_path(wc.bed_id),
    output:
        mapping_file = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_mapping_file.tsv",
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_raw.tsv",
    params:
        names     = lambda wc: _quoted(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]),
        bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]]),
        outprefix = lambda wc: (_cohort_outdir + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_coverage/gene_coverage"),
        title     = lambda wc: config.get("gene_quant_title") or (str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (max coverage)"),
        script    = workflow.basedir + "/scripts/quantify_gene_expression.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "plot_relative_gene_coverage", config["threads"])
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_coverage_quantification.log"
    shell:
        """
        mkdir -p $(dirname {output.mapping_file}) $(dirname {log})

        paste <(printf '%s\\n' {params.names}) \\
              <(printf '%s\\n' {params.bams})  \\
        > {output.mapping_file}

        python -u {params.script} \\
            --mapping-file {output.mapping_file} \\
            --bed          {input.bed} \\
            --metric       coverage \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """


rule _10C_run_amalgam:
    # Placeholder for a future independent-tool-based quantification
    # (e.g. featureCounts/Salmon). Intentionally a no-op for now -- writes a
    # marker file explaining that, rather than actually quantifying
    # anything, so the rule has a real trackable output and downstream
    # rules/all_outputs() can depend on it without special-casing "not
    # implemented yet".
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]],
        bed  = lambda wc: bed_path(wc.bed_id),
    output:
        marker = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/NOT_YET_IMPLEMENTED.txt",
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_amalgam_quantification.log"
    shell:
        """
        mkdir -p $(dirname {output.marker}) $(dirname {log})
        echo "_10C_run_amalgam is a placeholder and does not run any quantification yet." > {output.marker}
        echo "_10C_run_amalgam: no-op placeholder, nothing to do" 2>&1 | tee {log}
        """
