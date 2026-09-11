"""
rules/8_cohort_qc.smk
Cohort-wide QC, run once per BED panel -- pooling every sample on that panel
across all its sample_types (unlike rules/6_merge_hits.smk and
rules/7_cohort_junction_analysis.smk, which both operate per (bed,
sample_type) group). Combines what used to be two separate rule files
(sequencing QC plots + sample-type validation) into one, since both operate
at the same bed-wide scope. See docs/rules/8_cohort_qc.md for details.

Four rules:
  1. _8A_get_on_target_rates -- per-sample mapping rate and on-target rate
                                (fraction of mapped reads falling in that
                                sample's own BED panel), via
                                scripts/plot_on_target_rates.py.
  2. _8B_get_read_attributes -- per-sample read length distributions, split
                                by on-target / off-target / unmapped status,
                                via scripts/plot_read_attributes.py. Computed
                                directly from each BAM's primary alignments
                                (no FASTQ input required -- see the script's
                                docstring for the approximation this makes).
  3. _8C1_build_group_junction_matrix / _8C2_validate_sample_types -- builds
                                each (bed, sample_type) group's junction
                                count matrix, then pools every group sharing
                                this BED panel and compares them against
                                GTEx reference tissues (distance heatmap +
                                PCA), via scripts/validate_sample_type.py.
                                Requesting _8C2's output pulls _8C1 along
                                with it.
  4. _8D_get_full_length_ratio -- per-gene, per-sample "full-length ratio"
                                (FLR): the average, across every read
                                overlapping that gene's canonical
                                transcript in the sample's own BAM, of the
                                fraction of the transcript's exonic
                                positions the read's alignment spans via a
                                non-N CIGAR operation. Flags samples/genes
                                with unusually truncated read coverage
                                (degraded RNA, partial-length library prep,
                                etc.), via scripts/get_full_length_ratio.py.
                                Reads directly from each sample's own BAM
                                with one indexed fetch() per gene against a
                                single open file handle per sample (NOT
                                phase_reads.py's per-gene bulk BAM files --
                                see the script's own docstring for why).
                                Genes are split by median (across samples)
                                read count (config["full_length_ratio_min_reads"],
                                default 100) into two stacked subplots
                                within one figure -- low-read-count genes
                                have noisy avgFLR, so they get their own
                                panel instead of showing up as blank/
                                unreliable cells scattered through one big
                                heatmap. Both panels also carry a per-gene
                                median-read-count bar (log-scaled), a
                                per-sample avgFLR boxplot, and a per-sample
                                total-read-count bar (log-scaled, summed
                                across every gene, not just the ones in
                                that particular panel) -- the boxplot
                                points and the read-count bars are colored
                                by sample_type, same colors as
                                _8C2_validate_sample_types.

_8A/_8B/_8D each build their own "name, bam[, ...]" mapping file from
bed_samples(wc.cohort_id, wc.bed_id) inline in the shell command (mirroring
rules/7_cohort_junction_analysis.smk's cohort_mapping construction) rather
than depending on a Snakemake input file, since the mapping file's content
is fully determined by the run config already in memory. _8D's mapping
file also includes sample_type and a pre-resolved color column (via the
Snakefile's own sample_type_color(), same colors _8C2 uses -- see
scripts/get_full_length_ratio.py).
"""

# NOTE: output:/log: path templates below use string concatenation, not
# f-strings, to combine a config value with a literal Snakemake wildcard
# placeholder like "{bed_id}" -- an f-string's "{{bed_id}}" escape (to
# produce a literal "{bed_id}") does not survive Snakemake's own rule
# parsing and raises a NameError at load time.
_cohort_outdir = config["output_dir"] + "/{cohort_id}"


rule _8A_get_on_target_rates:
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
        beds = lambda wc: [SAMPLES[s]["bed"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
    output:
        mapping_file = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates_mapping_file.tsv",
        tsv          = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates.tsv",
        mapping_pdf  = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates_mapping.pdf",
        ontarget_pdf = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates_ontarget.pdf",
    params:
        names     = lambda wc: _quoted(bed_samples(wc.cohort_id, wc.bed_id)),
        bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        beds      = lambda wc: _quoted([SAMPLES[s]["bed"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        groups    = lambda wc: _quoted([SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        outprefix = lambda wc: (config["output_dir"] + '/' + str(wc.cohort_id) + '/' + str(wc.bed_id) + '/output/cohort_qc/on_target_rates/' + str(wc.bed_id) + '_on_target_rates'),
        title     = lambda wc: config.get("cohort_qc_title", (str(wc.cohort_id) + ' ' + str(wc.bed_id) + ' on-target rates')),
        script    = workflow.basedir + "/scripts/plot_on_target_rates.py",
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "get_on_target_rates", config["threads"])
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(bed_samples(wc.cohort_id, wc.bed_id))),
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


rule _8B_get_read_attributes:
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
        beds = lambda wc: [SAMPLES[s]["bed"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
    output:
        mapping_file       = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_mapping_file.tsv",
        read_attributes    = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_attributes.tsv",
        summary            = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_summary.tsv",
        length_boxplot     = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_lengths_boxplot.pdf",
        length_violin      = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes_read_lengths_violin.pdf",
    params:
        names     = lambda wc: _quoted(bed_samples(wc.cohort_id, wc.bed_id)),
        bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        beds      = lambda wc: _quoted([SAMPLES[s]["bed"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        outprefix = lambda wc: (config["output_dir"] + '/' + str(wc.cohort_id) + '/' + str(wc.bed_id) + '/output/cohort_qc/read_attributes/' + str(wc.bed_id) + '_read_attributes'),
        title     = lambda wc: config.get("cohort_qc_title", (str(wc.cohort_id) + ' ' + str(wc.bed_id) + ' read attributes')),
        script    = workflow.basedir + "/scripts/plot_read_attributes.py",
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "get_read_attributes", config["threads"])
    resources:
        # Reads through every alignment in every BAM (twice, when a bed is
        # given: once for on-target IDs, once for length) -- scale with
        # cohort size like _8A, with a higher floor since this rule is the
        # more memory-hungry of the two (holds every read's length in
        # memory per worker before down-sampling to 100k/group).
        mem_mb  = lambda wc, attempt: attempt * 1024 * 2 * max(8, len(bed_samples(wc.cohort_id, wc.bed_id))),
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


def _group_junction_matrix_inputs(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/junction_counts/' + str(s) + '_junction_count_matrix.tsv')
            for s in GROUPS[group_id]]


rule _8C1_build_group_junction_matrix:
    input:
        matrices = lambda wc: _group_junction_matrix_inputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/junction_count_matrix.tsv",
    params:
        samples = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        script  = workflow.basedir + "/scripts/build_group_junction_matrix.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "build_group_junction_matrix", config["threads"])
    resources:
        # Scales with group size (the matrix has one column per sample).
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/group_junction_matrix.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --infiles      {input.matrices} \\
            --sample-names {params.samples} \\
            --outfile      {output.matrix} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """


rule _8C2_validate_sample_types:
    input:
        query_matrices = lambda wc: [
            (str(group_outdir(gid)) + '/merged_junction_analysis/junction_count_matrix.tsv')
            for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]
        ],
        gtex_matrices = lambda wc: [_gtex_file(t) for t in config["validate_ref_tissues"]],
        bed = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
    output:
        heatmap = _cohort_outdir + "/{bed_id}/output/cohort_qc/validate_sample_types/{bed_id}_distance_heatmap.pdf",
        pca     = _cohort_outdir + "/{bed_id}/output/cohort_qc/validate_sample_types/{bed_id}_PCA.pdf",
    params:
        outprefix    = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/validate_sample_types/' + str(wc.bed_id)),
        ref_names    = lambda wc: _quoted(config["validate_ref_tissues"]),
        ref_colors   = lambda wc: _quoted(config["validate_ref_colors"]),
        query_names  = lambda wc: _quoted([GROUP_SAMPLE_TYPE[gid] for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]]),
        query_colors = lambda wc: _quoted([sample_type_color(GROUP_SAMPLE_TYPE[gid]) for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]]),
        script       = workflow.basedir + "/scripts/validate_sample_type.py",
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "validate_sample_types", 1)
    resources:
        # Pools every sample_type sharing this BED panel; the reference GTEx
        # matrices dominate baseline memory use, hence the high floor -- 256GB
        # covers that regardless of cohort size, and scales up further only
        # once the local cohort itself exceeds that many samples.
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(256, len(bed_samples(wc.cohort_id, wc.bed_id))),
        runtime = 1440,   # 1 day, matches the old bash wrapper's --time=1-00:00:00
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_validate_sample_types.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --matrix-refs  {input.gtex_matrices} \\
            --ref-names    {params.ref_names} \\
            --ref-colors   {params.ref_colors} \\
            --matrix-query {input.query_matrices} \\
            --query-names  {params.query_names} \\
            --query-colors {params.query_colors} \\
            --bed          {input.bed} \\
            --outprefix    {params.outprefix} \\
        2>&1 | tee {log}
        """


rule _8D_get_full_length_ratio:
    input:
        # Reads each sample's own original BAM directly (fetch()-ing each
        # gene's canonical-transcript region) rather than phase_reads.py's
        # per-gene bulk BAM files -- see scripts/get_full_length_ratio.py's
        # module docstring for why. No dependency on phase_reads.py.
        bams = lambda wc: [SAMPLES[s]["bam"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
        bed  = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
        gtf  = config["annotation"],
    output:
        mapping_file = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_mapping_file.tsv",
        matrix       = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_matrix.tsv",
        read_counts  = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_read_counts_matrix.tsv",
        heatmap      = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_heatmap.pdf",
    params:
        names        = lambda wc: _quoted(bed_samples(wc.cohort_id, wc.bed_id)),
        bams         = lambda wc: _quoted([SAMPLES[s]["bam"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        sample_types = lambda wc: _quoted([SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        # Same colors _8C2_validate_sample_types uses -- resolved here via
        # the Snakefile's own sample_type_color() and passed through as a
        # plain column, rather than the script inventing its own palette.
        sample_colors = lambda wc: _quoted([sample_type_color(SAMPLES[s]["sample_type"]) for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        outprefix = lambda wc: (config["output_dir"] + '/' + str(wc.cohort_id) + '/' + str(wc.bed_id) + '/output/cohort_qc/full_length_ratio/' + str(wc.bed_id) + '_full_length_ratio'),
        title     = lambda wc: config.get("cohort_qc_title", (str(wc.cohort_id) + ' ' + str(wc.bed_id) + ' full-length ratio')),
        min_reads = config.get("full_length_ratio_min_reads", 100),
        script    = workflow.basedir + "/scripts/get_full_length_ratio.py",
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "get_full_length_ratio", config["threads"])
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(bed_samples(wc.cohort_id, wc.bed_id))),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_full_length_ratio.log"
    shell:
        """
        mkdir -p $(dirname {output.mapping_file})
        mkdir -p $(dirname {log})

        paste <(printf '%s\\n' {params.names}) \\
              <(printf '%s\\n' {params.bams})  \\
              <(printf '%s\\n' {params.sample_types}) \\
              <(printf '%s\\n' {params.sample_colors}) \\
        > {output.mapping_file}

        python -u {params.script} \\
            --mapping-file {output.mapping_file} \\
            --bed          {input.bed} \\
            --gtf          {input.gtf} \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --min-reads    {params.min_reads} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """
