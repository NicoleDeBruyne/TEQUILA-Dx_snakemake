"""
rules/10_quantify_genes.smk
Approximate relative gene expression across a cohort, run once per (bed,
sample_type) group -- one gene expression matrix + one boxplot-per-gene PDF
set per quantification method. See docs/rules/10_quantify_genes.md.

Four rules:
  1. _10A_plot_relative_gene_count    -- read-count-based proxy (number of
                                         distinct reads overlapping each
                                         gene's BED region), via
                                         scripts/quantify_gene_expression.py
                                         --metric count.
  2. _10B_plot_relative_gene_coverage -- max-coverage-based proxy (peak
                                         per-base pileup depth in each
                                         gene's BED region), via the same
                                         script with --metric coverage.
  3. _10C1-_10C6 (AMALGAM)             -- isoform discovery/quantification
                                         via AMALGAM (github.com/RNA-ROB/amalgam),
                                         run as its own dedicated 6-rule
                                         sub-pipeline: per-sample StringTie
                                         assembly (_10C1), group-level
                                         GffCompare merge against the
                                         reference annotation (_10C2),
                                         AMALGAM's own transcript filtering
                                         (_10C3) and (optional, currently
                                         unused downstream -- see _10C4's
                                         docstring) ORF annotation (_10C4),
                                         per-sample transcript quantification
                                         against the group's filtered
                                         transcriptome (_10C5), and a final
                                         cohort-wide matrix aggregation
                                         (_10C6). See setup.sh's AMALGAM
                                         section for how the tool itself and
                                         its dedicated conda env get
                                         installed.
  4. _10D_plot_relative_gene_by_assignment -- splice-site-sharing-based
                                         proxy: assigns each alignment to
                                         whichever GTF gene (genome-wide, not
                                         just genes on the BED panel) it
                                         shares the most annotated splice
                                         sites with, via
                                         scripts/quantify_gene_by_assignment.py.
                                         See that script's module docstring
                                         for the assignment rules (spliced
                                         vs. unspliced alignments, ties, and
                                         the genome-wide gene index used to
                                         avoid comparing every read against
                                         every gene in the GTF).

_10A and _10B report two things per gene per sample: a raw value (read
count, or max depth) and a CPTM ("counts per target million") value =
raw_value / (sum of every gene's raw value for that sample, i.e. every gene
on this BED panel) * 1e6. _10D reports the analogous raw assigned-read count
and CPTM, but ALSO writes a raw matrix across every GTF gene that received
>=1 assigned read anywhere in the cohort (assignment isn't restricted to the
BED panel, even though CPTM/boxplots are). CPTM is what the matrices and
boxplots primarily report -- it's the "relative expression" number,
comparable across samples regardless of depth, and normalized against the
panel's own total signal (not overall sequencing depth) since this is
targeted, not whole-transcriptome, sequencing.

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


# ---------------------------------------------------------------------------
# AMALGAM sub-pipeline helpers (_10C1-_10C6). AMALGAM's own tool + assets +
# dedicated conda env live under config["amalgam_dir"] (see setup.sh's
# AMALGAM section) -- referenced as-is throughout, same convention already
# used for config["annovar_dir"] elsewhere in this pipeline.
# ---------------------------------------------------------------------------
def _amalgam_sample_group_id(sample):
    """Which (bed_id, sample_type) group a sample belongs to -- needed
    because StringTie (_10C1) and transcript quantification (_10C5) are
    per-sample steps, but depend on / feed into the group-level merged
    transcriptome (_10C2-_10C4) the same way _10A/_10B/_10D's group-level
    matrices work."""
    s = SAMPLES[sample]
    return _group_id(s["bed"], s["sample_type"])


def _amalgam_group_dir(bed_id, sample_type):
    return _cohort_outdir + "/" + str(bed_id) + "/output/sample_types/" + str(sample_type) + "/output/gene_quantification/by_amalgam"


def _amalgam_stringtie_gtf(bed_id, sample_type, sample):
    return _amalgam_group_dir(bed_id, sample_type) + "/stringtie/" + sample + "_stringtie.gtf"


def _amalgam_quantification_tsv(bed_id, sample_type, sample):
    return _amalgam_group_dir(bed_id, sample_type) + "/quantification/" + sample + "_transcript_quantification.tsv"


rule _10C1_amalgam_stringtie:
    # Step 1 of AMALGAM's own pipeline (see its README): de-novo
    # transcriptome assembly per sample, short-read mode/default settings.
    # Lives under the GROUP's by_amalgam/stringtie/ directory (not the
    # sample's own directory) even though StringTie itself only looks at
    # one BAM at a time -- bed_id/sample_type/sample are all real
    # wildcards here so every step of this sub-pipeline's output lands in
    # one place: by_amalgam/{stringtie,annotation,quantification}/.
    input:
        bam = lambda wc: SAMPLES[wc.sample]["bam"],
    output:
        gtf = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/stringtie/{sample}_stringtie.gtf",
    params:
        amalgam_env = config["amalgam_dir"] + "/conda_env",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample}_amalgam_stringtie.log"
    shell:
        """
        mkdir -p $(dirname {output.gtf}) $(dirname {log})
        (
            export PATH="{params.amalgam_env}/bin:$PATH"
            stringtie -o {output.gtf} {input.bam}
            echo -e "\\nFinished running StringTie on {wildcards.sample}."
        ) 2>&1 | tee {log}
        """


rule _10C2_amalgam_merge_gtfs:
    # Step 2: GffCompare merges every sample's StringTie GTF in this
    # group with the reference annotation into one combined transcript
    # set. The reference annotation MUST be first in the input list (see
    # AMALGAM's README) -- gtf_list.tsv is built in that order below.
    input:
        annotation  = config["annotation"],
        sample_gtfs = lambda wc: [
            _amalgam_stringtie_gtf(wc.bed_id, wc.sample_type, s)
            for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]
        ],
    output:
        combined_gtf = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/merged.combined.gtf",
        tracking      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/merged.tracking",
    params:
        amalgam_env = config["amalgam_dir"] + "/conda_env",
        outprefix   = lambda wc: _amalgam_group_dir(wc.bed_id, wc.sample_type) + "/annotation/merged",
        gtf_list    = lambda wc: _amalgam_group_dir(wc.bed_id, wc.sample_type) + "/annotation/gtf_list.tsv",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 32,
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_merge_gtfs.log"
    shell:
        """
        mkdir -p $(dirname {output.combined_gtf}) $(dirname {log})
        (
            export PATH="{params.amalgam_env}/bin:$PATH"
            echo "{input.annotation}" > {params.gtf_list}
            for f in {input.sample_gtfs}; do echo "$f" >> {params.gtf_list}; done
            gffcompare -i {params.gtf_list} -T -o {params.outprefix}
            echo "Finished merging GTFs."
        ) 2>&1 | tee {log}
        """


rule _10C3_amalgam_build_transcriptome:
    # Step 3: Build_Transcriptome.py identifies high-confidence,
    # full-length transcripts from the GffCompare merge, using AMALGAM's
    # bundled RefTSS/PolyASite reference BED files (config["amalgam_dir"]/
    # assets/, cloned alongside the tool itself -- see setup.sh). Keeps
    # BOTH the uncompressed filtered.gtf (needed as-is by _10C4 below,
    # matching AMALGAM's own README) and the sorted/bgzip/tabix-indexed
    # filtered.gtf.gz (needed by _10C5) as real tracked outputs.
    input:
        combined_gtf = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/merged.combined.gtf",
        tracking      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/merged.tracking",
        annotation    = config["annotation"],
        genome        = config["genome"],
    output:
        filtered_gtf    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/filtered.gtf",
        filtered_gtf_gz = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/filtered.gtf.gz",
        filtered_tbi    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/filtered.gtf.gz.tbi",
    params:
        amalgam_env    = config["amalgam_dir"] + "/conda_env",
        main_env       = workflow.basedir + "/envs/conda_env",  # bgzip/tabix (htslib) -- not part of AMALGAM's own env
        amalgam_dir    = config["amalgam_dir"],
        merge_prefix   = lambda wc: _amalgam_group_dir(wc.bed_id, wc.sample_type) + "/annotation/merged",
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_build_transcriptome.log"
    shell:
        """
        mkdir -p $(dirname {output.filtered_gtf}) $(dirname {log})
        (
            export PATH="{params.amalgam_env}/bin:{params.main_env}/bin:$PATH"
            python -u {params.amalgam_dir}/scripts/Build_Transcriptome.py \\
                -i {params.merge_prefix} \\
                -g {input.annotation} \\
                -f {input.genome} \\
                -x {params.amalgam_dir}/assets/human.refTSS_v4.1.hg38.bed.gz \\
                -y {params.amalgam_dir}/assets/atlas.clusters.2.0.GRCh38.bed.gz \\
                -o {output.filtered_gtf}
            echo "Finished filtering GTF."
            sort -k1,1V -k4,4g -k5,5g {output.filtered_gtf} | bgzip > {output.filtered_gtf_gz}
            tabix -p gff {output.filtered_gtf_gz}
            echo "Finished sorting and indexing GTF."
        ) 2>&1 | tee {log}
        """


rule _10C4_amalgam_annotate_orf:
    # Step 4 (OPTIONAL per AMALGAM's own README): Annotate_ORF.py adds
    # open-reading-frame annotations to the filtered transcriptome.
    # NOTE: this step's output is currently NOT consumed by anything
    # downstream -- _10C5 (transcript quantification) runs against
    # _10C3's filtered.gtf.gz, not this rule's annotated.gtf.gz, matching
    # exactly how the group's own AMALGAM submission script was written
    # (Quantify_Transcripts.py -g pointed at step3's output, not step4's).
    # Kept as a real rule (rather than dropped) since the group's script
    # ran it unconditionally, but flagging here in case that was meant to
    # feed step 5 and didn't due to an oversight in the original script.
    input:
        filtered_gtf = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/filtered.gtf",
        annotation   = config["annotation"],
        genome       = config["genome"],
    output:
        annotated_gtf    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/annotated.gtf",
        annotated_gtf_gz = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/annotated.gtf.gz",
        annotated_tbi    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/annotated.gtf.gz.tbi",
    params:
        amalgam_env = config["amalgam_dir"] + "/conda_env",
        main_env    = workflow.basedir + "/envs/conda_env",
        amalgam_dir = config["amalgam_dir"],
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,  # AMALGAM's README: ~8GB observed
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_annotate_orf.log"
    shell:
        """
        mkdir -p $(dirname {output.annotated_gtf}) $(dirname {log})
        (
            export PATH="{params.amalgam_env}/bin:{params.main_env}/bin:$PATH"
            python -u {params.amalgam_dir}/scripts/Annotate_ORF.py \\
                -i {input.filtered_gtf} \\
                -a {input.annotation} \\
                -f {input.genome} \\
                -o {output.annotated_gtf}
            sort -k1,1V -k4,4g -k5,5g {output.annotated_gtf} | bgzip > {output.annotated_gtf_gz}
            tabix -p gff {output.annotated_gtf_gz}
        ) 2>&1 | tee {log}
        """


rule _10C5_amalgam_quantify_transcripts:
    # Step 5: Quantify_Transcripts.py, per sample, against the GROUP's
    # filtered transcriptome from _10C3 (every sample in a group shares
    # the same transcriptome; only the BAM being quantified differs).
    # bed_id/sample_type are real wildcards here (not resolved from
    # sample via a helper) since output now lives under the group's
    # by_amalgam/quantification/ directory rather than the sample's own.
    input:
        bam    = lambda wc: SAMPLES[wc.sample]["bam"],
        gtf_gz = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/filtered.gtf.gz",
        tbi    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/filtered.gtf.gz.tbi",
    output:
        tsv = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/{sample}_transcript_quantification.tsv",
    params:
        amalgam_env = config["amalgam_dir"] + "/conda_env",
        amalgam_dir = config["amalgam_dir"],
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,  # AMALGAM's README: ~3.5GB observed
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample}_amalgam_quantify_transcripts.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        (
            export PATH="{params.amalgam_env}/bin:$PATH"
            python -u {params.amalgam_dir}/scripts/Quantify_Transcripts.py \\
                -i {input.bam} \\
                -g {input.gtf_gz} \\
                -o {output.tsv}
        ) 2>&1 | tee {log}
        """


rule _10C6_amalgam_aggregate_matrices:
    # Step 6 (not part of AMALGAM itself -- the group's own aggregation
    # step from its submission script): combine every sample's
    # transcript-level quantification in this group into cohort-wide
    # transcript and gene matrices. Extracted into its own script
    # (scripts/aggregate_amalgam_matrices.py) rather than kept as inline
    # Python, matching this repo's convention elsewhere.
    input:
        tsvs = lambda wc: [
            _amalgam_quantification_tsv(wc.bed_id, wc.sample_type, s)
            for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]
        ],
    output:
        transcript_matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_transcript_matrix.tsv",
        gene_matrix        = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_gene_matrix.tsv",
    params:
        samples   = lambda wc: GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)],
        outprefix = lambda wc: _amalgam_group_dir(wc.bed_id, wc.sample_type) + "/quantification/gene_amalgam",
        script    = workflow.basedir + "/scripts/aggregate_amalgam_matrices.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_aggregate_matrices.log"
    shell:
        """
        mkdir -p $(dirname {output.gene_matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.tsvs} \\
            --samples   {params.samples} \\
            --outprefix {params.outprefix} \\
        2>&1 | tee {log}
        """


rule _10D_plot_relative_gene_by_assignment:
    input:
        bams = lambda wc: [SAMPLES[s]["bam"] for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]],
        bed  = lambda wc: bed_path(wc.bed_id),
        gtf  = config["annotation"],
    output:
        mapping_file = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_mapping_file.tsv",
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw.tsv",
        matrix_raw_all_genes = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw_all_genes.tsv",
    params:
        names     = lambda wc: _quoted(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]),
        bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]]),
        outprefix = lambda wc: (_cohort_outdir + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_assignment/gene_assignment"),
        title     = lambda wc: config.get("gene_quant_title") or (str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (splice-site assignment)"),
        script    = workflow.basedir + "/scripts/quantify_gene_by_assignment.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "plot_relative_gene_by_assignment", config["threads"])
    resources:
        # Each worker parses the full GTF once (genome-wide gene index, not
        # just the BED panel), on top of the per-sample BAM scan -- higher
        # floor than _10A/_10B to account for that shared parsing cost.
        mem_mb  = lambda wc, attempt: attempt * 2048 * max(8, len(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_assignment_quantification.log"
    shell:
        """
        mkdir -p $(dirname {output.mapping_file}) $(dirname {log})

        paste <(printf '%s\\n' {params.names}) \\
              <(printf '%s\\n' {params.bams})  \\
        > {output.mapping_file}

        python -u {params.script} \\
            --mapping-file {output.mapping_file} \\
            --bed          {input.bed} \\
            --gtf          {input.gtf} \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """
