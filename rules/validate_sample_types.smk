"""
rules/validate_sample_types.smk
Runs once per BED panel, across every sample_type that shares that panel:
1. build_group_junction_matrix  – combine per-sample junction count matrices into
                                   one matrix per (bed, sample_type) group
2. validate_sample_types        – PSI-based tissue-identity check of every group's
                                   samples against GTEx reference tissues

Output is nested as {merged_outdir}/{bed_id}/{sample_type}/... (see merge_hits.smk),
so build_group_junction_matrix also uses bed_id/sample_type as two separate
wildcards and reconstructs group_id via _group_id_from_ids() as needed.

Replaces the "Validate sample types" section of the old bash wrapper.
(final_merge_hits_by_bed, formerly the third rule in this file, has moved to
rules/merge_hits.smk and been renamed final_merge.)
"""

def _group_junction_matrix_inputs(group_id):
    return [f"{SAMPLES[s]['outdir']}/junction_analysis/junction_counts/{s}_junction_count_matrix.tsv"
            for s in GROUPS[group_id]]


rule build_group_junction_matrix:
    input:
        matrices = lambda wc: _group_junction_matrix_inputs(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        matrix = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/junction_analysis/junction_count_matrix.tsv",
    params:
        samples = lambda wc: GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)],
        script  = workflow.basedir + "/scripts/build_group_junction_matrix.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "build_group_junction_matrix", 1)
    resources:
        # Scales with group size since the matrix this rule builds has one
        # column per sample in the group -- large groups (e.g. many
        # fibroblast samples) need proportionally more memory. 1GB/sample
        # with an 8GB floor for small groups by default; override per-group
        # via `groups: <group_id>: build_group_junction_matrix_mem_gb` in
        # the run config if this still under/over-shoots in practice.
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "build_group_junction_matrix",
            max(8, len(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)])))),
        runtime = config["time"],
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/group_junction_matrix.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --infiles      {input.matrices} \\
            --sample-names {params.samples} \\
            --outfile      {output.matrix} \\
        2>&1 | tee {log}
        """


rule validate_sample_types:
    input:
        query_matrices = lambda wc: [
            f"{group_outdir(gid)}/junction_analysis/junction_count_matrix.tsv"
            for gid in BED_GROUPS[wc.bed_id]
        ],
        gtex_matrices = lambda wc: [_gtex_file(t) for t in config["validate_ref_tissues"]],
        bed = lambda wc: bed_path(wc.bed_id),
    output:
        heatmap = f"{config['merged_outdir']}/{{bed_id}}/validate_sample_types/{{bed_id}}_distance_heatmap.pdf",
        pca     = f"{config['merged_outdir']}/{{bed_id}}/validate_sample_types/{{bed_id}}_PCA.pdf",
    params:
        outprefix    = lambda wc: f"{bed_outdir(wc.bed_id)}/validate_sample_types/{wc.bed_id}",
        ref_names    = lambda wc: _quoted(config["validate_ref_tissues"]),
        ref_colors   = lambda wc: _quoted(config["validate_ref_colors"]),
        query_names  = lambda wc: _quoted([GROUP_SAMPLE_TYPE[gid] for gid in BED_GROUPS[wc.bed_id]]),
        query_colors = lambda wc: _quoted([sample_type_color(GROUP_SAMPLE_TYPE[gid]) for gid in BED_GROUPS[wc.bed_id]]),
        script       = workflow.basedir + "/scripts/validate_sample_type.py",
    threads: lambda wc: _group_threads(wc.bed_id, "validate_sample_types", 1)
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(wc.bed_id, "validate_sample_types", 256)),
        runtime = 1440,   # 1 day, matches the old bash wrapper's --time=1-00:00:00
    log:
        f"{config['merged_outdir']}/{{bed_id}}/logs/{{bed_id}}_validate_sample_types.log"
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
