
_cohort_outdir = config["output_dir"] + "/{cohort_id}"

def _group_gene_bam_mapping_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/phased_reads/' + str(s) + '_gene_bam_mapping_file.tsv') for s in GROUPS[group_id]]


_cja_manifest_path = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{bed_id}_{sample_type}_gene_manifest.tsv"


rule _8A_cohort_junction_analysis:
    input:
        mapping_files = lambda wc: _group_gene_bam_mapping_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        bed           = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
    output:
        cohort_mapping = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{bed_id}_{sample_type}_gene_bam_mapping_file.tsv",
        manifest = _cja_manifest_path,
        note = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{bed_id}_{sample_type}_note.txt",
    params:
        raw_outdir  = lambda wc: (str(group_outdir(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))) + '/cohort_junction_analysis/' + str(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)) + '_raw'),
        genome      = config["genome"],
        cov_thr     = config["sample_coverage_threshold"],
        phasing_thr = config["cohort_jxn_phasing_threshold"],
        min_reads   = config["cohort_jxn_min_reads"],
        min_samples = config["cohort_jxn_min_samples"],
        script      = workflow.basedir + "/scripts/cohort_junction_analysis.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "cohort_junction_analysis", config["threads"])
    resources:
        mem_mb     = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime    = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/cohort_junction_analysis.log"
    shell:


rule _8B_identify_cohort_junction_outliers:
    input:
        manifest = _cja_manifest_path,
        bed      = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
    output:
        outliers          = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{cohort_id}_{bed_id}_{sample_type}_{thr_label}/{cohort_id}_{bed_id}_{sample_type}_outliers.tsv",
        outliers_filtered = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{cohort_id}_{bed_id}_{sample_type}_{thr_label}/{cohort_id}_{bed_id}_{sample_type}_outliers_filtered.tsv",
        outliers_alias    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{cohort_id}_{bed_id}_{sample_type}_{thr_label}/{cohort_id}_{bed_id}_{sample_type}_outliers_alias.tsv",
    params:
        outprefix  = lambda wc: (str(group_outdir(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))) + '/cohort_junction_analysis/' + str(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))),
        has_ipa    = "--has-ipa" if config["genome"] else "",
        thr_flag   = lambda wc: _cja_thr_flag(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        gtf        = config["annotation"],
        cov_thr    = config["sample_coverage_threshold"],
        n_thr      = lambda wc: _cja_n_threshold(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/identify_cohort_junction_outliers.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "identify_cohort_junction_outliers", config["threads"])
    resources:
        mem_mb     = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime    = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/identify_cohort_junction_outliers_{thr_label}.log"
    shell:
