from math import ceil

_cohort_outdir  = config["output_dir"] + "/{cohort_id}"
_JXN_FILTER_SETS = ["default", "stringent"]

_variant_tsv = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_variant_calling/{filter_set}/all_candidate_variants.tsv"
_ase_tsv     = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_ase_analysis/{filter_set}/outlier_ase.tsv"

def _variant_sample_fraction(filter_set):
    if filter_set == "stringent":
        return config.get("merge_variant_sample_fraction_stringent", 0.01)
    return config["merge_variant_sample_fraction"]

def _ase_sample_fraction(filter_set):
    if filter_set == "stringent":
        return config.get("merge_ase_sample_fraction_stringent", 0.01)
    return config["merge_ase_sample_fraction"]

def _jxn_sample_fraction(filter_set):
    if filter_set == "stringent":
        return config.get("merge_jxn_sample_fraction_stringent", 0.01)
    return config["merge_jxn_sample_fraction"]

_junction_final = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/{filter_set}/gtex_{tissue}/outlier_junctions_gtex_{tissue}_final.tsv"
_junction_significant = (_cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/{filter_set}/gtex_{tissue}/outlier_junctions_gtex_{tissue}_"
                         + str(int(config["merge_jxn_coverage_threshold"])) + "jxncov_"
                         + str(float(config["merge_jxn_padj_threshold"])) + "padj_"
                         + str(float(config["merge_delta_psi_threshold"])) + "deltaPSI.tsv")
_all_hits_tsv   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/all_hits.tsv"

def _group_variant_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/variant_calling/compiled_variants/' + str(s) + '_filtered_variants.tsv') for s in GROUPS[group_id]]

def _group_ase_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/ase_analysis/' + str(s) + '_binomial_ase_results.tsv') for s in GROUPS[group_id]]

def _group_tissue_samples(group_id, tissue):
    return [s for s in GROUPS[group_id] if tissue in sample_tissues(s)]

def _group_tissue_junction_files(group_id, tissue):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/gtex_' + str(tissue) + '/' + str(s) + '_gtex_' + str(tissue) + '_all_junctions.tsv')
            for s in _group_tissue_samples(group_id, tissue)]

def _group_junction_outprefix(group_id, tissue, filter_set):
    return (str(group_outdir(group_id)) + '/merged_junction_analysis/' + str(filter_set) + '/gtex_' + str(tissue) + '/outlier_junctions_gtex_' + str(tissue))

def _group_junction_final_path(group_id, tissue, filter_set):
    return (str(_group_junction_outprefix(group_id, tissue, filter_set)) + '_final.tsv')

def _group_junction_significant_path(group_id, tissue, filter_set):
    return (str(_group_junction_outprefix(group_id, tissue, filter_set)) + '_'
            + str(int(config['merge_jxn_coverage_threshold'])) + 'jxncov_'
            + str(float(config['merge_jxn_padj_threshold'])) + 'padj_'
            + str(float(config['merge_delta_psi_threshold'])) + 'deltaPSI.tsv')

def _group_junction_source_glob(group_id, tissue, filter_set):
    n_tissue_samples = len(_group_tissue_samples(group_id, tissue))
    n = ceil(n_tissue_samples * _jxn_sample_fraction(filter_set))
    base = ((str(_group_junction_outprefix(group_id, tissue, filter_set)) + '_')
            + (str(config['merge_jxn_coverage_threshold']) + 'jxncov_')
            + (str(config['merge_jxn_padj_threshold']) + 'padj_')
            + (str(config['merge_delta_psi_threshold']) + 'deltaPSI_event'))
    return (str(base) + '*_' + str(n) + 'samples.tsv')


def _sample_qc_file(sample, name):
    return str(SAMPLES[sample]['outdir']) + '/qc/' + str(sample) + '_' + str(name) + '.tsv'

def _bed_qc_files(cohort_id, bed_id, name):
    return [_sample_qc_file(s, name) for s in bed_samples(cohort_id, bed_id)]


rule _10A_merge_on_target_rates:
    input:
        infiles = lambda wc: _bed_qc_files(wc.cohort_id, wc.bed_id, "on_target"),
    output:
        tsv          = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates.tsv",
        mapping_pdf  = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_mapping_rates.pdf",
        ontarget_pdf = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates.pdf",
        ontarget_pdf_alias = _cohort_outdir + "/{bed_id}/output/cohort_qc/on_target_rates/{bed_id}_on_target_rates_alias.pdf",
    params:
        groups     = lambda wc: _quoted([SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        outprefix  = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/on_target_rates/' + str(wc.bed_id)),
        title      = lambda wc: config.get("cohort_qc_title", (str(wc.cohort_id) + ' ' + str(wc.bed_id) + ' on-target rates')),
        alias_args = lambda wc: _quoted(alias_map_args(bed_samples(wc.cohort_id, wc.bed_id))),
        colors     = lambda wc: _quoted([sample_type_color(SAMPLES[s]["sample_type"]) for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        script     = workflow.basedir + "/scripts/_10A_merge_on_target_rates.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_10A_merge_on_target_rates.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --groups    {params.groups} \\
            --colors    {params.colors} \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _10B_merge_read_attributes:
    input:
        infiles = lambda wc: _bed_qc_files(wc.cohort_id, wc.bed_id, "read_attributes"),
    output:
        tsv           = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes.tsv",
        length_boxplot = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_lengths.pdf",
        length_boxplot_alias = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_lengths_alias.pdf",
    params:
        outprefix  = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/read_attributes/' + str(wc.bed_id)),
        title      = lambda wc: config.get("cohort_qc_title", (str(wc.cohort_id) + ' ' + str(wc.bed_id) + ' read attributes')),
        alias_args = lambda wc: _quoted(alias_map_args(bed_samples(wc.cohort_id, wc.bed_id))),
        script     = workflow.basedir + "/scripts/_10B_merge_read_attributes.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_10B_merge_read_attributes.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


def _group_junction_matrix_inputs(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/junction_counts/' + str(s) + '_junction_count_matrix.tsv')
            for s in GROUPS[group_id]]


rule _10C1_build_group_junction_matrix:
    input:
        matrices = lambda wc: _group_junction_matrix_inputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        bed      = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
    output:
        matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/junction_count_matrix.tsv",
    params:
        samples = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        script  = workflow.basedir + "/scripts/_10C1_build_group_junction_matrix.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "build_group_junction_matrix", config["threads"])
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_10C1_build_group_junction_matrix.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --infiles      {input.matrices} \\
            --sample-names {params.samples} \\
            --outfile      {output.matrix} \\
            --bed          {input.bed} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """


rule _10C2_validate_sample_types:
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
        pca_coords = _cohort_outdir + "/{bed_id}/output/cohort_qc/validate_sample_types/{bed_id}_PCA_coords.tsv",
        heatmap_alias = _cohort_outdir + "/{bed_id}/output/cohort_qc/validate_sample_types/{bed_id}_distance_heatmap_alias.pdf",
    params:
        outprefix    = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/validate_sample_types/' + str(wc.bed_id)),
        ref_names    = lambda wc: _quoted(config["validate_ref_tissues"]),
        ref_colors   = lambda wc: _quoted([validate_ref_color(t) for t in config["validate_ref_tissues"]]),
        query_names  = lambda wc: _quoted([GROUP_SAMPLE_TYPE[gid] for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]]),
        query_colors = lambda wc: _quoted([sample_type_color(GROUP_SAMPLE_TYPE[gid]) for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]]),
        alias_args   = lambda wc: _quoted(alias_map_args(bed_samples(wc.cohort_id, wc.bed_id))),
        script       = workflow.basedir + "/scripts/_10C2_validate_sample_types.py",
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "validate_sample_types", 1)
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(32, len(bed_samples(wc.cohort_id, wc.bed_id))),
        runtime = 1440,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_10C2_validate_sample_types.log"
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
            --alias-map    {params.alias_args} \\
        2>&1 | tee {log}
        """


_AVG_COV_DIR = _cohort_outdir + "/{bed_id}/output/cohort_qc/avg_coverage/{bed_id}"

rule _10D_merge_avg_coverage:
    input:
        infiles         = lambda wc: _bed_qc_files(wc.cohort_id, wc.bed_id, "avg_coverage"),
        on_target_files = lambda wc: _bed_qc_files(wc.cohort_id, wc.bed_id, "on_target"),
    output:
        avg_matrix          = _AVG_COV_DIR + "_avg_coverage_matrix.tsv",
        avg_matrix_alias    = _AVG_COV_DIR + "_avg_coverage_matrix_alias.tsv",
        avg_heatmap         = _AVG_COV_DIR + "_avg_coverage_heatmap.pdf",
        avg_heatmap_alias   = _AVG_COV_DIR + "_avg_coverage_heatmap_alias.pdf",
        median_matrix        = _AVG_COV_DIR + "_median_coverage_matrix.tsv",
        median_matrix_alias  = _AVG_COV_DIR + "_median_coverage_matrix_alias.tsv",
        median_heatmap       = _AVG_COV_DIR + "_median_coverage_heatmap.pdf",
        median_heatmap_alias = _AVG_COV_DIR + "_median_coverage_heatmap_alias.pdf",
        breadth_matrix        = _AVG_COV_DIR + "_breadth_coverage_matrix.tsv",
        breadth_matrix_alias  = _AVG_COV_DIR + "_breadth_coverage_matrix_alias.tsv",
        breadth_heatmap       = _AVG_COV_DIR + "_breadth_coverage_heatmap.pdf",
        breadth_heatmap_alias = _AVG_COV_DIR + "_breadth_coverage_heatmap_alias.pdf",
    params:
        samples      = lambda wc: _quoted(bed_samples(wc.cohort_id, wc.bed_id)),
        sample_types = lambda wc: _quoted([SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        colors       = lambda wc: _quoted([sample_type_color(SAMPLES[s]["sample_type"]) for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        outprefix    = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/avg_coverage/' + str(wc.bed_id)),
        title        = lambda wc: config.get("cohort_qc_title", (str(wc.cohort_id) + ' ' + str(wc.bed_id) + ' canonical transcript coverage')),
        center       = config.get("coverage_heatmap_center", 50),
        alias_args   = lambda wc: _quoted(alias_map_args(bed_samples(wc.cohort_id, wc.bed_id))),
        script       = workflow.basedir + "/scripts/_10D_merge_avg_coverage.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_10D_merge_avg_coverage.log"
    shell:
        """
        mkdir -p $(dirname {output.avg_matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles         {input.infiles} \\
            --on-target-files {input.on_target_files} \\
            --samples         {params.samples} \\
            --sample-types    {params.sample_types} \\
            --colors          {params.colors} \\
            --outprefix       {params.outprefix} \\
            --title           {params.title:q} \\
            --center          {params.center} \\
            --alias-map       {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _10E_merge_full_length_ratio:
    input:
        infiles = lambda wc: _bed_qc_files(wc.cohort_id, wc.bed_id, "full_length_ratio"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_matrix.tsv",
        read_counts  = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_read_counts_matrix.tsv",
        heatmap      = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_heatmap.pdf",
        matrix_alias = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_matrix_alias.tsv",
    params:
        sample_types = lambda wc: _quoted([SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        colors     = lambda wc: _quoted([sample_type_color(SAMPLES[s]["sample_type"]) for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        outprefix  = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/full_length_ratio/' + str(wc.bed_id) + '_full_length_ratio'),
        title      = lambda wc: config.get("cohort_qc_title", (str(wc.cohort_id) + ' ' + str(wc.bed_id) + ' full-length ratio')),
        min_reads  = config.get("full_length_ratio_min_reads", 100),
        alias_args = lambda wc: _quoted(alias_map_args(bed_samples(wc.cohort_id, wc.bed_id))),
        script     = workflow.basedir + "/scripts/_10E_merge_full_length_ratio.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_10E_merge_full_length_ratio.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles      {input.infiles} \\
            --sample-types {params.sample_types} \\
            --colors       {params.colors} \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --min-reads    {params.min_reads} \\
            --alias-map    {params.alias_args} \\
        2>&1 | tee {log}
        """


def _sample_quant_file(sample, name):
    return str(SAMPLES[sample]['outdir']) + '/gene_quantification/' + str(sample) + '_' + str(name) + '.tsv'

def _group_quant_files(group_id, name):
    return [_sample_quant_file(s, name) for s in GROUPS[group_id]]


rule _10F_merge_gene_count:
    input:
        infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_count"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_cptm.tsv",
        matrix_motr  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_motr.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_raw.tsv",
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_count/gene_count"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (read count)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/_10_merge_gene_expression.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_10F_merge_gene_count.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --metric    count \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _10G_merge_gene_coverage:
    input:
        infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_coverage"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_cptm.tsv",
        matrix_motr  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_motr.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_raw.tsv",
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_coverage/gene_coverage"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (max coverage)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/_10_merge_gene_expression.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_10G_merge_gene_coverage.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --metric    coverage \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _10H_merge_gene_expression_by_assignment:
    input:
        infiles       = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_assignment"),
        stats_infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "read_outcomes"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_cptm.tsv",
        matrix_motr  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_motr.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw.tsv",
        matrix_raw_all_genes = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw_all_genes.tsv",
        assignment_stats     = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_read_outcomes.tsv",
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_assignment/gene_assignment"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (splice-site assignment)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/_10_merge_gene_expression_by_assignment.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_10H_merge_gene_expression_by_assignment.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles       {input.infiles} \\
            --stats-infiles {input.stats_infiles} \\
            --outprefix     {params.outprefix} \\
            --title         {params.title:q} \\
            --alias-map     {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _10I_amalgam_normalize_matrix:
    input:
        gene_matrix_raw = rules._9E_amalgam_aggregate_matrices.output.gene_matrix_raw,
    output:
        matrix_cptm = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_cptm.tsv",
        matrix_motr = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_motr.tsv",
        matrix_cptm_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/gene_amalgam",
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (AMALGAM)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/_10I_amalgam_normalize_matrix.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_10I_amalgam_normalize_matrix.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix_cptm}) $(dirname {log})
        python -u {params.script} \\
            --raw-matrix  {input.gene_matrix_raw} \\
            --outprefix   {params.outprefix} \\
            --title       {params.title:q} \\
            --alias-map   {params.alias_args} \\
        2>&1 | tee {log}
        """


_GQ_METHODS = {
    "count":      ("by_count",      "gene_count",      "read count"),
    "coverage":   ("by_coverage",   "gene_coverage",   "max coverage"),
    "assignment": ("by_assignment", "gene_assignment", "assigned reads"),
    "amalgam":    ("by_amalgam",    "gene_amalgam",    "AMALGAM"),
}

def _gq_matrix(cohort_id, bed_id, sample_type, method, norm):
    subdir, prefix, _ = _GQ_METHODS[method]
    return (config["output_dir"] + "/" + str(cohort_id) + "/" + str(bed_id) + "/output/sample_types/"
            + str(sample_type) + "/output/gene_quantification/" + subdir + "/" + prefix + "_matrix_" + norm + ".tsv")

def _gq_types(cohort_id, bed_id):
    return [GROUP_SAMPLE_TYPE[gid] for gid in BED_GROUPS[(cohort_id, bed_id)]]


rule _10J_gene_quantification_boxplots:
    input:
        cptm = lambda wc: [_gq_matrix(wc.cohort_id, wc.bed_id, st, wc.method, "cptm") for st in _gq_types(wc.cohort_id, wc.bed_id)],
        motr = lambda wc: [_gq_matrix(wc.cohort_id, wc.bed_id, st, wc.method, "motr") for st in _gq_types(wc.cohort_id, wc.bed_id)],
    output:
        outdir = directory(_cohort_outdir + "/{bed_id}/output/cohort_qc/gene_quantification/by_{method}"),
    wildcard_constraints:
        method = "count|coverage|assignment|amalgam",
    params:
        cptm_args = lambda wc: _quoted([st + ":" + _gq_matrix(wc.cohort_id, wc.bed_id, st, wc.method, "cptm") for st in _gq_types(wc.cohort_id, wc.bed_id)]),
        motr_args = lambda wc: _quoted([st + ":" + _gq_matrix(wc.cohort_id, wc.bed_id, st, wc.method, "motr") for st in _gq_types(wc.cohort_id, wc.bed_id)]),
        colors    = lambda wc: _quoted([st + ":" + sample_type_color(st) for st in _gq_types(wc.cohort_id, wc.bed_id)]),
        label     = lambda wc: _GQ_METHODS[wc.method][2],
        script    = workflow.basedir + "/scripts/_10J_gene_quantification_boxplots.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 4,
        runtime = 120,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_{method}_10J_gene_quantification_boxplots.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --cptm         {params.cptm_args} \\
            --motr         {params.motr_args} \\
            --colors       {params.colors} \\
            --metric-label {params.label:q} \\
            --outdir       {output.outdir} \\
        2>&1 | tee {log}
        """


rule _10K_merge_group_variants:
    input:
        variant_files = lambda wc: _group_variant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        tsv = _variant_tsv,
    params:
        group_id     = lambda wc: _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type),
        n            = lambda wc: ceil(len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])
                                        * _variant_sample_fraction(wc.filter_set)),
        n_rna        = lambda wc: ceil(len([s for s in GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]
                                            if not SAMPLES[s].get("wgs_vcf")])
                                        * _variant_sample_fraction(wc.filter_set)),
        outprefix    = lambda wc, output: output.tsv[:-len(".tsv")],
        num_callers_snv   = config["merge_num_callers_threshold_snv"],
        num_callers_indel = config["merge_num_callers_threshold_indel"],
        min_dp_snv        = config["merge_min_dp_snv"],
        min_dp_indel      = config["merge_min_dp_indel"],
        script       = workflow.basedir + "/scripts/_10K_merge_group_variants.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "merge_group_variants", 1)
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_{filter_set}_10K_merge_group_variants.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles {input.variant_files} \\
            --outprefix {params.outprefix} \\
            --num-callers-threshold-SNV {params.num_callers_snv} \\
            --num-callers-threshold-indel {params.num_callers_indel} \\
            --min-DP-SNV {params.min_dp_snv} \\
            --min-DP-indel {params.min_dp_indel} \\
            --sample-number-threshold {params.n} \\
            --sample-number-threshold-RNA {params.n_rna} \\
            --plot \\
            --plot-variant-type SNV indel \\
            --title "{params.group_id} ({wildcards.filter_set}) Variant Counts" \\
        2>&1 | tee {log}
        """


rule _10L_merge_group_ase:
    input:
        ase_files = lambda wc: _group_ase_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        tsv = _ase_tsv,
    params:
        group_id     = lambda wc: _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type),
        n            = lambda wc: ceil(len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])
                                        * _ase_sample_fraction(wc.filter_set)),
        outprefix    = lambda wc, output: output.tsv[:-len(".tsv")],
        min_hap_ratio       = config["merge_min_haplotype_ratio"],
        minor_hap_freq_thr  = config["merge_minor_haplotype_frequency_threshold"],
        ase_padj_thr        = config["merge_ase_padj_threshold"],
        script       = workflow.basedir + "/scripts/_10L_merge_group_ase.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "merge_group_ase", 1)
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_{filter_set}_10L_merge_group_ase.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles {input.ase_files} \\
            --outprefix {params.outprefix} \\
            --min-haplotype-ratio {params.min_hap_ratio} \\
            --minor-haplotype-frequency-threshold {params.minor_hap_freq_thr} \\
            --padj-threshold {params.ase_padj_thr} \\
            --plot \\
            --title "{params.group_id} ({wildcards.filter_set}): Number of Genes with Allele-specific Expression by Sample" \\
            --sample-number-threshold {params.n} \\
        2>&1 | tee {log}
        """


rule _10M_merge_group_junctions:
    input:
        junction_files = lambda wc: _group_tissue_junction_files(
            _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue),
        bed            = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
        gtf            = config["annotation"],
    output:
        tsv         = _junction_final,
        significant = _junction_significant,
    params:
        group_id  = lambda wc: _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type),
        n         = lambda wc: ceil(len(_group_tissue_samples(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue))
                                     * _jxn_sample_fraction(wc.filter_set)),
        outprefix = lambda wc: _group_junction_outprefix(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue, wc.filter_set),
        source_glob = lambda wc: _group_junction_source_glob(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue, wc.filter_set),
        jxn_cov_thr   = config["merge_jxn_coverage_threshold"],
        jxn_padj_thr  = config["merge_jxn_padj_threshold"],
        delta_psi_thr = config["merge_delta_psi_threshold"],
        script    = workflow.basedir + "/scripts/_10M_merge_group_junctions.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "merge_group_junctions", 1)
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(_group_tissue_samples(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue)) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_{filter_set}_{tissue}_10M_merge_group_junctions.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles {input.junction_files} \\
            --outprefix {params.outprefix} \\
            --jxn-coverage-threshold {params.jxn_cov_thr} \\
            --padj-threshold {params.jxn_padj_thr} \\
            --delta-PSI-threshold {params.delta_psi_thr} \\
            --bed {input.bed} \\
            --gtf {input.gtf} \\
            --event-types exon_skipping_approx exon_inclusion_approx alt_5ss_approx alt_3ss_approx complex_approx \\
            --sample-number-threshold {params.n} \\
            --plot \\
            --title "{params.group_id} ({wildcards.filter_set}): Number of Genes with Outlier Junctions by Sample" \\
        2>&1 | tee {log}
        SRC=$(ls {params.source_glob} 2>/dev/null | head -1)
        if [ -z "$SRC" ]; then
            echo "WARNING: no output file matched glob {params.source_glob}" | tee -a {log} >&2
            exit 1
        fi
        cp "$SRC" {output.tsv}
        echo "Copied $SRC -> {output.tsv}" >> {log}
        """


def _group_all_hits_path(group_id, filter_set):
    return str(group_outdir(group_id)) + "/merged_hits/" + str(filter_set) + "/all_hits.tsv"


def _group_avg_coverage_files(group_id):
    return [_sample_qc_file(s, "avg_coverage") for s in GROUPS[group_id]]

def _group_gtex_all_junction_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/gtex_' + str(t) + '/' + str(s)
             + '_gtex_' + str(t) + '_all_junctions.tsv')
            for s in GROUPS[group_id] for t in sample_tissues(s)]

def _group_cohort_results_dir(group_id):
    return (str(group_outdir(group_id)) + '/cohort_junction_analysis/' + str(group_id)
            + '_results_' + _cja_method_for_group(group_id))

def _group_cohort_scoring_sentinel(cohort_id, bed_id, sample_type):
    return _cja_scoring_sentinel.format(cohort_id=cohort_id, bed_id=bed_id, sample_type=sample_type)

rule _10N1_merge_group_hits:
    input:
        bed            = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
        gtf            = config["annotation"],
        avg_coverage_files = lambda wc: (
            _group_avg_coverage_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ) if flag("qc") else [],
        gtex_all_junction_files = lambda wc: (
            _group_gtex_all_junction_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ) if flag("junction_analysis") else [],
        cohort_scoring_sentinel = lambda wc: (
            _group_cohort_scoring_sentinel(wc.cohort_id, wc.bed_id, wc.sample_type)
        ) if cja_enabled(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)) else [],
        variant_tsv    = _variant_tsv,
        ase_tsv        = _ase_tsv,
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), t, wc.filter_set)
            for t in group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ],
        cohort_junction_tsv = lambda wc: (
            _cja_outliers_filtered_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ) if cja_enabled(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)) else [],
        gene_expression_matrix = lambda wc: (
            config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
            + "/output/gene_quantification/by_assignment/gene_assignment_matrix_cptm.tsv"
        ) if config.get("gene_quantification") else [],
        gene_expression_matrix_motr = lambda wc: (
            config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
            + "/output/gene_quantification/by_assignment/gene_assignment_matrix_motr.tsv"
        ) if config.get("gene_quantification") else [],
        amelie_files = lambda wc: [amelie_tsv(s) for s in GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]
                                   if sample_hpo(s)],
        isoform_outlier_files = lambda wc: (
            ctx_outputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ) if ctx_enabled(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)) else [],
    output:
        all_hits = _all_hits_tsv,
    params:
        samples      = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        cohort_junction_flag = lambda wc, input: (
            "--cohort-junction-tsv " + str(input.cohort_junction_tsv)
        ) if input.cohort_junction_tsv else "",
        gene_expression_flag = lambda wc, input: (
            "--gene-expression-matrix " + str(input.gene_expression_matrix)
        ) if config.get("gene_quantification") else "",
        gene_expression_motr_flag = lambda wc, input: (
            "--gene-expression-matrix-motr " + str(input.gene_expression_matrix_motr)
        ) if config.get("gene_quantification") else "",
        avg_coverage_flag = lambda wc, input: (
            "--avg-coverage-files " + " ".join(input.avg_coverage_files)
        ) if flag("qc") else "",
        gtex_all_junctions_flag = lambda wc, input: (
            "--gtex-all-junction-files " + " ".join(input.gtex_all_junction_files)
        ) if flag("junction_analysis") else "",
        cohort_results_flag = lambda wc: (
            "--cohort-results-dir " + _group_cohort_results_dir(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ) if cja_enabled(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)) else "",
        amelie_flag  = lambda wc, input: ("--amelie-files " + " ".join(input.amelie_files)) if input.amelie_files else "",
        isoform_flag = lambda wc, input: (
            "--isoform-outlier-files " + " ".join(input.isoform_outlier_files)) if input.isoform_outlier_files else "",
        script       = workflow.basedir + "/scripts/_10N1_merge_group_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4 * 1024,
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_{filter_set}_10N1_merge_group_hits.log"
    shell:
        """
        mkdir -p $(dirname {output.all_hits}) $(dirname {log})
        python -u {params.script} \\
            --outfile     {output.all_hits} \\
            --variant-tsv {input.variant_tsv} \\
            --ase-tsv     {input.ase_tsv} \\
            --tissues     {params.tissues} \\
            --junction-files {input.junction_files} \\
            --samples     {params.samples} \\
            {params.omim_flag} \\
            {params.cohort_junction_flag} \\
            {params.gene_expression_flag} \\
            {params.gene_expression_motr_flag} \\
            --bed         {input.bed} \\
            --gtf         {input.gtf} \\
            {params.avg_coverage_flag} \\
            {params.gtex_all_junctions_flag} \\
            {params.cohort_results_flag} \\
            {params.amelie_flag} \\
            {params.isoform_flag} \\
        2>&1 | tee {log}
        """


rule _10N2_plot_group_hits:
    input:
        all_hits = lambda wc: _group_all_hits_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.filter_set),
    output:
        pathogenic_bar  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/genes_with_pathogenic_variant_barplot.pdf",
        ase_bar         = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/genes_with_ASE_barplot.pdf",
        junction_bar    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/genes_with_outlier_junction_barplot.pdf",
        dysreg_bar      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/genes_with_RNA_dysregulation_barplot.pdf",
    params:
        outdir = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}",
        script = workflow.basedir + "/scripts/_10N2_plot_group_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 4 * 1024,
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_{filter_set}_10N2_plot_group_hits.log"
    shell:
        """
        python -u {params.script} \\
            --infile {input.all_hits} \\
            --outdir {params.outdir} \\
        2>&1 | tee {log}
        """


def _method_comparison_dir(group_id, filter_set):
    return str(group_outdir(group_id)) + "/merged_hits/" + str(filter_set) + "/method_comparison"

def _method_comparison_gtex_inputs(group_id):
    return [(s, t, str(SAMPLES[s]['outdir']) + '/junction_analysis/gtex_' + str(t) + '/' + str(s)
             + '_gtex_' + str(t) + '_all_junctions.tsv')
            for s in GROUPS[group_id] for t in sample_tissues(s)]

rule _10N3_compare_junction_methods:
    input:
        bed               = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
        gtex_final        = lambda wc: [_group_junction_final_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), t, wc.filter_set)
                                        for t in group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))],
        gtex_significant  = lambda wc: [_group_junction_significant_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), t, wc.filter_set)
                                        for t in group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))],
        gtex_all          = lambda wc: [x[2] for x in _method_comparison_gtex_inputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))],
        cohort_sentinel   = lambda wc: _group_cohort_scoring_sentinel(wc.cohort_id, wc.bed_id, wc.sample_type),
        cohort_outliers   = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)).replace("_outliers_filtered.tsv", "_outliers.tsv"),
        cohort_filtered   = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        isoform_outliers  = lambda wc: (
            ctx_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))["isoform_outliers"]
        ) if ctx_enabled(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)) else [],
        landscape_outliers = lambda wc: (
            ctx_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))["landscape_outliers"]
        ) if ctx_enabled(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)) else [],
    output:
        pdf   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/method_comparison/{sample_type}_method_comparison.pdf",
        genes = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/method_comparison/{sample_type}_gene_hits_by_method.tsv",
        jxns  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/{filter_set}/method_comparison/{sample_type}_junction_status.tsv.gz",
    params:
        group_id   = lambda wc: _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type),
        samples    = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        tissues    = lambda wc: group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        thresholds = lambda wc: [ceil(len(_group_tissue_samples(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), t))
                                      * _jxn_sample_fraction(wc.filter_set))
                                 for t in group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))],
        all_samples = lambda wc: [x[0] for x in _method_comparison_gtex_inputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))],
        all_tissues = lambda wc: [x[1] for x in _method_comparison_gtex_inputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))],
        results_dir = lambda wc: _group_cohort_results_dir(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        iso_flag   = lambda wc, input: ("--isoform-outliers " + str(input.isoform_outliers)
                                        + " --landscape-outliers " + str(input.landscape_outliers)) if input.isoform_outliers else "",
        outdir     = lambda wc: _method_comparison_dir(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.filter_set),
        min_reads  = config["jxn_min_reads"],
        cov_thr    = config["sample_coverage_threshold"],
        n_thr      = lambda wc: _cja_n_threshold(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        script     = workflow.basedir + "/scripts/_10N3_compare_junction_methods.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 4),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_{filter_set}_10N3_compare_junction_methods.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        python -u {params.script} \\
            --samples                   {params.samples} \\
            --bed                       {input.bed} \\
            --tissues                   {params.tissues} \\
            --gtex-final-files          {input.gtex_final} \\
            --gtex-significant-files    {input.gtex_significant} \\
            --gtex-sample-thresholds    {params.thresholds} \\
            --gtex-all-junction-files   {input.gtex_all} \\
            --gtex-all-junction-samples {params.all_samples} \\
            --gtex-all-junction-tissues {params.all_tissues} \\
            --cohort-results-dir        {params.results_dir} \\
            --cohort-outliers           {input.cohort_outliers} \\
            --cohort-outliers-filtered  {input.cohort_filtered} \\
            {params.iso_flag} \\
            --min-jxn-reads             {params.min_reads} \\
            --coverage-threshold        {params.cov_thr} \\
            --n-threshold               {params.n_thr} \\
            --title                     "{params.group_id} ({wildcards.filter_set})" \\
            --outdir                    {params.outdir} \\
            --prefix                    {wildcards.sample_type} \\
        2>&1 | tee {log}
        """


rule _10O1_final_merge:
    input:
        all_hits = lambda wc: [_group_all_hits_path(gid, wc.filter_set) for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]],
    output:
        merged = _cohort_outdir + "/{bed_id}/output/{filter_set}/merged_all_hits.tsv",
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "final_merge", 1)
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_{filter_set}_10O1_final_merge.log"
    shell:
        """
        mkdir -p $(dirname {output.merged}) $(dirname {log})
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.all_hits} > {output.merged} 2> {log}
        echo "Finished final merge to {output.merged}." >> {log}
        """


rule _10O2_plot_hits_upset:
    input:
        all_hits = _cohort_outdir + "/{bed_id}/output/{filter_set}/merged_all_hits.tsv",
    output:
        pdf_density      = _cohort_outdir + "/{bed_id}/output/{filter_set}/hits_upset_density.pdf",
        tsv              = _cohort_outdir + "/{bed_id}/output/{filter_set}/hits_upset_counts.tsv",
        tier_pdf_density = _cohort_outdir + "/{bed_id}/output/{filter_set}/tier_gene_counts_density.pdf",
        tier_tsv         = _cohort_outdir + "/{bed_id}/output/{filter_set}/tier_gene_counts.tsv",
    params:
        samples      = lambda wc: bed_samples(wc.cohort_id, wc.bed_id),
        sample_types = lambda wc: [SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
        outdir       = _cohort_outdir + "/{bed_id}/output/{filter_set}",
        title        = lambda wc: f"{wc.bed_id} ({wc.filter_set}): Candidate Gene Hit Categories by Sample",
        colors       = lambda wc: _quoted([sample_type_color(SAMPLES[s]["sample_type"]) for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        script       = workflow.basedir + "/scripts/_10O2_plot_hits_upset.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 4,
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_{filter_set}_10O2_plot_hits_upset.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        python -u {params.script} \\
            --infile      {input.all_hits} \\
            --samples     {params.samples} \\
            --sample-types {params.sample_types} \\
            --colors      {params.colors} \\
            --outdir      {params.outdir} \\
            --title       "{params.title}" \\
        2>&1 | tee {log}
        """


def _bed_cohort_manifest_args(cohort_id, bed_id):
    pairs = []
    for gid in BED_GROUPS[(cohort_id, bed_id)]:
        st = GROUP_SAMPLE_TYPE[gid]
        manifest = (config["output_dir"] + "/" + str(cohort_id) + "/" + str(bed_id)
                    + "/output/sample_types/" + str(st)
                    + "/output/cohort_junction_analysis/" + str(bed_id) + "_" + str(st) + "_gene_manifest.tsv")
        pairs.append(str(st) + ':' + manifest)
    return pairs

def _bed_cja_manifests(cohort_id, bed_id):
    files = []
    for gid in BED_GROUPS[(cohort_id, bed_id)]:
        st = GROUP_SAMPLE_TYPE[gid]
        files.append(config["output_dir"] + "/" + str(cohort_id) + "/" + str(bed_id)
                     + "/output/sample_types/" + str(st)
                     + "/output/cohort_junction_analysis/" + str(bed_id) + "_" + str(st) + "_gene_manifest.tsv")
    return files


def _bed_expression_matrices(wc):
    return [(GROUP_SAMPLE_TYPE[gid], str(bed_outdir(wc.cohort_id, wc.bed_id)) + "/sample_types/" + str(GROUP_SAMPLE_TYPE[gid])
             + "/output/gene_quantification/by_assignment/gene_assignment_matrix_motr.tsv")
            for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]]


def _cohort_qc_file(wc, rel):
    return str(bed_outdir(wc.cohort_id, wc.bed_id)) + "/cohort_qc/" + rel.format(b=wc.bed_id)


def _bed_wgs_samples(cohort_id, bed_id):
    return [s for s in bed_samples(cohort_id, bed_id) if SAMPLES[s].get("wgs_vcf")]


def _wgs_concordance_summary(wc):
    if not _bed_wgs_samples(wc.cohort_id, wc.bed_id):
        return []
    return _cohort_qc_file(wc, "wgs_concordance/{b}_wgs_concordance_summary.tsv")


def _sample_gene_bam_mapping(sample):
    return str(SAMPLES[sample]['outdir']) + '/phased_reads/' + str(sample) + '_gene_bam_mapping_file.tsv'


rule _10P_check_wgs_concordance:
    input:
        nanots_vcfs = lambda wc: [str(SAMPLES[s]['outdir']) + '/variant_calling/nanoTS/' + str(s) + '_nanoTS_norm.vcf.gz'
                                  for s in bed_samples(wc.cohort_id, wc.bed_id)],
        wgs_vcfs    = lambda wc: [str(SAMPLES[s]['outdir']) + '/variant_calling/wgs/' + str(s) + '_wgs_panel.vcf.gz'
                                  for s in _bed_wgs_samples(wc.cohort_id, wc.bed_id)],
        bed         = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
    output:
        matrix  = _cohort_outdir + "/{bed_id}/output/cohort_qc/wgs_concordance/{bed_id}_wgs_concordance_matrix.tsv",
        summary = _cohort_outdir + "/{bed_id}/output/cohort_qc/wgs_concordance/{bed_id}_wgs_concordance_summary.tsv",
        heatmap = _cohort_outdir + "/{bed_id}/output/cohort_qc/wgs_concordance/{bed_id}_wgs_concordance_heatmap.pdf",
    params:
        samples     = lambda wc: _quoted(bed_samples(wc.cohort_id, wc.bed_id)),
        wgs_samples = lambda wc: _quoted(_bed_wgs_samples(wc.cohort_id, wc.bed_id)),
        min_dp      = config.get("wgs_concordance_min_dp", 20),
        min_sites   = config.get("wgs_concordance_min_sites", 20),
        threshold   = config.get("wgs_concordance_threshold", 0.8),
        outprefix   = lambda wc, output: output.summary[:-len("_summary.tsv")],
        script      = workflow.basedir + "/scripts/_10P_check_wgs_concordance.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_10P_check_wgs_concordance.log"
    shell:
        """
        mkdir -p $(dirname {output.summary}) $(dirname {log})
        python -u {params.script} \\
            --samples     {params.samples} \\
            --nanoTS-vcfs {input.nanots_vcfs} \\
            --wgs-samples {params.wgs_samples} \\
            --wgs-vcfs    {input.wgs_vcfs} \\
            --bed         {input.bed} \\
            --min-dp      {params.min_dp} \\
            --min-sites   {params.min_sites} \\
            --threshold   {params.threshold} \\
            --outprefix   {params.outprefix} \\
        2>&1 | tee {log}
        """


def _report_group(wc):
    return _group_id_from_ids(wc.cohort_id, wc.bed_id, SAMPLES[wc.sample]["sample_type"])


def _report_tx_args(wc):
    gid = _report_group(wc)
    if not ctx_enabled(gid):
        return "--tx-unavailable " + _quoted([ctx_unavailable_reason(gid)])[0]
    f = ctx_files(gid)
    return " ".join([
        "--tx-usage-matrix " + f["usage_matrix"],
        "--tx-gene-status " + f["gene_status"],
        "--tx-isoform-outliers " + f["isoform_outliers"],
        "--tx-landscape-outliers " + f["landscape_outliers"],
        "--tx-gtf " + f["gtf"],
        "--tx-method " + _cja_method_for_group(gid),
        "--tx-n-threshold " + str(_cja_n_threshold(gid)),
        "--tx-min-gene-reads " + str(config.get("cohort_tx_min_gene_reads", 50)),
        "--tx-min-usage " + str(config.get("cohort_tx_min_usage", 0.05)),
        "--tx-min-reads " + str(config.get("cohort_tx_min_reads", 10)),
    ])


rule _10Q_make_final_report:
    input:
        merged_hits   = _cohort_outdir + "/{bed_id}/output/default/merged_all_hits.tsv",
        bam           = lambda wc: SAMPLES[wc.sample]["bam"],
        mapping       = lambda wc: _sample_gene_bam_mapping(wc.sample) if flag("phase_reads") else [],
        bed           = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
        gtf           = config["annotation"],
        cja_manifests = lambda wc: _bed_cja_manifests(wc.cohort_id, wc.bed_id) if flag("cohort_junction_analysis") else [],
        on_target     = lambda wc: (_cohort_qc_file(wc, "on_target_rates/{b}_on_target_rates.tsv") if flag("qc") else []),
        read_attrs    = lambda wc: (_cohort_qc_file(wc, "read_attributes/{b}_read_attributes.tsv") if flag("qc") else []),
        flr_matrix    = lambda wc: (_cohort_qc_file(wc, "full_length_ratio/{b}_full_length_ratio_matrix.tsv") if flag("qc") else []),
        pca_coords    = lambda wc: (_cohort_qc_file(wc, "validate_sample_types/{b}_PCA_coords.tsv") if flag("qc") else []),
        coverage      = lambda wc: (_cohort_qc_file(wc, "avg_coverage/{b}_breadth_coverage_matrix.tsv") if flag("qc") else []),
        sample_cov    = lambda wc: ((SAMPLES[wc.sample]["outdir"] + "/qc/" + wc.sample + "_avg_coverage.tsv") if flag("qc") else []),
        expr_mats     = lambda wc: ([p for _, p in _bed_expression_matrices(wc)] if flag("gene_quantification") else []),
        wgs_concord   = _wgs_concordance_summary,
        hpo_terms     = lambda wc: hpo_terms_tsv(wc.sample) if sample_hpo(wc.sample) else [],
        tx_files      = lambda wc: (
            list(ctx_files(_report_group(wc)).values()) if ctx_enabled(_report_group(wc)) else []),
    output:
        pdf = _cohort_outdir + "/{bed_id}/output/default/final_reports/{sample}_report.pdf",
        html = _cohort_outdir + "/{bed_id}/output/default/final_reports/{sample}_review.html",
    params:
        genome           = config["genome"],
        mapping_arg      = lambda wc, input: ("--gene-bam-mapping " + str(input.mapping)) if input.mapping else "",
        qc_args          = lambda wc, input: " ".join(
            f"{opt} {val}" for opt, val in (("--on-target-tsv", input.on_target),
                                            ("--read-attributes-tsv", input.read_attrs),
                                            ("--flr-matrix", input.flr_matrix),
                                            ("--pca-coords", input.pca_coords),
                                            ("--coverage-matrix", input.coverage),
                                            ("--sample-coverage-tsv", input.sample_cov),
                                            ("--wgs-concordance-tsv", input.wgs_concord)) if val)
                                   + (f" --coverage-depth {config.get('merge_min_dp_snv', 50)}" if input.coverage else ""),
        cohort_manifests = lambda wc: _quoted(_bed_cohort_manifest_args(wc.cohort_id, wc.bed_id)) if flag("cohort_junction_analysis") else [],
        samples          = lambda wc: bed_samples(wc.cohort_id, wc.bed_id),
        sample_types     = lambda wc: [SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
        colors           = lambda wc: _quoted([sample_type_color(SAMPLES[s]["sample_type"]) for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        sample_bams      = lambda wc: _quoted([SAMPLES[s]["bam"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        expr_args        = lambda wc: (("--expression-matrices " + " ".join(_quoted(
            [t + ":" + p for t, p in _bed_expression_matrices(wc)]))) if flag("gene_quantification") else ""),
        cov_thr          = config["sample_coverage_threshold"],
        top_n            = config.get("hit_report_top_n", 0),
        max_tier         = config.get("hit_report_max_tier", 5),
        max_reads        = config.get("hit_report_max_reads", 1000),
        hpo_arg          = lambda wc, input: ("--hpo-terms-tsv " + str(input.hpo_terms)) if input.hpo_terms else "",
        tx_args          = _report_tx_args,
        script           = workflow.basedir + "/scripts/_10Q_make_final_report.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/logs/final_reports/{sample}_10Q_make_final_report.log"
    shell:
        """
        mkdir -p $(dirname {output.pdf}) $(dirname {log})
        python -u {params.script} \\
            --sample             {wildcards.sample} \\
            --hits-tsv           {input.merged_hits} \\
            --bam                {input.bam} \\
            {params.mapping_arg} \\
            {params.qc_args} \\
            {params.hpo_arg} \\
            {params.tx_args} \\
            --bed                {input.bed} \\
            --gtf                {input.gtf} \\
            --genome             {params.genome} \\
            --cohort-manifests   {params.cohort_manifests} \\
            --samples            {params.samples} \\
            --sample-types       {params.sample_types} \\
            --colors             {params.colors} \\
            --sample-bams        {params.sample_bams} \\
            {params.expr_args} \\
            --coverage-threshold {params.cov_thr} \\
            --top-n              {params.top_n} \\
            --max-tier           {params.max_tier} \\
            --max-reads          {params.max_reads} \\
            --outfile            {output.pdf} \\
            --html-out           {output.html} \\
            --previous-reviews   "$(dirname {output.html})/{wildcards.sample}_reviews.tsv" \\
        2>&1 | tee {log}
        """
