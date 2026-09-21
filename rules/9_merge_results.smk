
from math import ceil

_cohort_outdir  = config["output_dir"] + "/{cohort_id}"
_variant_tsv    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_variant_calling/all_candidate_variants.tsv"
_ase_tsv        = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_ase_analysis/outlier_ase.tsv"
_junction_final = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/gtex_{tissue}/outlier_junctions_gtex_{tissue}_final.tsv"
_all_hits_tsv   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/all_hits.tsv"

def _group_variant_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/variant_calling/compiled_variants/' + str(s) + '_filtered_variants.tsv') for s in GROUPS[group_id]]

def _group_ase_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/ase_analysis/' + str(s) + '_binomial_ase_results.tsv') for s in GROUPS[group_id]]

def _group_tissue_samples(group_id, tissue):
    return [s for s in GROUPS[group_id] if tissue in sample_tissues(s)]

def _group_tissue_junction_files(group_id, tissue):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/gtex_' + str(tissue) + '/' + str(s) + '_gtex_' + str(tissue) + '_all_junctions.tsv')
            for s in _group_tissue_samples(group_id, tissue)]

def _group_junction_outprefix(group_id, tissue):
    return (str(group_outdir(group_id)) + '/merged_junction_analysis/gtex_' + str(tissue) + '/outlier_junctions_gtex_' + str(tissue))

def _group_junction_final_path(group_id, tissue):
    return (str(_group_junction_outprefix(group_id, tissue)) + '_final.tsv')

def _group_junction_source_glob(group_id, tissue):
    n_tissue_samples = len(_group_tissue_samples(group_id, tissue))
    n = ceil(n_tissue_samples * config["merge_jxn_sample_fraction"])
    base = ((str(_group_junction_outprefix(group_id, tissue)) + '_')
            + (str(config['merge_jxn_coverage_threshold']) + 'jxncov_')
            + (str(config['merge_jxn_padj_threshold']) + 'padj_')
            + (str(config['merge_delta_psi_threshold']) + 'deltaPSI_event'))
    return (str(base) + '*_' + str(n) + 'samples.tsv')



rule _9A_merge_group_variants:
    input:
        variant_files = lambda wc: _group_variant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        tsv = _variant_tsv,
    params:
        group_id     = lambda wc: _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type),
        n            = lambda wc: ceil(len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])
                                        * config["merge_variant_sample_fraction"]),
        outprefix    = lambda wc, output: output.tsv[:-len(".tsv")],
        num_callers_snv   = config["merge_num_callers_threshold_snv"],
        num_callers_indel = config["merge_num_callers_threshold_indel"],
        min_dp_snv        = config["merge_min_dp_snv"],
        min_dp_indel      = config["merge_min_dp_indel"],
        script       = workflow.basedir + "/scripts/merge_and_filter_variants.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "merge_group_variants", 1)
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_variants.log"
    shell:


rule _9B_merge_group_ase:
    input:
        ase_files = lambda wc: _group_ase_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        tsv = _ase_tsv,
    params:
        group_id     = lambda wc: _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type),
        n            = lambda wc: ceil(len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])
                                        * config["merge_ase_sample_fraction"]),
        outprefix    = lambda wc, output: output.tsv[:-len(".tsv")],
        min_hap_ratio       = config["merge_min_haplotype_ratio"],
        delta_hap_ratio_thr = config["merge_delta_haplotype_ratio_threshold"],
        ase_padj_thr        = config["merge_ase_padj_threshold"],
        script       = workflow.basedir + "/scripts/merge_and_filter_ase_results.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "merge_group_ase", 1)
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_ase.log"
    shell:


rule _9C_merge_group_junctions:
    input:
        junction_files = lambda wc: _group_tissue_junction_files(
            _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue),
    output:
        tsv = _junction_final,
    params:
        group_id  = lambda wc: _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type),
        n         = lambda wc: ceil(len(_group_tissue_samples(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue))
                                     * config["merge_jxn_sample_fraction"]),
        outprefix = lambda wc: _group_junction_outprefix(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue),
        source_glob = lambda wc: _group_junction_source_glob(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue),
        jxn_cov_thr   = config["merge_jxn_coverage_threshold"],
        jxn_padj_thr  = config["merge_jxn_padj_threshold"],
        delta_psi_thr = config["merge_delta_psi_threshold"],
        script    = workflow.basedir + "/scripts/merge_and_filter_junction_results.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "merge_group_junctions", 1)
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(_group_tissue_samples(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue)) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_junctions_{tissue}.log"
    shell:


rule _9D1_merge_group_hits_preliminary:
    input:
        variant_tsv    = _variant_tsv,
        ase_tsv        = _ase_tsv,
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), t)
            for t in group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ],
    output:
        all_hits = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/all_hits_preliminary.tsv",
    params:
        samples      = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        cohort_junction_tsv = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        script       = workflow.basedir + "/scripts/merge_group_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_hits_preliminary.log"
    shell:


rule _9D2_merge_group_hits_with_cohort_junctions:
    input:
        variant_tsv    = _variant_tsv,
        ase_tsv        = _ase_tsv,
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), t)
            for t in group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))
        ],
        cohort_junction_tsv = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        gene_expression_matrix = lambda wc: (
            config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
            + "/output/gene_quantification/by_assignment/gene_assignment_matrix_cptm.tsv"
        ) if config.get("gene_quantification") else [],
        gene_expression_matrix_motr = lambda wc: (
            config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
            + "/output/gene_quantification/by_assignment/gene_assignment_matrix_motr.tsv"
        ) if config.get("gene_quantification") else [],
        gene_expression_zscores = lambda wc: (
            config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
            + "/output/gene_quantification/by_assignment/gene_assignment_zscores_cptm.tsv"
        ) if config.get("gene_quantification") else [],
    output:
        all_hits = _all_hits_tsv,
    params:
        samples      = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        gene_expression_flag = lambda wc, input: (
            "--gene-expression-matrix " + str(input.gene_expression_matrix)
        ) if config.get("gene_quantification") else "",
        gene_expression_motr_flag = lambda wc, input: (
            "--gene-expression-matrix-motr " + str(input.gene_expression_matrix_motr)
        ) if config.get("gene_quantification") else "",
        gene_expression_zscore_flag = lambda wc, input: (
            "--gene-expression-zscores " + str(input.gene_expression_zscores)
            + " --gene-expression-outlier-threshold " + str(config.get("gene_outlier_zscore_threshold", 3.0))
        ) if config.get("gene_quantification") else "",
        script       = workflow.basedir + "/scripts/merge_group_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_hits.log"
    shell:


rule _9E_plot_group_hits:
    input:
        all_hits = lambda wc: _group_all_hits_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        pathogenic_bar  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_pathogenic_variant_barplot.pdf",
        pathogenic_box  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_pathogenic_variant_boxplot.pdf",
        ase_bar         = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_ASE_barplot.pdf",
        ase_box         = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_ASE_boxplot.pdf",
        junction_bar    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_outlier_junction_barplot.pdf",
        junction_box    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_outlier_junction_boxplot.pdf",
        dysreg_bar      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_RNA_dysregulation_barplot.pdf",
        dysreg_box      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_RNA_dysregulation_boxplot.pdf",
    params:
        outdir = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits",
        script = workflow.basedir + "/scripts/plot_candidate_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/plot_group_hits.log"
    shell:


def _group_all_hits_path(group_id):
    fname = "all_hits.tsv" if config.get("merge_results_include_cohort_junctions", True) else "all_hits_preliminary.tsv"
    return str(group_outdir(group_id)) + "/merged_hits/" + fname


rule _9F_final_merge:
    input:
        all_hits = lambda wc: [_group_all_hits_path(gid) for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]],
    output:
        merged          = _cohort_outdir + "/{bed_id}/output/merged_all_hits.tsv",
        simplified      = _cohort_outdir + "/{bed_id}/output/merged_all_hits_simplified.tsv",
        simplified_alias = _cohort_outdir + "/{bed_id}/output/merged_all_hits_simplified_alias.tsv",
    params:
        script     = workflow.basedir + "/scripts/simplify_all_hits.py",
        alias_args = lambda wc: _quoted(alias_map_args(bed_samples(wc.cohort_id, wc.bed_id))),
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "final_merge", 1)
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_final_merge.log"
    shell:


rule _9F_plot_hits_upset:
    input:
        all_hits = _cohort_outdir + "/{bed_id}/output/merged_all_hits.tsv",
    output:
        pdf_density     = _cohort_outdir + "/{bed_id}/output/hits_upset_density.pdf",
        pdf_density_log = _cohort_outdir + "/{bed_id}/output/hits_upset_density_log.pdf",
        tsv = _cohort_outdir + "/{bed_id}/output/hits_upset_counts.tsv",
    params:
        samples      = lambda wc: bed_samples(wc.cohort_id, wc.bed_id),
        sample_types = lambda wc: [SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)],
        outdir       = _cohort_outdir + "/{bed_id}/output",
        title        = lambda wc: f"{wc.bed_id}: Candidate Gene Hit Categories by Sample",
        script       = workflow.basedir + "/scripts/plot_hits_upset.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * 2),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_hits_upset.log"
    shell:



def _sample_qc_file(sample, name):
    return str(SAMPLES[sample]['outdir']) + '/qc/' + str(sample) + '_' + str(name) + '.tsv'

def _bed_qc_files(cohort_id, bed_id, name):
    return [_sample_qc_file(s, name) for s in bed_samples(cohort_id, bed_id)]


rule _9G_merge_on_target_rates:
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
        script     = workflow.basedir + "/scripts/plot_on_target_rates.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_on_target_rates.log"
    shell:


rule _9H_merge_read_attributes:
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
        script     = workflow.basedir + "/scripts/plot_read_attributes.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_read_attributes.log"
    shell:


def _group_junction_matrix_inputs(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/junction_counts/' + str(s) + '_junction_count_matrix.tsv')
            for s in GROUPS[group_id]]


rule _9I1_build_group_junction_matrix:
    input:
        matrices = lambda wc: _group_junction_matrix_inputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/junction_count_matrix.tsv",
    params:
        samples = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        script  = workflow.basedir + "/scripts/build_group_junction_matrix.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "build_group_junction_matrix", config["threads"])
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/group_junction_matrix.log"
    shell:


rule _9I2_validate_sample_types:
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
        heatmap_alias = _cohort_outdir + "/{bed_id}/output/cohort_qc/validate_sample_types/{bed_id}_distance_heatmap_alias.pdf",
    params:
        outprefix    = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/validate_sample_types/' + str(wc.bed_id)),
        ref_names    = lambda wc: _quoted(config["validate_ref_tissues"]),
        ref_colors   = lambda wc: _quoted(config["validate_ref_colors"]),
        query_names  = lambda wc: _quoted([GROUP_SAMPLE_TYPE[gid] for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]]),
        query_colors = lambda wc: _quoted([sample_type_color(GROUP_SAMPLE_TYPE[gid]) for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]]),
        alias_args   = lambda wc: _quoted(alias_map_args(bed_samples(wc.cohort_id, wc.bed_id))),
        script       = workflow.basedir + "/scripts/validate_sample_type.py",
    threads: lambda wc: _group_threads(str(wc.cohort_id) + "_" + str(wc.bed_id), "validate_sample_types", 1)
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(256, len(bed_samples(wc.cohort_id, wc.bed_id))),
        runtime = 1440,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_validate_sample_types.log"
    shell:


rule _9J_merge_full_length_ratio:
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
        script     = workflow.basedir + "/scripts/get_full_length_ratio.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(bed_samples(wc.cohort_id, wc.bed_id)) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_full_length_ratio.log"
    shell:



def _sample_quant_file(sample, name):
    return str(SAMPLES[sample]['outdir']) + '/gene_quantification/' + str(sample) + '_' + str(name) + '.tsv'

def _group_quant_files(group_id, name):
    return [_sample_quant_file(s, name) for s in GROUPS[group_id]]


def _quoted_outlier_args(cid=None):
    return (
        "--outlier-pseudocount " + str(config.get("gene_outlier_pseudocount", 1.0))
        + " --outlier-shrinkage-k " + str(config.get("gene_outlier_shrinkage_k", 10.0))
        + " --outlier-min-mad " + str(config.get("gene_outlier_min_mad", 0.1))
        + " --outlier-zscore-threshold " + str(config.get("gene_outlier_zscore_threshold", 3.0))
    )


rule _9K_merge_gene_count:
    input:
        infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_count"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_cptm.tsv",
        matrix_motr  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_motr.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_raw.tsv",
        zscores_cptm = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_zscores_cptm.tsv",
        zscores_motr = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_zscores_motr.tsv",
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_count/gene_count"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (read count)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        outlier_args = _quoted_outlier_args(),
        script     = workflow.basedir + "/scripts/quantify_gene_expression.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_count_quantification.log"
    shell:


rule _9L_merge_gene_coverage:
    input:
        infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_coverage"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_cptm.tsv",
        matrix_motr  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_motr.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_raw.tsv",
        zscores_cptm = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_zscores_cptm.tsv",
        zscores_motr = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_zscores_motr.tsv",
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_coverage/gene_coverage"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (max coverage)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        outlier_args = _quoted_outlier_args(),
        script     = workflow.basedir + "/scripts/quantify_gene_expression.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_coverage_quantification.log"
    shell:


rule _9M_merge_gene_by_assignment:
    input:
        infiles       = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_assignment"),
        stats_infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "read_outcomes"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_cptm.tsv",
        matrix_motr  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_motr.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw.tsv",
        matrix_raw_all_genes = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw_all_genes.tsv",
        assignment_stats     = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_read_outcomes.tsv",
        assignment_stats_pdf = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_read_outcomes.pdf",
        zscores_cptm = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_zscores_cptm.tsv",
        zscores_motr = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_zscores_motr.tsv",
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_assignment/gene_assignment"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (splice-site assignment)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        outlier_args = _quoted_outlier_args(),
        script     = workflow.basedir + "/scripts/quantify_gene_by_assignment.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_assignment_quantification.log"
    shell:


def _amalgam_group_dir(cohort_id, bed_id, sample_type):
    return config["output_dir"] + "/" + str(cohort_id) + "/" + str(bed_id) + "/output/sample_types/" + str(sample_type) + "/output/gene_quantification/by_amalgam"


def _amalgam_stringtie_gtf(sample):
    return str(SAMPLES[sample]['outdir']) + "/gene_quantification/amalgam/" + str(sample) + "_stringtie.gtf"


def _amalgam_quantification_tsv(cohort_id, bed_id, sample_type, sample):
    return _amalgam_group_dir(cohort_id, bed_id, sample_type) + "/quantification/" + sample + "_transcript_quantification.tsv"


rule _9N1_amalgam_merge_gtfs:
    input:
        annotation  = config["annotation"],
        sample_gtfs = lambda wc: [
            _amalgam_stringtie_gtf(s)
            for s in GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]
        ],
    output:
        combined_gtf = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/merged.combined.gtf",
        tracking      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/annotation/merged.tracking",
    params:
        amalgam_env = config["amalgam_dir"] + "/conda_env",
        outprefix   = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/annotation/merged",
        gtf_list    = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/annotation/gtf_list.tsv",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 32,
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_merge_gtfs.log"
    shell:


rule _9N2_amalgam_build_transcriptome:
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
        main_env       = workflow.basedir + "/envs/conda_env",
        amalgam_dir    = config["amalgam_dir"],
        merge_prefix   = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/annotation/merged",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(32, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_build_transcriptome.log"
    shell:


rule _9N3_amalgam_annotate_orf:
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
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_annotate_orf.log"
    shell:


rule _9N4_amalgam_quantify_transcripts:
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
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample}_amalgam_quantify_transcripts.log"
    shell:


rule _9N5_amalgam_aggregate_matrices:
    input:
        tsvs = lambda wc: [
            _amalgam_quantification_tsv(wc.cohort_id, wc.bed_id, wc.sample_type, s)
            for s in GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]
        ],
    output:
        transcript_matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_transcript_matrix.tsv",
        gene_matrix        = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_gene_matrix.tsv",
        gene_matrix_alias  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_gene_matrix_alias.tsv",
    params:
        samples    = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        outprefix  = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/quantification/gene_amalgam",
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/aggregate_amalgam_matrices.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_aggregate_matrices.log"
    shell:


rule _9N6_amalgam_normalize_matrix:
    input:
        gene_matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_gene_matrix.tsv",
        bed = lambda wc: SAMPLES[GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)][0]]["bed"],
        gtf = config["annotation"],
    output:
        matrix_raw  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_raw.tsv",
        matrix_cptm = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_cptm.tsv",
        matrix_motr = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_motr.tsv",
        zscores_cptm = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_zscores_cptm.tsv",
        zscores_motr = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_zscores_motr.tsv",
        matrix_cptm_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_cptm_alias.tsv",
        matrix_motr_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_motr_alias.tsv",
    params:
        outprefix  = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/gene_amalgam",
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (AMALGAM)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        outlier_args = _quoted_outlier_args(),
        script     = workflow.basedir + "/scripts/normalize_amalgam_matrix.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_normalize_matrix.log"
    shell:
