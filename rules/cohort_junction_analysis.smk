"""
rules/cohort_junction_analysis.smk
Cohort-level splice junction outlier analysis: for each (bed, sample_type) group
(see GROUPS in the Snakefile), combines every sample's own
{sample}_gene_bam_mapping_file.tsv (written by phase_reads.py -- see
rules/phase_reads.smk) into a single cohort-wide gene/sample/BAM mapping file,
then runs run_cohort_junction_outlier_analysis.py to fit per-junction Beta
distributions across the whole cohort's bulk BAMs and flag per-sample outliers
(rather than against GTEx reference tissues, as rules/junction_analysis.smk does).

Only depends on phasing having finished for every sample in the group -- not on
merge_hits or the per-sample GTEx-based junction analysis.
"""

def _group_gene_bam_mapping_files(group_id):
    return [f"{SAMPLES[s]['outdir']}/phased_reads/{s}_gene_bam_mapping_file.tsv" for s in GROUPS[group_id]]


rule run_cohort_junction_outlier_analysis:
    input:
        mapping_files = lambda wc: _group_gene_bam_mapping_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        bed           = lambda wc: bed_path(wc.bed_id),
    output:
        # Rebuilt from every sample's own gene_bam_mapping_file.tsv -- one
        # "gene\tsample\tbulk_bam\thap1_bam\thap2_bam" row per (sample, gene).
        cohort_mapping = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/cohort_junction_analysis/{{bed_id}}_{{sample_type}}_gene_bam_mapping_file.tsv",
        # Path mirrors run_cohort_junction_outlier_analysis.py's own naming
        # (outdir/{prefix_name}_padj{padj}_delta{delta}/{prefix_name}_outliers.tsv)
        # for the single padj:delta pair passed via --thresholds below.
        outliers = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/cohort_junction_analysis/{{bed_id}}_{{sample_type}}_padj{config['padj_threshold']}_delta{config['delta_psi_threshold']}/{{bed_id}}_{{sample_type}}_outliers.tsv",
    params:
        outprefix   = lambda wc: f"{group_outdir(_group_id_from_ids(wc.bed_id, wc.sample_type))}/cohort_junction_analysis/{_group_id_from_ids(wc.bed_id, wc.sample_type)}",
        genome      = config["genome"],
        gtf         = config["annotation"],
        cov_thr     = config["sample_coverage_threshold"],
        phasing_thr = config["cohort_jxn_phasing_threshold"],
        n_thr       = config["cohort_jxn_n_threshold"],
        min_reads   = config["cohort_jxn_min_reads"],
        padj_thr    = config["padj_threshold"],
        dpsi_thr    = config["delta_psi_threshold"],
        script      = workflow.basedir + "/scripts/run_cohort_junction_outlier_analysis.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "run_cohort_junction_outlier_analysis", config["threads"])
    resources:
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "run_cohort_junction_outlier_analysis", threads * 4)),
        runtime    = config["time"],
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/cohort_junction_outlier_analysis.log"
    shell:
        """
        mkdir -p $(dirname {output.cohort_mapping})
        mkdir -p $(dirname {log})

        echo -e "gene\\tsample\\tbulk_bam\\thap1_bam\\thap2_bam" > {output.cohort_mapping}
        for f in {input.mapping_files}; do
            tail -n +2 "$f" | awk -F'\\t' 'BEGIN {{OFS="\\t"}} {{print $3, $1, $4, $5, $6}}'
        done >> {output.cohort_mapping}

        python -u {params.script} \\
            --mapping-file       {output.cohort_mapping} \\
            --bed                {input.bed} \\
            --outprefix          {params.outprefix} \\
            --genome             {params.genome} \\
            --gtf                {params.gtf} \\
            --coverage-threshold {params.cov_thr} \\
            --phasing-threshold  {params.phasing_thr} \\
            --n-threshold        {params.n_thr} \\
            --min-jxn-reads      {params.min_reads} \\
            --thresholds         {params.padj_thr}:{params.dpsi_thr} \\
            --threads            {threads} \\
        2>&1 | tee {log}
        """
