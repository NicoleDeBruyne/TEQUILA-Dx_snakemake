
# Compare each sample's junction usage against the rest of its own group (rather than GTEx)
_cohort_outdir = config["output_dir"] + "/{cohort_id}"

def _group_gene_bam_mapping_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/phased_reads/' + str(s) + '_gene_bam_mapping_file.tsv') for s in GROUPS[group_id]]


_cja_manifest_path = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{bed_id}_{sample_type}_gene_manifest.tsv"


def _cja_method(group_id):
    """Return 'beta_binomial' or 'modified_zscore' for the given group."""
    # Use beta-binomial when the group has enough samples (same logic as _cja_thr_flag)
    flag = _cja_thr_flag(group_id)
    if "--bb-thresholds" in flag:
        return "beta_binomial"
    return "modified_zscore"


# Compute per-gene, per-sample junction coverage/usage metrics across the whole group
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
        min_reads   = config["cohort_jxn_min_reads"],
        script      = workflow.basedir + "/scripts/cohort_junction_analysis.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "cohort_junction_analysis", config["threads"])
    resources:
        mem_mb     = lambda wc, attempt: attempt * 1024 * max(32, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime    = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/cohort_junction_analysis.log"
    shell:
        """
        mkdir -p $(dirname {output.cohort_mapping})
        mkdir -p {params.raw_outdir}
        mkdir -p $(dirname {log})

        echo -e "gene\\tsample\\tbulk_bam\\thap1_bam\\thap2_bam" > {output.cohort_mapping}
        for f in {input.mapping_files}; do
            tail -n +2 "$f" | awk -F'\\t' 'BEGIN {{OFS="\\t"}} {{print $3, $1, $4, $5, $6}}'
        done >> {output.cohort_mapping}

        python -u {params.script} \\
            --mapping-file  {output.cohort_mapping} \\
            --bed           {input.bed} \\
            --outdir        {params.raw_outdir} \\
            --manifest      {output.manifest} \\
            --note          {output.note} \\
            --genome        {params.genome} \\
            --min-jxn-reads {params.min_reads} \\
            --threads       {threads} \\
        2>&1 | tee {log}
        """


# Fit distributions, run statistical tests, FDR correction, write combined scored TSV
_cja_scoring_sentinel = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{cohort_id}_{bed_id}_{sample_type}_scoring.done"

rule _8B_fit_and_score_cohort_junctions:
    input:
        manifest = _cja_manifest_path,
        bed      = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
    output:
        scored_sentinel = _cja_scoring_sentinel,
    params:
        outprefix  = lambda wc: (str(group_outdir(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))) + '/cohort_junction_analysis/' + str(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))),
        has_ipa    = "--has-ipa" if config["genome"] else "",
        method     = lambda wc: _cja_method(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        gtf        = config["annotation"],
        cov_thr    = config["sample_coverage_threshold"],
        n_thr      = lambda wc: _cja_n_threshold(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        script     = workflow.basedir + "/scripts/fit_and_score_cohort_junctions.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "fit_and_score_cohort_junctions", config["threads"])
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 4),
        runtime    = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/fit_and_score_cohort_junctions.log"
    shell:
        """
        mkdir -p $(dirname {log})

        python -u {params.script} \\
            --manifest           {input.manifest} \\
            --bed                {input.bed} \\
            --outprefix          {params.outprefix} \\
            {params.has_ipa} \\
            --method             {params.method} \\
            --gtf                {params.gtf} \\
            --coverage-threshold {params.cov_thr} \\
            --n-threshold        {params.n_thr} \\
            --threads            {threads} \\
        2>&1 | tee {log}
        """


# Read scored TSV, apply thresholds, write outlier outputs
rule _8C_identify_cohort_junction_outliers:
    input:
        scored_sentinel = _cja_scoring_sentinel,
        bed             = lambda wc: bed_path(wc.cohort_id, wc.bed_id),
    output:
        outliers          = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{cohort_id}_{bed_id}_{sample_type}_{thr_label}/{cohort_id}_{bed_id}_{sample_type}_outliers.tsv",
        outliers_filtered = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{cohort_id}_{bed_id}_{sample_type}_{thr_label}/{cohort_id}_{bed_id}_{sample_type}_outliers_filtered.tsv",
        outliers_alias    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_junction_analysis/{cohort_id}_{bed_id}_{sample_type}_{thr_label}/{cohort_id}_{bed_id}_{sample_type}_outliers_alias.tsv",
    params:
        outprefix    = lambda wc: (str(group_outdir(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))) + '/cohort_junction_analysis/' + str(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))),
        has_ipa      = "--has-ipa" if config["genome"] else "",
        thr_flag     = lambda wc: _cja_thr_flag(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        gtf          = config["annotation"],
        cov_thr      = config["sample_coverage_threshold"],
        phasing_thr  = config["cohort_jxn_phasing_threshold"],
        alias_args   = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script       = workflow.basedir + "/scripts/identify_cohort_junction_outliers.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "identify_cohort_junction_outliers", config["threads"])
    resources:
        mem_mb = lambda wc, attempt: attempt * 1024 * 16,
        runtime    = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/identify_cohort_junction_outliers_{thr_label}.log"
    shell:
        """
        mkdir -p $(dirname {log})

        python -u {params.script} \\
            --bed                {input.bed} \\
            --outprefix          {params.outprefix} \\
            {params.has_ipa} \\
            {params.thr_flag} \\
            --gtf                {params.gtf} \\
            --coverage-threshold {params.cov_thr} \\
            --phasing-threshold  {params.phasing_thr} \\
            --alias-map          {params.alias_args} \\
            --threads            {threads} \\
        2>&1 | tee {log}
        """
