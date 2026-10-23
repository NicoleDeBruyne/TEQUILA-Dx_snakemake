_cohort_outdir = config["output_dir"] + "/{cohort_id}"


def _amalgam_group_dir(cohort_id, bed_id, sample_type):
    return config["output_dir"] + "/" + str(cohort_id) + "/" + str(bed_id) + "/output/sample_types/" + str(sample_type) + "/output/gene_quantification/by_amalgam"


def _amalgam_stringtie_gtf(sample):
    return str(SAMPLES[sample]['outdir']) + "/gene_quantification/amalgam/" + str(sample) + "_stringtie.gtf"


def _amalgam_quantification_tsv(cohort_id, bed_id, sample_type, sample):
    return _amalgam_group_dir(cohort_id, bed_id, sample_type) + "/quantification/" + sample + "_transcript_quantification.tsv"


rule _9A_amalgam_merge_gtfs:
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
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_9A_amalgam_merge_gtfs.log"
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


rule _9B_amalgam_build_transcriptome:
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
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_9B_amalgam_build_transcriptome.log"
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


rule _9C_amalgam_annotate_orf:
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
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_9C_amalgam_annotate_orf.log"
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


rule _9D_amalgam_quantify_transcripts:
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
        mem_mb  = lambda wc, attempt: attempt * 1024 * 16,
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_{sample}_9D_amalgam_quantify_transcripts.log"
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


rule _9E_amalgam_aggregate_matrices:
    input:
        tsvs = lambda wc: [
            _amalgam_quantification_tsv(wc.cohort_id, wc.bed_id, wc.sample_type, s)
            for s in GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]
        ],
        bed = lambda wc: SAMPLES[GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)][0]]["bed"],
        gtf = config["annotation"],
    output:
        gene_matrix_raw_all    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_raw_all_genes.tsv",
        gene_matrix_raw        = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/gene_amalgam_matrix_raw.tsv",
        tx_matrix_raw_all      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/transcript_amalgam_matrix_raw_all_genes.tsv",
        tx_matrix_raw          = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/transcript_amalgam_matrix_raw.tsv",
    params:
        samples    = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        outdir     = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/_9E_amalgam_aggregate_matrices.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_9E_amalgam_aggregate_matrices.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.tsvs} \\
            --samples   {params.samples} \\
            --outdir    {params.outdir} \\
            --bed       {input.bed} \\
            --gtf       {input.gtf} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _9F_fit_and_score_transcript_usage:
    input:
        tx_matrix = rules._9E_amalgam_aggregate_matrices.output.tx_matrix_raw,
    output:
        isoform_scores   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_transcript_usage_scores.tsv.gz",
        landscape_scores = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_gene_isoform_landscape_scores.tsv.gz",
        usage_matrix     = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_transcript_usage_matrix.tsv",
        gene_status      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_transcript_usage_gene_status.tsv",
    params:
        outprefix      = lambda wc: _ctx_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/" + wc.sample_type,
        method         = _ctx_method,
        n_thr          = lambda wc: _cja_n_threshold(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        min_gene_reads = config.get("cohort_tx_min_gene_reads", 50),
        min_usage      = config.get("cohort_tx_min_usage", 0.05),
        min_reads      = config.get("cohort_tx_min_reads", 10),
        rescale        = config.get("psi_rescale_factor", 0.001),
        dm_sims        = config.get("cohort_tx_dm_simulations", 10000),
        script         = workflow.basedir + "/scripts/_9F_fit_and_score_transcript_usage.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "cohort_tx", 4)
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_9F_fit_and_score_transcript_usage.log"
    shell:
        """
        mkdir -p $(dirname {output.isoform_scores}) $(dirname {log})
        python -u {params.script} \\
            --tx-matrix      {input.tx_matrix} \\
            --outprefix      {params.outprefix} \\
            --method         {params.method} \\
            --n-threshold    {params.n_thr} \\
            --min-gene-reads {params.min_gene_reads} \\
            --min-usage      {params.min_usage} \\
            --min-reads      {params.min_reads} \\
            --rescale-factor {params.rescale} \\
            --dm-simulations {params.dm_sims} \\
            --threads        {threads} \\
        2>&1 | tee {log}
        """


rule _9G_identify_transcript_usage_outliers:
    input:
        isoform_scores   = rules._9F_fit_and_score_transcript_usage.output.isoform_scores,
        landscape_scores = rules._9F_fit_and_score_transcript_usage.output.landscape_scores,
    output:
        isoform_outliers   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_transcript_usage_outliers_{tx_thr}.tsv",
        landscape_outliers = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_gene_isoform_landscape_outliers_{tx_thr}.tsv",
        isoform_alias      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_transcript_usage_outliers_{tx_thr}_alias.tsv",
        landscape_alias    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/cohort_transcript_analysis/{sample_type}_gene_isoform_landscape_outliers_{tx_thr}_alias.tsv",
    wildcard_constraints:
        tx_thr = "(z|padj)[0-9.]+_delta[0-9.]+_tvd[0-9.]+",
    params:
        method      = _ctx_method,
        filter_args = _ctx_filter_args,
        alias_args  = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script      = workflow.basedir + "/scripts/_9G_identify_transcript_usage_outliers.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 8,
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_9G_identify_transcript_usage_outliers_{tx_thr}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --isoform-scores    {input.isoform_scores} \\
            --landscape-scores  {input.landscape_scores} \\
            --method            {params.method} \\
            {params.filter_args} \\
            --isoform-outfile   {output.isoform_outliers} \\
            --landscape-outfile {output.landscape_outliers} \\
            --isoform-alias-outfile   {output.isoform_alias} \\
            --landscape-alias-outfile {output.landscape_alias} \\
            --alias-map         {params.alias_args} \\
        2>&1 | tee {log}
        """
