def _softclip_alu_arg(wc):
    path = config.get("softclip_alu_annotation") or ""
    return ("--alu " + shlex.quote(str(path))) if path else ""


def _softclip_genes_arg(wc):
    genes = config.get("softclip_genes") or []
    if isinstance(genes, str):
        genes = [g for g in re.split(r"[,;\s]+", genes) if g]
    return ("--genes " + " ".join(shlex.quote(str(g)) for g in genes)) if genes else ""


def _softclip_mapping(wc):
    return _sample_gene_bam_mapping(wc.sample) if flag("phase_reads") else []


rule _11A_classify_read_ends:
    input:
        bam     = lambda wc: SAMPLES[wc.sample]["bam"],
        bed     = lambda wc: SAMPLES[wc.sample]["bed"],
        gtf     = config["annotation"],
        mapping = _softclip_mapping,
    output:
        clusters = "{outdir}/softclip_analysis/{sample}_softclip_clusters.tsv",
        counts   = "{outdir}/softclip_analysis/{sample}_read_end_counts.tsv",
    params:
        genome      = config["genome"],
        alu_arg     = _softclip_alu_arg,
        genes_arg   = _softclip_genes_arg,
        mapping_arg = lambda wc, input: ("--gene-bam-mapping " + str(input.mapping)) if input.mapping else "",
        min_clip    = config.get("softclip_min_clip", 50),
        cons_len    = config.get("softclip_consensus_len", 50),
        script      = workflow.basedir + "/scripts/_11A_classify_read_ends.py",
    threads: lambda wc: _rule_threads(wc, "softclip")
    resources:
        mem_mb  = lambda wc, threads, attempt: max(8192, int(attempt * threads * 2 * 1024)),
        runtime = config["time"],
    log:
        "{outdir}/../logs/{sample}_11A_classify_read_ends.log"
    shell:
        """
        mkdir -p $(dirname {output.clusters}) $(dirname {log})
        python -u {params.script} \\
            --sample           {wildcards.sample} \\
            --bam              {input.bam} \\
            --bed              {input.bed} \\
            --gtf              {input.gtf} \\
            --genome           {params.genome} \\
            {params.alu_arg} \\
            {params.genes_arg} \\
            {params.mapping_arg} \\
            --min-clip         {params.min_clip} \\
            --consensus-len    {params.cons_len} \\
            --threads          {threads} \\
            --clusters-outfile {output.clusters} \\
            --counts-outfile   {output.counts} \\
        2>&1 | tee {log}
        """


def _softclip_group_samples(wc):
    return GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]


def _softclip_dir():
    return _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/softclip_analysis"


rule _11B_find_softclip_outliers:
    input:
        clusters = lambda wc: [str(SAMPLES[s]["outdir"]) + "/softclip_analysis/" + str(s) + "_softclip_clusters.tsv"
                               for s in _softclip_group_samples(wc)],
    output:
        outliers = _softclip_dir() + "/{sample_type}_softclip_outliers.tsv",
        all_tsv  = _softclip_dir() + "/{sample_type}_softclip_clusters_all.tsv",
        pdf      = _softclip_dir() + "/{sample_type}_softclip_outliers.pdf",
    params:
        samples         = _softclip_group_samples,
        min_reads       = config.get("softclip_min_reads", 10),
        min_fraction    = config.get("softclip_min_fraction", 0.2),
        min_identity    = config.get("softclip_min_identity", 0.75),
        present_frac    = config.get("softclip_present_fraction", 0.05),
        max_sample_frac = config.get("softclip_sample_fraction", 0.05),
        title           = lambda wc: f"{wc.bed_id} {wc.sample_type}: unexplained soft-clip clusters",
        alias_args      = lambda wc: _quoted(alias_map_args(_softclip_group_samples(wc))),
        script          = workflow.basedir + "/scripts/_11B_find_softclip_outliers.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * 4,
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample_type}_11B_find_softclip_outliers.log"
    shell:
        """
        mkdir -p $(dirname {output.outliers}) $(dirname {log})
        python -u {params.script} \\
            --samples             {params.samples} \\
            --cluster-files       {input.clusters} \\
            --min-reads           {params.min_reads} \\
            --min-fraction        {params.min_fraction} \\
            --min-identity        {params.min_identity} \\
            --present-fraction    {params.present_frac} \\
            --max-sample-fraction {params.max_sample_frac} \\
            --title               {params.title:q} \\
            --all-outfile         {output.all_tsv} \\
            --outliers-outfile    {output.outliers} \\
            --plot-outfile        {output.pdf} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """
