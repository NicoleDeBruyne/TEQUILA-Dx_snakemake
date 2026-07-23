"""
rules/7_cohort_junction_analysis.smk
Cohort-level splice junction outlier analysis: for each (bed, sample_type) group
(see GROUPS in the Snakefile), combines every sample's own
{sample}_gene_bam_mapping_file.tsv (written by phase_reads.py -- see
rules/3_phase_reads.smk) into a single cohort-wide gene/sample/BAM mapping file,
then computes per-junction, per-sample coverage/usage metrics and flags
per-sample outliers against the whole cohort's own bulk BAMs (rather than
against GTEx reference tissues, as rules/5_junction_analysis.smk does), using
either beta-binomial testing or a modified z-score.

Which method a given group uses is decided by _cja_method_for_group() in the
main Snakefile (config["cohort_jxn_method"]: "auto" by default -- beta_binomial
for groups with >= config["cohort_jxn_beta_min_samples"] bulk samples,
modified_zscore for smaller groups -- or force one method for every group via
"beta_binomial"/"modified_zscore"). Since a group's method can only be known
once GROUPS is built, and Snakemake output: paths must be static wildcard
templates (not a value chosen per-group at parse time), the outlier output
path below uses an extra {thr_label} wildcard rather than picking between two
different literal path templates; _7B's own params re-derive the method from
{bed_id}/{sample_type} directly (not by inspecting {thr_label}) to keep a
single source of truth with the Snakefile's target-building in all_outputs().

Split into two rules:
  1. _7A_cohort_junction_analysis          -- per-gene metric computation
                                             (BAM reading, no statistics),
                                             via scripts/cohort_junction_analysis.py.
                                             Writes one raw TSV per gene plus
                                             a manifest mapping every BED
                                             gene to its result path (or
                                             "None"). Groups with fewer than
                                             config["cohort_jxn_min_samples"]
                                             samples are skipped entirely
                                             (see output.note for why).
  2. _7B_identify_cohort_junction_outliers -- reads the manifest, fits a
                                             per-junction reference
                                             distribution per gene and scores
                                             every sample against it, then
                                             identifies + classifies outlier
                                             junctions against the configured
                                             threshold, via
                                             scripts/identify_cohort_junction_outliers.py.

Only depends on phasing having finished for every sample in the group -- not on
merge_hits or the per-sample GTEx-based junction analysis.
"""

# NOTE: output:/input:/log: path templates below use string concatenation,
# not f-strings, to combine a config value with a literal Snakemake
# wildcard placeholder like "{bed_id}" -- an f-string's "{{bed_id}}" escape
# (to produce a literal "{bed_id}") does not survive Snakemake's own rule
# parsing and raises a NameError at load time.
_merged_outdir = config["merged_outdir"]

def _group_gene_bam_mapping_files(group_id):
    return [f"{SAMPLES[s]['outdir']}/phased_reads/{s}_gene_bam_mapping_file.tsv" for s in GROUPS[group_id]]


_cja_manifest_path = _merged_outdir + "/{bed_id}/{sample_type}/cohort_junction_analysis/{bed_id}_{sample_type}_gene_manifest.tsv"


rule _7A_cohort_junction_analysis:
    input:
        mapping_files = lambda wc: _group_gene_bam_mapping_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        bed           = lambda wc: bed_path(wc.bed_id),
    output:
        # Rebuilt from every sample's own gene_bam_mapping_file.tsv -- one
        # "gene\tsample\tbulk_bam\thap1_bam\thap2_bam" row per (sample, gene).
        cohort_mapping = _merged_outdir + "/{bed_id}/{sample_type}/cohort_junction_analysis/{bed_id}_{sample_type}_gene_bam_mapping_file.tsv",
        # One row per gene in the BED file: gene name + path to that gene's
        # raw per-junction metrics TSV (under params.raw_outdir), or the
        # literal string "None" if the gene produced no output (including
        # when the whole group was skipped for having too few samples --
        # see output.note). Presence of this file marks the rule complete
        # for Snakemake's purposes.
        manifest = _cja_manifest_path,
        # Always written, whether the analysis ran or was skipped -- explains
        # which, and why. Check this file first if a group's outlier results
        # look unexpectedly empty.
        note = _merged_outdir + "/{bed_id}/{sample_type}/cohort_junction_analysis/{bed_id}_{sample_type}_note.txt",
    params:
        raw_outdir  = lambda wc: f"{group_outdir(_group_id_from_ids(wc.bed_id, wc.sample_type))}/cohort_junction_analysis/{_group_id_from_ids(wc.bed_id, wc.sample_type)}_raw",
        genome      = config["genome"],
        cov_thr     = config["sample_coverage_threshold"],
        phasing_thr = config["cohort_jxn_phasing_threshold"],
        min_reads   = config["cohort_jxn_min_reads"],
        min_samples = config["cohort_jxn_min_samples"],
        script      = workflow.basedir + "/scripts/cohort_junction_analysis.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "cohort_junction_analysis", config["threads"])
    resources:
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "cohort_junction_analysis", threads * 4)),
        runtime    = config["time"],
    log:
        _merged_outdir + "/{bed_id}/{sample_type}/logs/cohort_junction_analysis.log"
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
            --mapping-file       {output.cohort_mapping} \\
            --bed                {input.bed} \\
            --outdir             {params.raw_outdir} \\
            --manifest           {output.manifest} \\
            --note               {output.note} \\
            --min-samples        {params.min_samples} \\
            --genome             {params.genome} \\
            --coverage-threshold {params.cov_thr} \\
            --phasing-threshold  {params.phasing_thr} \\
            --min-jxn-reads      {params.min_reads} \\
            --threads            {threads} \\
        2>&1 | tee {log}
        """


rule _7B_identify_cohort_junction_outliers:
    input:
        manifest = _cja_manifest_path,
        bed      = lambda wc: bed_path(wc.bed_id),
    output:
        # {thr_label} is an extra wildcard (rather than a value chosen once
        # for the whole file, as the rest of this path is) since which
        # method -- and therefore which label, e.g. "padj0.05_delta0.1" vs
        # "z3.5" -- applies varies per group. Snakemake binds it from
        # whatever concrete path was actually requested (built by
        # all_outputs() via _cja_thr_label()); it isn't otherwise used
        # below, since params.thr_flag re-derives the method independently
        # from {bed_id}/{sample_type} rather than parsing this wildcard.
        outliers = _merged_outdir + "/{bed_id}/{sample_type}/cohort_junction_analysis/{bed_id}_{sample_type}_{thr_label}/{bed_id}_{sample_type}_outliers.tsv",
    params:
        outprefix = lambda wc: f"{group_outdir(_group_id_from_ids(wc.bed_id, wc.sample_type))}/cohort_junction_analysis/{_group_id_from_ids(wc.bed_id, wc.sample_type)}",
        has_ipa   = "--has-ipa" if config["genome"] else "",
        thr_flag  = lambda wc: _cja_thr_flag(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        gtf       = config["annotation"],
        cov_thr   = config["sample_coverage_threshold"],
        n_thr     = config["cohort_jxn_n_threshold"],
        script    = workflow.basedir + "/scripts/identify_cohort_junction_outliers.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "identify_cohort_junction_outliers", config["threads"])
    resources:
        mem_mb     = lambda wc, threads, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "identify_cohort_junction_outliers", threads * 4)),
        runtime    = config["time"],
    log:
        _merged_outdir + "/{bed_id}/{sample_type}/logs/identify_cohort_junction_outliers_{thr_label}.log"
    shell:
        """
        mkdir -p $(dirname {log})

        python -u {params.script} \\
            --manifest           {input.manifest} \\
            --bed                {input.bed} \\
            --outprefix          {params.outprefix} \\
            {params.has_ipa} \\
            {params.thr_flag} \\
            --gtf                {params.gtf} \\
            --coverage-threshold {params.cov_thr} \\
            --n-threshold        {params.n_thr} \\
            --threads            {threads} \\
        2>&1 | tee {log}
        """
