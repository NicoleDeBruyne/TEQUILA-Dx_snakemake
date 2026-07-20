"""
rules/merge_hits.smk
1. merge_hits_by_sample_type – cross-sample merge-and-filter stage. Runs once per
                                (bed, sample_type) group (see GROUPS in the
                                Snakefile): merges & filters variant calls, ASE
                                results, and outlier splice junctions across every
                                sample in the group, then merges those three hit
                                types per sample into one cohort-level
                                candidate-hits table.
2. final_merge                – concatenates merge_hits_by_sample_type's
                                 all_hits.tsv across every sample_type sharing a
                                 BED panel into one merged_all_hits.tsv per panel.

Output of merge_hits_by_sample_type is nested as
{merged_outdir}/{bed_id}/{sample_type}/... rather than a flat
{merged_outdir}/{bed_id}_{sample_type}/..., so the rule uses bed_id and sample_type
as two separate wildcards and reconstructs group_id via _group_id_from_ids() to
look back up into GROUPS/GROUP_* wherever needed.

Replaces the old merge_candidate_hits.sh, which did the same thing manually via
bash + directory globbing after all per-sample jobs had finished.
"""

def _group_variant_files(group_id):
    return [f"{SAMPLES[s]['outdir']}/variant_calling/compiled_variants/{s}_filtered_variants.tsv" for s in GROUPS[group_id]]

def _group_ase_files(group_id):
    return [f"{SAMPLES[s]['outdir']}/ase_analysis/{s}_binomial_ase_results.tsv" for s in GROUPS[group_id]]

def _group_junction_manifest(group_id):
    """One 'tissue|sample|filepath' entry per (sample, tissue) pair in this group."""
    entries = []
    for s in GROUPS[group_id]:
        od = SAMPLES[s]["outdir"]
        for t in sample_tissues(s):
            entries.append(f"{t}:::{s}:::{od}/junction_analysis/gtex_{t}/{s}_gtex_{t}_all_junctions.tsv")
    return entries

def _group_junction_files(group_id):
    """Flat list of junction file paths (for Snakemake dependency tracking)."""
    return [e.split(":::", 2)[2] for e in _group_junction_manifest(group_id)]


rule merge_hits_by_sample_type:
    input:
        variant_files  = lambda wc: _group_variant_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        ase_files      = lambda wc: _group_ase_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        junction_files = lambda wc: _group_junction_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        all_hits = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/all_hits.tsv",
    params:
        group_id             = lambda wc: _group_id_from_ids(wc.bed_id, wc.sample_type),
        outdir                = lambda wc: group_outdir(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        samples               = lambda wc: GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)],
        junction_manifest     = lambda wc: _group_junction_manifest(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        omim                 = config["omim_file"],
        num_callers_snv      = config["merge_num_callers_threshold_snv"],
        num_callers_indel    = config["merge_num_callers_threshold_indel"],
        min_dp_snv           = config["merge_min_dp_snv"],
        min_dp_indel         = config["merge_min_dp_indel"],
        variant_frac         = config["merge_variant_sample_fraction"],
        min_hap_ratio        = config["merge_min_haplotype_ratio"],
        delta_hap_ratio_thr  = config["merge_delta_haplotype_ratio_threshold"],
        ase_padj_thr         = config["merge_ase_padj_threshold"],
        ase_frac             = config["merge_ase_sample_fraction"],
        jxn_cov_thr          = config["merge_jxn_coverage_threshold"],
        jxn_padj_thr         = config["merge_jxn_padj_threshold"],
        delta_psi_thr        = config["merge_delta_psi_threshold"],
        jxn_frac             = config["merge_jxn_sample_fraction"],
        script               = workflow.basedir + "/scripts/merge_and_filter_cohort_hits.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "merge_hits_by_sample_type", 1)
    resources:
        mem_mb     = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "merge_hits_by_sample_type", 16)),
        runtime    = config["time"],
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/merge_hits.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --group-id                          {params.group_id} \\
            --outdir                            {params.outdir} \\
            --omim                              {params.omim} \\
            --samples                           {params.samples} \\
            --variant-files                     {input.variant_files} \\
            --ase-files                          {input.ase_files} \\
            --junction-manifest                  {params.junction_manifest} \\
            --num-callers-threshold-snv         {params.num_callers_snv} \\
            --num-callers-threshold-indel       {params.num_callers_indel} \\
            --min-dp-snv                        {params.min_dp_snv} \\
            --min-dp-indel                      {params.min_dp_indel} \\
            --variant-sample-fraction           {params.variant_frac} \\
            --min-haplotype-ratio               {params.min_hap_ratio} \\
            --delta-haplotype-ratio-threshold   {params.delta_hap_ratio_thr} \\
            --ase-padj-threshold                {params.ase_padj_thr} \\
            --ase-sample-fraction               {params.ase_frac} \\
            --jxn-coverage-threshold            {params.jxn_cov_thr} \\
            --jxn-padj-threshold                {params.jxn_padj_thr} \\
            --delta-psi-threshold               {params.delta_psi_thr} \\
            --jxn-sample-fraction                {params.jxn_frac} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Final merge across all sample types sharing a BED panel
# (moved here from validate_and_finalize.smk; renamed from final_merge_hits_by_bed)
# ---------------------------------------------------------------------------
rule final_merge:
    input:
        all_hits = lambda wc: [f"{group_outdir(gid)}/merged_hits/all_hits.tsv" for gid in BED_GROUPS[wc.bed_id]],
    output:
        merged = f"{config['merged_outdir']}/{{bed_id}}/merged_all_hits.tsv",
    threads: lambda wc: _group_threads(wc.bed_id, "final_merge", 1)
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(wc.bed_id, "final_merge", 2)),
        runtime = 60,
    log:
        f"{config['merged_outdir']}/{{bed_id}}/logs/{{bed_id}}_final_merge.log"
    shell:
        """
        mkdir -p $(dirname {log})
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.all_hits} > {output.merged} 2> {log}
        echo "Finished final merge to {output.merged}." >> {log}
        """
