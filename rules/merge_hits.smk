"""
rules/merge_hits.smk
Cross-sample merge-and-filter stage. Runs once per (bed, sample_type) group (see
GROUPS in the Snakefile). Every stage below calls exactly one standalone script
directly -- no script here shells out to another script (that used to happen
inside merge_and_filter_cohort_hits.py, which has been removed; Snakemake itself
now does that orchestration, as it should).

Stages, in dependency order:
  1. merge_group_variants   -- merge_and_filter_variants.py,          once per group
  2. merge_group_ase        -- merge_and_filter_ase_results.py,       once per group
  3. merge_group_junctions  -- merge_and_filter_junction_results.py,  once per (group, tissue)
  4. merge_sample_hits      -- split_group_hits_by_sample.py + merge_hits.py,
                                once per sample (splits this sample's rows out of the
                                three group-level files above, then merges them)
  5. concat_group_hits      -- plain awk concat (no script needed),   once per group
  6. plot_group_hits        -- plot_candidate_hits.py,                once per group
  7. final_merge            -- plain awk concat across sample_types sharing a BED panel
                                (unchanged from before)

Output is nested as {merged_outdir}/{bed_id}/{sample_type}/... rather than a flat
{merged_outdir}/{bed_id}_{sample_type}/..., so most rules use bed_id and sample_type
as two separate wildcards and reconstruct group_id via _group_id_from_ids() to look
back up into GROUPS/GROUP_* wherever needed.
"""

from math import ceil

def _group_variant_files(group_id):
    return [f"{SAMPLES[s]['outdir']}/variant_calling/compiled_variants/{s}_filtered_variants.tsv" for s in GROUPS[group_id]]

def _group_ase_files(group_id):
    return [f"{SAMPLES[s]['outdir']}/ase_analysis/{s}_binomial_ase_results.tsv" for s in GROUPS[group_id]]

def _group_tissue_samples(group_id, tissue):
    """Samples in this group that have the given tissue configured."""
    return [s for s in GROUPS[group_id] if tissue in sample_tissues(s)]

def _group_tissue_junction_files(group_id, tissue):
    return [f"{SAMPLES[s]['outdir']}/junction_analysis/gtex_{tissue}/{s}_gtex_{tissue}_all_junctions.tsv"
            for s in _group_tissue_samples(group_id, tissue)]

def _group_junction_outprefix(group_id, tissue):
    return f"{group_outdir(group_id)}/junction_analysis/gtex_{tissue}/outlier_junctions_gtex_{tissue}"

def _group_junction_final_path(group_id, tissue):
    """Static, Snakemake-tracked output path for a (group, tissue)'s merged
    junction hits -- what merge_group_junctions below copies its real (dynamically
    named, see _group_junction_source_glob) output to, so downstream rules never
    need to know merge_and_filter_junction_results.py's internal naming scheme."""
    return f"{_group_junction_outprefix(group_id, tissue)}_final.tsv"

def _group_junction_source_glob(group_id, tissue):
    """Shell glob matching whatever filename merge_and_filter_junction_results.py
    actually produces for this (group, tissue), given the flags always passed
    below (--event-types, --filter-by-cohort-IQR, --sample-number-threshold).
    Mirrors that script's own internal stage-naming logic (stage1 -> stage3
    [event] -> stage5 [cohortIQR, only if >=8 samples] -> stage6 [Nsamples]) --
    but only closely enough to build an unambiguous glob, not an exact filename,
    so it doesn't need to be kept in lockstep with that script's internals.
    """
    n_tissue_samples = len(_group_tissue_samples(group_id, tissue))
    n = ceil(n_tissue_samples * config["merge_jxn_sample_fraction"])
    base = (f"{_group_junction_outprefix(group_id, tissue)}_"
            f"{config['merge_jxn_coverage_threshold']}jxncov_"
            f"{config['merge_jxn_padj_threshold']}padj_"
            f"{config['merge_delta_psi_threshold']}deltaPSI_event")
    return f"{base}*_{n}samples.tsv"


# ---------------------------------------------------------------------------
# 1. Merge & filter variant calls across all samples in a group
# ---------------------------------------------------------------------------
rule merge_group_variants:
    input:
        variant_files = lambda wc: _group_variant_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        tsv = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/variant_calling/all_candidate_variants.tsv",
    params:
        group_id     = lambda wc: _group_id_from_ids(wc.bed_id, wc.sample_type),
        n            = lambda wc: ceil(len(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)])
                                        * config["merge_variant_sample_fraction"]),
        outprefix    = lambda wc, output: output.tsv[:-len(".tsv")],
        num_callers_snv   = config["merge_num_callers_threshold_snv"],
        num_callers_indel = config["merge_num_callers_threshold_indel"],
        min_dp_snv        = config["merge_min_dp_snv"],
        min_dp_indel      = config["merge_min_dp_indel"],
        script       = workflow.basedir + "/scripts/merge_and_filter_variants.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "merge_group_variants", 1)
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "merge_group_variants", 8)),
        runtime = config["time"],
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/merge_group_variants.log"
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
            --plot \\
            --plot-variant-type SNV indel \\
            --title "{params.group_id} Variant Counts" \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 2. Merge & filter ASE results across all samples in a group
# ---------------------------------------------------------------------------
rule merge_group_ase:
    input:
        ase_files = lambda wc: _group_ase_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        tsv = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/ase_analysis/outlier_ase.tsv",
    params:
        group_id     = lambda wc: _group_id_from_ids(wc.bed_id, wc.sample_type),
        n            = lambda wc: ceil(len(GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)])
                                        * config["merge_ase_sample_fraction"]),
        outprefix    = lambda wc, output: output.tsv[:-len(".tsv")],
        min_hap_ratio       = config["merge_min_haplotype_ratio"],
        delta_hap_ratio_thr = config["merge_delta_haplotype_ratio_threshold"],
        ase_padj_thr        = config["merge_ase_padj_threshold"],
        script       = workflow.basedir + "/scripts/merge_and_filter_ase_results.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "merge_group_ase", 1)
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "merge_group_ase", 8)),
        runtime = config["time"],
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/merge_group_ase.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles {input.ase_files} \\
            --outprefix {params.outprefix} \\
            --min-haplotype-ratio {params.min_hap_ratio} \\
            --delta-haplotype-ratio-threshold {params.delta_hap_ratio_thr} \\
            --padj-threshold {params.ase_padj_thr} \\
            --plot \\
            --title "{params.group_id}: Number of Genes with Allele-specific Expression by Sample" \\
            --sample-number-threshold {params.n} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 3. Merge & filter outlier junctions, once per (group, tissue)
# ---------------------------------------------------------------------------
rule merge_group_junctions:
    input:
        junction_files = lambda wc: _group_tissue_junction_files(
            _group_id_from_ids(wc.bed_id, wc.sample_type), wc.tissue),
    output:
        # Static, wildcard-only path -- Snakemake's output: can't be a callable
        # (only input:/params: can), so the dynamically-named file the script
        # actually produces gets cp'd here at the end of the shell block below.
        tsv = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/junction_analysis/gtex_{{tissue}}/outlier_junctions_gtex_{{tissue}}_final.tsv",
    params:
        group_id  = lambda wc: _group_id_from_ids(wc.bed_id, wc.sample_type),
        n         = lambda wc: ceil(len(_group_tissue_samples(_group_id_from_ids(wc.bed_id, wc.sample_type), wc.tissue))
                                     * config["merge_jxn_sample_fraction"]),
        outprefix = lambda wc: _group_junction_outprefix(_group_id_from_ids(wc.bed_id, wc.sample_type), wc.tissue),
        source_glob = lambda wc: _group_junction_source_glob(_group_id_from_ids(wc.bed_id, wc.sample_type), wc.tissue),
        jxn_cov_thr   = config["merge_jxn_coverage_threshold"],
        jxn_padj_thr  = config["merge_jxn_padj_threshold"],
        delta_psi_thr = config["merge_delta_psi_threshold"],
        script    = workflow.basedir + "/scripts/merge_and_filter_junction_results.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.bed_id, wc.sample_type), "merge_group_junctions", 1)
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(
            _group_id_from_ids(wc.bed_id, wc.sample_type), "merge_group_junctions", 8)),
        runtime = config["time"],
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/merge_group_junctions_{{tissue}}.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles {input.junction_files} \\
            --outprefix {params.outprefix} \\
            --jxn-coverage-threshold {params.jxn_cov_thr} \\
            --padj-threshold {params.jxn_padj_thr} \\
            --delta-PSI-threshold {params.delta_psi_thr} \\
            --event-types exon_skipping exon_inclusion alt_ss1 alt_ss2 \\
            --sample-number-threshold {params.n} \\
            --filter-by-cohort-IQR \\
            --plot \\
            --title "{params.group_id}: Number of Genes with Outlier Junctions by Sample" \\
        2>&1 | tee {log}
        SRC=$(ls {params.source_glob} 2>/dev/null | head -1)
        if [ -z "$SRC" ]; then
            echo "ERROR: no output file matched glob {params.source_glob}" | tee -a {log} >&2
            exit 1
        fi
        cp "$SRC" {output.tsv}
        echo "Copied $SRC -> {output.tsv}" >> {log}
        """


# ---------------------------------------------------------------------------
# 4. Split this sample's rows out of the group-level merged variant/ASE/
#    junction files, then merge them into one candidate-hits table -- one
#    job per sample. Both scripts called here are fully standalone (no
#    script calls another script); this rule just runs them in sequence via
#    a plain shell block, which is where that kind of two-step orchestration
#    belongs.
#
#    (Combined into one rule, rather than a separate group-level "split"
#    rule producing per-sample outputs, because classic/non-checkpoint
#    Snakemake output: must be statically resolvable from wildcards -- it
#    can't be a runtime-computed list of per-sample paths. Each sample's
#    split step re-reads the group's already-small, already-filtered merged
#    files, which is cheap.)
# ---------------------------------------------------------------------------
rule merge_sample_hits:
    input:
        variant_tsv    = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/variant_calling/all_candidate_variants.tsv",
        ase_tsv        = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/ase_analysis/outlier_ase.tsv",
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.bed_id, wc.sample_type), t)
            for t in group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type))
        ],
    output:
        # Static path -- {sample} is an extra wildcard here (alongside bed_id/
        # sample_type) purely so this stays a static, Snakemake-legal output
        # pattern; only ever requested for (bed_id, sample_type, sample)
        # combinations that are actually valid per GROUPS (see concat_group_hits).
        tsv = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/by_sample/{{sample}}_hits.tsv",
    params:
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        stub_dir     = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/by_sample",
        omim         = config["omim_file"],
        split_script = workflow.basedir + "/scripts/split_group_hits_by_sample.py",
        merge_script = workflow.basedir + "/scripts/merge_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/{{sample}}_merge_hits.log"
    shell:
        """
        mkdir -p {params.stub_dir} $(dirname {log})
        python -u {params.split_script} \\
            --variant-tsv {input.variant_tsv} \\
            --ase-tsv     {input.ase_tsv} \\
            --tissues     {params.tissues} \\
            --junction-files {input.junction_files} \\
            --samples     {wildcards.sample} \\
            --outdir      {params.stub_dir} \\
        2>&1 | tee {log}
        python -u {params.merge_script} \\
            --outfile      {output.tsv} \\
            --sample-name  {wildcards.sample} \\
            --variant-hits {params.stub_dir}/{wildcards.sample}_variant_hits.tsv \\
            --ase-hits     {params.stub_dir}/{wildcards.sample}_ase_hits.tsv \\
            --junction-hits {params.stub_dir}/{wildcards.sample}_junction_hits.tsv \\
            --omim         {params.omim} \\
        2>&1 | tee -a {log}
        """


# ---------------------------------------------------------------------------
# 5. Concatenate every sample's hits.tsv into one group-level all_hits.tsv.
#    Plain awk header-dedup, same pattern as the final_merge rule below --
#    no script needed for a straight concat.
# ---------------------------------------------------------------------------
rule concat_group_hits:
    input:
        sample_hits = lambda wc: [
            f"{group_outdir(_group_id_from_ids(wc.bed_id, wc.sample_type))}/merged_hits/by_sample/{s}_hits.tsv"
            for s in GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)]
        ],
    output:
        all_hits = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/all_hits.tsv",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * 2),
        runtime = 60,
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/concat_group_hits.log"
    shell:
        """
        mkdir -p $(dirname {output.all_hits}) $(dirname {log})
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.sample_hits} > {output.all_hits} 2> {log}
        echo "Finished concat to {output.all_hits}." >> {log}
        """


# ---------------------------------------------------------------------------
# 6. Plot candidate hits for a group (calls plot_candidate_hits.py directly).
# ---------------------------------------------------------------------------
rule plot_group_hits:
    input:
        all_hits = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/all_hits.tsv",
    output:
        pathogenic  = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/genes_with_pathogenic_variant.pdf",
        ase         = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/genes_with_ASE.pdf",
        junction    = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/genes_with_outlier_junction.pdf",
        dysreg      = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits/genes_with_RNA_dysregulation.pdf",
    params:
        outdir = f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/merged_hits",
        script = workflow.basedir + "/scripts/plot_candidate_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = 60,
    log:
        f"{config['merged_outdir']}/{{bed_id}}/{{sample_type}}/logs/plot_group_hits.log"
    shell:
        """
        python -u {params.script} \\
            --infile {input.all_hits} \\
            --outdir {params.outdir} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 7. Final merge across all sample types sharing a BED panel (unchanged).
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
