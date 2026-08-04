"""
rules/6_merge_hits.smk
Cross-sample merge-and-filter stage. Runs once per (bed, sample_type) group.
See docs/rules/6_merge_hits.md for the full stage list, dependency order, and
the {bed_id}/output/sample_types/{sample_type} wildcard convention used throughout
this file.

Per-sample merge is split into two mutually-exclusive rules:
  _6D1  -- doesn't depend on cohort_junction_analysis (rule 7) at all, so it
          can run as soon as variant/ASE/junction results are ready.
  _6D2 -- the fully-informed version: does depend on rule 7, and picks up
          cohort-comparison junction results once they're available.

config['merge_hits_include_cohort_junctions'] (default True) controls which
one _6E (and therefore everything chained after it -- _6F/_6G/
merged_all_hits.tsv) actually consumes -- see _group_sample_hits_files()
below, the single place that decision is made. Only the selected branch is
ever requested, so only it ever runs; the two are never both computed in
the same run. Flip the config value and rerun snakemake to switch which one
merged_all_hits.tsv is built from (this also determines whether that run
waits on rule 7 or not).
"""

from math import ceil

# NOTE: output:/input:/log: path templates below use string concatenation,
# not f-strings, to combine a config value with a literal Snakemake
# wildcard placeholder like "{bed_id}" -- an f-string's "{{bed_id}}" escape
# (to produce a literal "{bed_id}") does not survive Snakemake's own rule
# parsing and raises a NameError at load time. Paths reused across more than
# one rule (e.g. all_candidate_variants.tsv, an output of _6A and an input
# of _6D1) are also factored out here so both rules stay in sync.
_cohort_outdir  = config["output_dir"] + "/cohort"
_variant_tsv    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_variant_calling/all_candidate_variants.tsv"
_ase_tsv        = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_ase_analysis/outlier_ase.tsv"
_junction_final = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/gtex_{tissue}/outlier_junctions_gtex_{tissue}_final.tsv"
_all_hits_tsv   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/all_hits.tsv"

def _group_sample_hits_files(group_id):
    """Every sample's per-sample hits.tsv for a group, from whichever of
    _6D1/_6D2 config['merge_hits_include_cohort_junctions'] selects (see
    that key's doc comment in config.yaml). _6D1 and _6D2 are mutually
    exclusive by construction: this is the *only* place that decides which
    one's output _6E (and therefore everything chained after it --
    _6F/_6G/merged_all_hits.tsv) actually consumes, so only that one branch
    ever gets pulled into the DAG for a given run."""
    god = group_outdir(group_id)
    subdir = "by_sample" if config.get("merge_hits_include_cohort_junctions", True) else "by_sample_preliminary"
    return [(str(god) + "/merged_hits/" + subdir + "/" + str(s) + "_hits.tsv") for s in GROUPS[group_id]]

def _group_variant_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/variant_calling/compiled_variants/' + str(s) + '_filtered_variants.tsv') for s in GROUPS[group_id]]

def _group_ase_files(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/ase_analysis/' + str(s) + '_binomial_ase_results.tsv') for s in GROUPS[group_id]]

def _group_tissue_samples(group_id, tissue):
    """Samples in this group that have the given tissue configured."""
    return [s for s in GROUPS[group_id] if tissue in sample_tissues(s)]

def _group_tissue_junction_files(group_id, tissue):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/gtex_' + str(tissue) + '/' + str(s) + '_gtex_' + str(tissue) + '_all_junctions.tsv')
            for s in _group_tissue_samples(group_id, tissue)]

def _group_junction_outprefix(group_id, tissue):
    return (str(group_outdir(group_id)) + '/merged_junction_analysis/gtex_' + str(tissue) + '/outlier_junctions_gtex_' + str(tissue))

def _group_junction_final_path(group_id, tissue):
    """Static, Snakemake-tracked output path for a (group, tissue)'s merged
    junction hits -- copied from merge_and_filter_junction_results.py's
    dynamically-named output at the end of the shell block below."""
    return (str(_group_junction_outprefix(group_id, tissue)) + '_final.tsv')

def _group_junction_source_glob(group_id, tissue):
    """Shell glob matching whatever filename merge_and_filter_junction_results.py
    actually produces for this (group, tissue). Mirrors that script's
    internal stage-naming closely enough for an unambiguous glob, without
    needing to track its exact filename."""
    n_tissue_samples = len(_group_tissue_samples(group_id, tissue))
    n = ceil(n_tissue_samples * config["merge_jxn_sample_fraction"])
    base = ((str(_group_junction_outprefix(group_id, tissue)) + '_')
            + (str(config['merge_jxn_coverage_threshold']) + 'jxncov_')
            + (str(config['merge_jxn_padj_threshold']) + 'padj_')
            + (str(config['merge_delta_psi_threshold']) + 'deltaPSI_event'))
    return (str(base) + '*_' + str(n) + 'samples.tsv')


# ---------------------------------------------------------------------------
# 6A. Merge & filter variant calls across all samples in a group
# ---------------------------------------------------------------------------
rule _6A_merge_group_variants:
    input:
        variant_files = lambda wc: _group_variant_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        tsv = _variant_tsv,
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
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_variants.log"
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
# 6B. Merge & filter ASE results across all samples in a group
# ---------------------------------------------------------------------------
rule _6B_merge_group_ase:
    input:
        ase_files = lambda wc: _group_ase_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        tsv = _ase_tsv,
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
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_ase.log"
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
# 6C. Merge & filter outlier junctions, once per (group, tissue)
# ---------------------------------------------------------------------------
rule _6C_merge_group_junctions:
    input:
        junction_files = lambda wc: _group_tissue_junction_files(
            _group_id_from_ids(wc.bed_id, wc.sample_type), wc.tissue),
    output:
        # Static, wildcard-only path -- the dynamically-named file the
        # script actually produces gets cp'd here at the end of the shell block.
        tsv = _junction_final,
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
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_junctions_{tissue}.log"
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
            echo "WARNING: no output file matched glob {params.source_glob}" | tee -a {log} >&2
            exit 1
        fi
        cp "$SRC" {output.tsv}
        echo "Copied $SRC -> {output.tsv}" >> {log}
        """


# ---------------------------------------------------------------------------
# 6D1. Split this sample's rows out of the group-level merged variant/ASE/
#    junction files, then merge them into one candidate-hits table -- one
#    job per sample. Combined into one rule (rather than a separate
#    group-level "split" rule) because classic Snakemake output: must be
#    statically resolvable from wildcards, not a runtime-computed list of
#    per-sample paths. Each sample's split step re-reads the group's
#    already-small, already-filtered merged files, which is cheap.
#
# PRELIMINARY / fast path: deliberately does NOT depend on
# cohort_junction_analysis (rules/7_cohort_junction_analysis.smk), which can
# take a while. cohort-junction-tsv is still passed to the scripts below as
# a params (not input) path: if that file already happens to exist on disk
# (e.g. a previous run already completed cohort_junction_analysis for this
# group), it's picked up opportunistically; if not, merge_hits.py's existing
# graceful "no cohort data" fallback applies (same as its --omim handling).
# Writes to merged_hits/by_sample_preliminary/, NOT the canonical
# merged_hits/by_sample/ path _6D2 writes to -- see _6D2 below.
# Only runs at all when config['merge_hits_include_cohort_junctions'] is
# False, since that's what determines which of _6D1/_6D2 _6E requests -- see
# _group_sample_hits_files() near the top of this file.
# ---------------------------------------------------------------------------
rule _6D1_merge_sample_hits_preliminary:
    input:
        variant_tsv    = _variant_tsv,
        ase_tsv        = _ase_tsv,
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.bed_id, wc.sample_type), t)
            for t in group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type))
        ],
    output:
        # {sample} is an extra wildcard purely so this stays a static,
        # Snakemake-legal output pattern; only ever requested for valid
        # (bed_id, sample_type, sample) combinations per GROUPS.
        tsv = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/by_sample_preliminary/{sample}_hits.tsv",
    params:
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        stub_dir     = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/by_sample_preliminary",
        # Not a Snakemake input: see the rule docstring above -- used only if
        # it already happens to exist when this rule actually runs.
        cohort_junction_tsv = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        # Empty string when omim_file isn't configured, so the --omim flag
        # is simply omitted from the shell command below (merge_hits.py
        # treats a missing --omim as "skip phenotype/inheritance annotation"
        # rather than requiring it).
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        split_script = workflow.basedir + "/scripts/split_group_hits_by_sample.py",
        merge_script = workflow.basedir + "/scripts/merge_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample}_merge_hits_preliminary.log"
    shell:
        """
        mkdir -p {params.stub_dir} $(dirname {log})
        python -u {params.split_script} \\
            --variant-tsv {input.variant_tsv} \\
            --ase-tsv     {input.ase_tsv} \\
            --tissues     {params.tissues} \\
            --junction-files {input.junction_files} \\
            --cohort-junction-tsv {params.cohort_junction_tsv} \\
            --samples     {wildcards.sample} \\
            --outdir      {params.stub_dir} \\
        2>&1 | tee {log}
        python -u {params.merge_script} \\
            --outfile      {output.tsv} \\
            --sample-name  {wildcards.sample} \\
            --variant-hits {params.stub_dir}/{wildcards.sample}_variant_hits.tsv \\
            --ase-hits     {params.stub_dir}/{wildcards.sample}_ase_hits.tsv \\
            --junction-hits {params.stub_dir}/{wildcards.sample}_junction_hits.tsv \\
            --cohort-junction-hits {params.stub_dir}/{wildcards.sample}_cohort_junction_hits.tsv \\
            {params.omim_flag} \\
        2>&1 | tee -a {log}
        """


# ---------------------------------------------------------------------------
# 6D2. The fully-informed version of the same merge: identical to _6D1 above,
#    except cohort_junction_tsv IS a real input here, so Snakemake waits for
#    cohort_junction_analysis to finish and reruns this rule (and everything
#    downstream: _6E/_6F/_6G) whenever its output changes.
#    Only runs at all when config['merge_hits_include_cohort_junctions'] is
#    True (the default) -- see _group_sample_hits_files() near the top of
#    this file, and _6D1's comment above.
# ---------------------------------------------------------------------------
rule _6D2_merge_sample_hits_with_cohort_junctions:
    input:
        variant_tsv    = _variant_tsv,
        ase_tsv        = _ase_tsv,
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.bed_id, wc.sample_type), t)
            for t in group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type))
        ],
        cohort_junction_tsv = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        tsv = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/by_sample/{sample}_hits.tsv",
    params:
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        stub_dir     = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/by_sample",
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        split_script = workflow.basedir + "/scripts/split_group_hits_by_sample.py",
        merge_script = workflow.basedir + "/scripts/merge_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/{sample}_merge_hits.log"
    shell:
        """
        mkdir -p {params.stub_dir} $(dirname {log})
        python -u {params.split_script} \\
            --variant-tsv {input.variant_tsv} \\
            --ase-tsv     {input.ase_tsv} \\
            --tissues     {params.tissues} \\
            --junction-files {input.junction_files} \\
            --cohort-junction-tsv {input.cohort_junction_tsv} \\
            --samples     {wildcards.sample} \\
            --outdir      {params.stub_dir} \\
        2>&1 | tee {log}
        python -u {params.merge_script} \\
            --outfile      {output.tsv} \\
            --sample-name  {wildcards.sample} \\
            --variant-hits {params.stub_dir}/{wildcards.sample}_variant_hits.tsv \\
            --ase-hits     {params.stub_dir}/{wildcards.sample}_ase_hits.tsv \\
            --junction-hits {params.stub_dir}/{wildcards.sample}_junction_hits.tsv \\
            --cohort-junction-hits {params.stub_dir}/{wildcards.sample}_cohort_junction_hits.tsv \\
            {params.omim_flag} \\
        2>&1 | tee -a {log}
        """


# ---------------------------------------------------------------------------
# 6E. Concatenate every sample's hits.tsv into one group-level all_hits.tsv.
#    Reads from whichever of _6D1/_6D2 _group_sample_hits_files() selects --
#    see that function's docstring above. This is the point where the two
#    branches become mutually exclusive: only the selected one is ever
#    requested, so only it ever runs.
# ---------------------------------------------------------------------------
rule _6E_concat_group_hits:
    input:
        sample_hits = lambda wc: _group_sample_hits_files(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        all_hits = _all_hits_tsv,
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * 2),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/concat_group_hits.log"
    shell:
        """
        mkdir -p $(dirname {output.all_hits}) $(dirname {log})
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.sample_hits} > {output.all_hits} 2> {log}
        echo "Finished concat to {output.all_hits}." >> {log}
        """


# ---------------------------------------------------------------------------
# 6F. Plot candidate hits for a group (calls plot_candidate_hits.py directly).
# ---------------------------------------------------------------------------
rule _6F_plot_group_hits:
    input:
        all_hits = _all_hits_tsv,
    output:
        pathogenic  = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_pathogenic_variant.pdf",
        ase         = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_ASE.pdf",
        junction    = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_outlier_junction.pdf",
        dysreg      = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/genes_with_RNA_dysregulation.pdf",
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
        """
        python -u {params.script} \\
            --infile {input.all_hits} \\
            --outdir {params.outdir} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 6G. Final merge across all sample types sharing a BED panel (unchanged).
# ---------------------------------------------------------------------------
rule _6G_final_merge:
    input:
        all_hits = lambda wc: [(str(group_outdir(gid)) + '/merged_hits/all_hits.tsv') for gid in BED_GROUPS[wc.bed_id]],
    output:
        merged = _cohort_outdir + "/{bed_id}/output/merged_all_hits.tsv",
    threads: lambda wc: _group_threads(wc.bed_id, "final_merge", 1)
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 1024 * _group_mem_gb(wc.bed_id, "final_merge", 2)),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_final_merge.log"
    shell:
        """
        mkdir -p $(dirname {log})
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.all_hits} > {output.merged} 2> {log}
        echo "Finished final merge to {output.merged}." >> {log}
        """
