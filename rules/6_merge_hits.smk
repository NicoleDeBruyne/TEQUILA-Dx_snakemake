"""
rules/6_merge_hits.smk
Cross-sample merge-and-filter stage. Runs once per (bed, sample_type) group.
See docs/rules/6_merge_hits.md for the full stage list, dependency order, and
the {bed_id}/output/sample_types/{sample_type} wildcard convention used throughout
this file.

Per-group merge (builds every sample's ranked hits AND concatenates them
into that group's all_hits.tsv, in one job) is split into two mutually-
exclusive rules:
  _6D1  -- doesn't depend on cohort_junction_analysis (rule 7) at all, so it
          can run as soon as variant/ASE/junction results are ready.
  _6D2 -- the fully-informed version: does depend on rule 7, and picks up
          cohort-comparison junction results once they're available.

config['merge_hits_include_cohort_junctions'] (default True) controls which
one produces all_hits.tsv (and therefore everything chained after it --
_6E/_6F/merged_all_hits.tsv) -- see the `subdir`/output path picked in each
rule below. Only the selected branch is ever requested, so only it ever
runs; the two are never both computed in the same run. Flip the config
value and rerun snakemake to switch which one all_hits.tsv is built from
(this also determines whether that run waits on rule 7 or not).

_6D1/_6D2 build every sample's ranked hit table directly from the group's
already-merged variant/ASE/junction/cohort-junction tables (via
scripts/merge_group_hits.py, using scripts/merge_hits.py as a library) and
write the group's all_hits.tsv in the same job -- one job per group, not
one job per sample plus a separate concat step. (Older versions of this
pipeline split this into a per-sample split+merge rule and a separate
_6E concat rule; that's been folded in here since the per-sample split
served no purpose once _6A/_6B/_6C already merge across the whole group --
see git history / conversation notes if you need the old per-sample
version for reference.)

_6G pools every all_hits.tsv across sample_types (via _6F's
merged_all_hits.tsv) into an UpSet-style category-combination plot -- see
that rule's docstring.
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
# 6D1. Build every sample's ranked candidate-hits table for a group,
#    directly from the group's already-merged variant/ASE/junction tables
#    (scripts/merge_group_hits.py, using scripts/merge_hits.py as a
#    library), and write this group's all_hits.tsv -- one job per group.
#
# PRELIMINARY / fast path: deliberately does NOT depend on
# cohort_junction_analysis (rules/7_cohort_junction_analysis.smk), which can
# take a while. --cohort-junction-tsv is still passed to the script as a
# params (not input) path: if that file already happens to exist on disk
# (e.g. a previous run already completed cohort_junction_analysis for this
# group), it's picked up opportunistically; if not, merge_group_hits.py's
# existing graceful "no cohort data" fallback applies (same as its --omim
# handling).
# Writes to merged_hits/all_hits_preliminary.tsv, NOT the canonical
# merged_hits/all_hits.tsv path _6D2 writes to -- see _6D2 below.
# Only runs at all when config['merge_hits_include_cohort_junctions'] is
# False -- see _6D2's docstring for how that's decided.
# ---------------------------------------------------------------------------
rule _6D1_merge_group_hits_preliminary:
    input:
        variant_tsv    = _variant_tsv,
        ase_tsv        = _ase_tsv,
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.bed_id, wc.sample_type), t)
            for t in group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type))
        ],
    output:
        all_hits = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_hits/all_hits_preliminary.tsv",
    params:
        samples      = lambda wc: GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)],
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        # Not a Snakemake input: see the rule docstring above -- used only if
        # it already happens to exist when this rule actually runs.
        cohort_junction_tsv = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        # Empty string when omim_file isn't configured, so the --omim flag
        # is simply omitted from the shell command below (merge_group_hits.py
        # treats a missing --omim as "skip phenotype/inheritance annotation"
        # rather than requiring it).
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        script       = workflow.basedir + "/scripts/merge_group_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_hits_preliminary.log"
    shell:
        """
        mkdir -p $(dirname {output.all_hits}) $(dirname {log})
        python -u {params.script} \\
            --outfile     {output.all_hits} \\
            --variant-tsv {input.variant_tsv} \\
            --ase-tsv     {input.ase_tsv} \\
            --tissues     {params.tissues} \\
            --junction-files {input.junction_files} \\
            --cohort-junction-tsv {params.cohort_junction_tsv} \\
            --samples     {params.samples} \\
            {params.omim_flag} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 6D2. The fully-informed version of the same group-level merge: identical
#    to _6D1 above, except cohort_junction_tsv IS a real input here, so
#    Snakemake waits for cohort_junction_analysis to finish and reruns this
#    rule (and everything downstream: _6E/_6F) whenever its output changes.
#    Only runs at all when config['merge_hits_include_cohort_junctions'] is
#    True (the default). This produces the canonical all_hits.tsv path
#    (_all_hits_tsv) that _6E/_6F consume; _6D1's all_hits_preliminary.tsv
#    is a dead end otherwise -- flip merge_hits_include_cohort_junctions to
#    False if you want a run built from _6D1 instead.
# ---------------------------------------------------------------------------
rule _6D2_merge_group_hits_with_cohort_junctions:
    input:
        variant_tsv    = _variant_tsv,
        ase_tsv        = _ase_tsv,
        junction_files = lambda wc: [
            _group_junction_final_path(_group_id_from_ids(wc.bed_id, wc.sample_type), t)
            for t in group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type))
        ],
        cohort_junction_tsv = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.bed_id, wc.sample_type)),
    output:
        all_hits = _all_hits_tsv,
    params:
        samples      = lambda wc: GROUPS[_group_id_from_ids(wc.bed_id, wc.sample_type)],
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.bed_id, wc.sample_type)),
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        script       = workflow.basedir + "/scripts/merge_group_hits.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: max(4096, attempt * 4 * 1024),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/merge_group_hits.log"
    shell:
        """
        mkdir -p $(dirname {output.all_hits}) $(dirname {log})
        python -u {params.script} \\
            --outfile     {output.all_hits} \\
            --variant-tsv {input.variant_tsv} \\
            --ase-tsv     {input.ase_tsv} \\
            --tissues     {params.tissues} \\
            --junction-files {input.junction_files} \\
            --cohort-junction-tsv {input.cohort_junction_tsv} \\
            --samples     {params.samples} \\
            {params.omim_flag} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 6E. Plot candidate hits for a group (calls plot_candidate_hits.py directly).
#    (Formerly _6F -- shifted down a letter when the old per-sample _6D1/
#    _6D2 + separate _6E concat rule were folded into single group-level
#    _6D1/_6D2 rules above, to keep letters matching execution order.)
# ---------------------------------------------------------------------------
rule _6E_plot_group_hits:
    input:
        all_hits = _all_hits_tsv,
    output:
        # plot_candidate_hits.py's plot() writes two files per category
        # (a stacked barplot and a boxplot-with-dots), never a single file
        # at the bare category name -- these must match its actual
        # filenames exactly or Snakemake reports the rule as failed even
        # when the script ran fine (each is f"{prefix}_barplot{ext}" /
        # f"{prefix}_boxplot{ext}" off the *args* passed below, i.e. off
        # the bare "genes_with_X.pdf" names, not off these output paths).
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
        """
        python -u {params.script} \\
            --infile {input.all_hits} \\
            --outdir {params.outdir} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 6F. Final merge across all sample types sharing a BED panel (unchanged).
#    (Formerly _6G -- shifted down a letter, see _6E's comment above.)
# ---------------------------------------------------------------------------
rule _6F_final_merge:
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


# ---------------------------------------------------------------------------
# 6G. UpSet-style boxplot of merged_all_hits.tsv's four hit categories
#    (variant, ASE, outlier_junction, cohort_outlier_junction): for every
#    non-empty combination of those categories, a boxplot of "how many
#    genes did this sample have in exactly this combination", one dot per
#    sample colored by sample_type, plus the standard UpSet
#    combination-membership matrix underneath. Pools across every
#    sample_type sharing this BED panel, same scope as merged_all_hits.tsv
#    itself (bed_samples(), not group-scoped -- see rules/9_plot_cohort_info.smk
#    for the same pooling pattern).
# ---------------------------------------------------------------------------
rule _6G_plot_hits_upset:
    input:
        all_hits = _cohort_outdir + "/{bed_id}/output/merged_all_hits.tsv",
    output:
        pdf = _cohort_outdir + "/{bed_id}/output/hits_upset_plot.pdf",
        tsv = _cohort_outdir + "/{bed_id}/output/hits_upset_counts.tsv",
    params:
        samples      = lambda wc: bed_samples(wc.bed_id),
        sample_types = lambda wc: [SAMPLES[s]["sample_type"] for s in bed_samples(wc.bed_id)],
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
        """
        mkdir -p {params.outdir} $(dirname {log})
        python -u {params.script} \\
            --infile      {input.all_hits} \\
            --samples     {params.samples} \\
            --sample-types {params.sample_types} \\
            --outdir      {params.outdir} \\
            --title       "{params.title}" \\
        2>&1 | tee {log}
        """
