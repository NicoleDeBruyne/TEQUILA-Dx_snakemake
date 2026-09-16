"""
rules/9_merge_results.smk
Cross-sample merge stage: everything that combines per-sample (rules 1-7)
and per-group cohort-junction (rule 8) results into cohort-level tables and
plots. This is the only stage that touches multiple samples' results at
once -- no rule in this file reads a BAM directly, which is what makes it
cheap to rerun when cohort membership (or an alias) changes without
rerunning any of the expensive per-sample BAM scans in rules 6/7.

Three groups of rules, in the order they appear below:

  A. Hit ranking (_9A-_9F, unchanged from this pipeline's original
     rules/6_merge_hits.smk): merges variant/ASE/GTEx-junction results per
     (bed, sample_type) group, ranks every gene into a tiered candidate-hit
     table (scripts/merge_hits.py's tier logic, see conversation notes for
     the full decision tree), and produces the final cross-group merge +
     upset plot. See the docstring immediately above _9A below for the
     _9D1/_9D2 preliminary-vs-cohort-junction-informed split.

  B. QC merge (_9G-_9J): combines rule 6's per-sample QC TSVs
     (on_target/read_attributes/full_length_ratio) into cohort-wide
     matrices and plots. _9I1/_9I2 (sample-type validation) are unchanged
     from the old rules/8_cohort_qc.smk -- they were already cohort-level,
     reading step 5's per-sample junction matrices directly rather than
     any per-sample BAM data.

  C. Gene-quantification merge (_9K-_9N5): combines rule 7's per-sample
     gene-quantification TSVs (count/coverage/assignment) and StringTie
     GTFs into cohort-wide matrices, boxplots, and the AMALGAM
     sub-pipeline's remaining cohort-level steps (transcript-discovery
     merge through final aggregation -- see _9N1's docstring for why only
     StringTie itself, not the rest of AMALGAM, could move to rule 7).
"""

from math import ceil

# NOTE: output:/input:/log: path templates below use string concatenation,
# not f-strings, to combine a config value with a literal Snakemake
# wildcard placeholder like "{bed_id}" -- an f-string's "{{bed_id}}" escape
# (to produce a literal "{bed_id}") does not survive Snakemake's own rule
# parsing and raises a NameError at load time. Paths reused across more than
# one rule (e.g. all_candidate_variants.tsv, an output of _9A and an input
# of _9D1) are also factored out here so both rules stay in sync.
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


# ===========================================================================
# A. Hit ranking (_9A-_9F)
# ===========================================================================

# ---------------------------------------------------------------------------
# 9A. Merge & filter variant calls across all samples in a group
# ---------------------------------------------------------------------------
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
# 9B. Merge & filter ASE results across all samples in a group
# ---------------------------------------------------------------------------
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
# 9C. Merge & filter outlier junctions, once per (group, tissue)
# ---------------------------------------------------------------------------
rule _9C_merge_group_junctions:
    input:
        junction_files = lambda wc: _group_tissue_junction_files(
            _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), wc.tissue),
    output:
        # Static, wildcard-only path -- the dynamically-named file the
        # script actually produces gets cp'd here at the end of the shell block.
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
# 9D1. Build every sample's ranked candidate-hits table for a group,
#    directly from the group's already-merged variant/ASE/junction tables
#    (scripts/merge_group_hits.py, using scripts/merge_hits.py as a
#    library), and write this group's all_hits.tsv -- one job per group.
#
# PRELIMINARY / fast path: deliberately does NOT depend on
# cohort_junction_analysis (rules/8_cohort_junction_analysis.smk), which can
# take a while. --cohort-junction-tsv is still passed to the script as a
# params (not input) path: if that file already happens to exist on disk
# (e.g. a previous run already completed cohort_junction_analysis for this
# group), it's picked up opportunistically; if not, merge_group_hits.py's
# existing graceful "no cohort data" fallback applies (same as its --omim
# handling).
# Writes to merged_hits/all_hits_preliminary.tsv, NOT the canonical
# merged_hits/all_hits.tsv path _9D2 writes to -- see _9D2 below.
# Only runs at all when config['merge_hits_include_cohort_junctions'] is
# False -- see _9D2's docstring for how that's decided.
# ---------------------------------------------------------------------------
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
        # Not a Snakemake input: see the rule docstring above -- used only if
        # it already happens to exist when this rule actually runs.
        cohort_junction_tsv = lambda wc: _cja_outliers_filtered_path(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
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
# 9D2. The fully-informed version of the same group-level merge: identical
#    to _9D1 above, except cohort_junction_tsv IS a real input here, so
#    Snakemake waits for cohort_junction_analysis (rule 8) to finish and
#    reruns this rule (and everything downstream: _9E/_9F) whenever its
#    output changes. Only runs at all when
#    config['merge_hits_include_cohort_junctions'] is True (the default).
#    This produces the canonical all_hits.tsv path (_all_hits_tsv) that
#    _9E/_9F consume; _9D1's all_hits_preliminary.tsv is a dead end
#    otherwise -- flip merge_hits_include_cohort_junctions to False if you
#    want a run built from _9D1 instead.
#
#    gene_expression_matrix comes from this same file's _9M rule below
#    (splice-site-assignment gene quantification) -- both live in
#    rules/9_merge_results.smk now, so there's no cross-file dependency to
#    track, just an ordinary Snakemake input.
# ---------------------------------------------------------------------------
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
            + "/output/gene_quantification/by_assignment/gene_assignment_matrix.tsv"
        ) if config.get("quantify_genes") else [],
    output:
        all_hits = _all_hits_tsv,
    params:
        samples      = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        tissues      = lambda wc: group_tissues(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
        omim_flag    = ("--omim " + config["omim_file"]) if config.get("omim_file") else "",
        gene_expression_flag = lambda wc, input: (
            "--gene-expression-matrix " + str(input.gene_expression_matrix)
        ) if config.get("quantify_genes") else "",
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
            {params.gene_expression_flag} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# 9E. Plot candidate hits for a group (calls plot_candidate_hits.py directly).
# ---------------------------------------------------------------------------
rule _9E_plot_group_hits:
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
# 9F. Final merge across all sample types sharing a BED panel, plus a
#    simplified companion file -- core columns only, plus every *_jxns
#    column merged and deduplicated into one "jxns" column -- written
#    right after, in the same rule. See scripts/simplify_all_hits.py's
#    module docstring for the simplification rules.
# ---------------------------------------------------------------------------
rule _9F_final_merge:
    input:
        all_hits = lambda wc: [(str(group_outdir(gid)) + '/merged_hits/all_hits.tsv') for gid in BED_GROUPS[(wc.cohort_id, wc.bed_id)]],
    output:
        merged          = _cohort_outdir + "/{bed_id}/output/merged_all_hits.tsv",
        simplified      = _cohort_outdir + "/{bed_id}/output/merged_all_hits_simplified.tsv",
        # Always produced alongside `simplified` (mirrors it verbatim when no
        # sample on this bed panel has an alias configured); only actually
        # *requested* by rule all via all_outputs() when bed_has_alias() is true.
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
        """
        mkdir -p $(dirname {log})
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.all_hits} > {output.merged} 2> {log}
        echo "Finished final merge to {output.merged}." >> {log}

        python -u {params.script} \\
            --infile  {output.merged} \\
            --outfile {output.simplified} \\
            --alias-outfile {output.simplified_alias} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee -a {log}
        """


# ---------------------------------------------------------------------------
# 9F. UpSet-style boxplot of merged_all_hits.tsv's four hit categories
#    (variant, ASE, outlier_junction, cohort_outlier_junction): for every
#    non-empty combination of those categories, a boxplot of "how many
#    genes did this sample have in exactly this combination", one dot per
#    sample colored by sample_type, plus the standard UpSet
#    combination-membership matrix underneath. Pools across every
#    sample_type sharing this BED panel, same scope as merged_all_hits.tsv
#    itself (bed_samples(), not group-scoped).
# ---------------------------------------------------------------------------
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


# ===========================================================================
# B. QC merge (_9G-_9J) -- combines rule 6's per-sample QC TSVs
# ===========================================================================

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
        # Always produced alongside `ontarget_pdf` (mirrors it verbatim when
        # no sample on this bed panel has an alias configured); only
        # actually *requested* by rule all via all_outputs() when
        # bed_has_alias() is true.
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
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --groups    {params.groups} \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _9H_merge_read_attributes:
    input:
        infiles = lambda wc: _bed_qc_files(wc.cohort_id, wc.bed_id, "read_attributes"),
    output:
        tsv           = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_attributes.tsv",
        length_boxplot = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_lengths.pdf",
        # Always produced alongside `length_boxplot` (mirrors it verbatim
        # when no sample on this bed panel has an alias configured); only
        # actually *requested* by rule all via all_outputs() when
        # bed_has_alias() is true.
        length_boxplot_alias = _cohort_outdir + "/{bed_id}/output/cohort_qc/read_attributes/{bed_id}_read_lengths_alias.pdf",
    params:
        outprefix  = lambda wc: (str(bed_outdir(wc.cohort_id, wc.bed_id)) + '/cohort_qc/read_attributes/' + str(wc.bed_id) + '_read_attributes'),
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
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


def _group_junction_matrix_inputs(group_id):
    return [(str(SAMPLES[s]['outdir']) + '/junction_analysis/junction_counts/' + str(s) + '_junction_count_matrix.tsv')
            for s in GROUPS[group_id]]


rule _9I1_build_group_junction_matrix:
    # Unchanged from the old rules/8_cohort_qc.smk: reads step 5's
    # per-sample junction count matrices directly, not any rule 6 output --
    # this analysis was already cohort-level with no BAM access of its own.
    input:
        matrices = lambda wc: _group_junction_matrix_inputs(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)),
    output:
        matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/merged_junction_analysis/junction_count_matrix.tsv",
    params:
        samples = lambda wc: GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)],
        script  = workflow.basedir + "/scripts/build_group_junction_matrix.py",
    threads: lambda wc: _group_threads(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "build_group_junction_matrix", config["threads"])
    resources:
        # Scales with group size (the matrix has one column per sample).
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(8, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/group_junction_matrix.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --infiles      {input.matrices} \\
            --sample-names {params.samples} \\
            --outfile      {output.matrix} \\
            --threads      {threads} \\
        2>&1 | tee {log}
        """


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
        # Always produced alongside `heatmap` (mirrors it verbatim when no
        # sample on this bed panel has an alias configured); only actually
        # *requested* by rule all via all_outputs() when bed_has_alias() is true.
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
        # Pools every sample_type sharing this BED panel; the reference GTEx
        # matrices dominate baseline memory use, hence the high floor -- 256GB
        # covers that regardless of cohort size, and scales up further only
        # once the local cohort itself exceeds that many samples.
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(256, len(bed_samples(wc.cohort_id, wc.bed_id))),
        runtime = 1440,   # 1 day, matches the old bash wrapper's --time=1-00:00:00
    log:
        _cohort_outdir + "/{bed_id}/logs/{bed_id}_validate_sample_types.log"
    shell:
        """
        mkdir -p $(dirname {log})
        python -u {params.script} \\
            --matrix-refs  {input.gtex_matrices} \\
            --ref-names    {params.ref_names} \\
            --ref-colors   {params.ref_colors} \\
            --matrix-query {input.query_matrices} \\
            --query-names  {params.query_names} \\
            --query-colors {params.query_colors} \\
            --bed          {input.bed} \\
            --outprefix    {params.outprefix} \\
            --alias-map    {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _9J_merge_full_length_ratio:
    input:
        infiles = lambda wc: _bed_qc_files(wc.cohort_id, wc.bed_id, "full_length_ratio"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_matrix.tsv",
        read_counts  = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_read_counts_matrix.tsv",
        heatmap      = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_heatmap.pdf",
        # Always produced alongside `matrix` (mirrors it verbatim when no
        # sample on this bed panel has an alias configured); only actually
        # *requested* by rule all via all_outputs() when bed_has_alias() is true.
        matrix_alias = _cohort_outdir + "/{bed_id}/output/cohort_qc/full_length_ratio/{bed_id}_full_length_ratio_matrix_alias.tsv",
    params:
        sample_types = lambda wc: _quoted([SAMPLES[s]["sample_type"] for s in bed_samples(wc.cohort_id, wc.bed_id)]),
        # Same colors _9I2_validate_sample_types uses -- resolved here via
        # the Snakefile's own sample_type_color() and passed through as a
        # plain column, rather than the script inventing its own palette.
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
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles      {input.infiles} \\
            --sample-types {params.sample_types} \\
            --colors       {params.colors} \\
            --outprefix    {params.outprefix} \\
            --title        {params.title:q} \\
            --min-reads    {params.min_reads} \\
            --alias-map    {params.alias_args} \\
        2>&1 | tee {log}
        """


# ===========================================================================
# C. Gene-quantification merge (_9K-_9N5) -- combines rule 7's per-sample
#    gene-quantification TSVs and StringTie GTFs
# ===========================================================================

def _sample_quant_file(sample, name):
    return str(SAMPLES[sample]['outdir']) + '/gene_quantification/' + str(sample) + '_' + str(name) + '.tsv'

def _group_quant_files(group_id, name):
    return [_sample_quant_file(s, name) for s in GROUPS[group_id]]


rule _9K_merge_gene_count:
    input:
        infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_count"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_raw.tsv",
        # Always produced alongside `matrix` (mirrors it verbatim when no
        # sample in this group has an alias configured); only actually
        # *requested* by rule all via all_outputs() when group_has_alias() is true.
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_count/gene_count_matrix_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_count/gene_count"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (read count)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/quantify_gene_expression.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_count_quantification.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --metric    count \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _9L_merge_gene_coverage:
    input:
        infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_coverage"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_raw.tsv",
        # Always produced alongside `matrix` (mirrors it verbatim when no
        # sample in this group has an alias configured); only actually
        # *requested* by rule all via all_outputs() when group_has_alias() is true.
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_coverage/gene_coverage_matrix_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_coverage/gene_coverage"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (max coverage)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/quantify_gene_expression.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_coverage_quantification.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.infiles} \\
            --metric    coverage \\
            --outprefix {params.outprefix} \\
            --title     {params.title:q} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """


rule _9M_merge_gene_by_assignment:
    input:
        infiles       = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "gene_assignment"),
        stats_infiles = lambda wc: _group_quant_files(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type), "read_outcomes"),
    output:
        matrix       = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix.tsv",
        matrix_raw   = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw.tsv",
        matrix_raw_all_genes = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_raw_all_genes.tsv",
        assignment_stats     = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_read_outcomes.tsv",
        assignment_stats_pdf = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_read_outcomes.pdf",
        # Always produced alongside `matrix` (mirrors it verbatim when no
        # sample in this group has an alias configured); only actually
        # *requested* by rule all via all_outputs() when group_has_alias() is true.
        matrix_alias = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_assignment/gene_assignment_matrix_alias.tsv",
    params:
        outprefix  = lambda wc: (config["output_dir"] + "/" + str(wc.cohort_id) + "/" + str(wc.bed_id) + "/output/sample_types/" + str(wc.sample_type)
                                 + "/output/gene_quantification/by_assignment/gene_assignment"),
        title      = lambda wc: config.get("gene_quant_title") or (str(wc.cohort_id) + " " + str(wc.bed_id) + " " + str(wc.sample_type) + " relative gene expression (splice-site assignment)"),
        alias_args = lambda wc: _quoted(alias_map_args(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        script     = workflow.basedir + "/scripts/quantify_gene_by_assignment.py",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(4, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]) // 8),
        runtime = 60,
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/gene_assignment_quantification.log"
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles       {input.infiles} \\
            --stats-infiles {input.stats_infiles} \\
            --outprefix     {params.outprefix} \\
            --title         {params.title:q} \\
            --alias-map     {params.alias_args} \\
        2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# AMALGAM sub-pipeline helpers (_9N1-_9N5). AMALGAM's own tool + assets +
# dedicated conda env live under config["amalgam_dir"] (see setup.sh's
# AMALGAM section) -- referenced as-is throughout, same convention already
# used for config["annovar_dir"] elsewhere in this pipeline.
# ---------------------------------------------------------------------------
def _amalgam_group_dir(cohort_id, bed_id, sample_type):
    return config["output_dir"] + "/" + str(cohort_id) + "/" + str(bed_id) + "/output/sample_types/" + str(sample_type) + "/output/gene_quantification/by_amalgam"


def _amalgam_stringtie_gtf(sample):
    """Lives under the SAMPLE's own outdir now (written by rule 7's _7D),
    not the group's by_amalgam/stringtie/ directory -- this is the one
    piece of the AMALGAM sub-pipeline that's genuinely per-sample with no
    cohort dependency, so it moved to rules/7_sample_gene_quantification.smk."""
    return str(SAMPLES[sample]['outdir']) + "/gene_quantification/amalgam/" + str(sample) + "_stringtie.gtf"


def _amalgam_quantification_tsv(cohort_id, bed_id, sample_type, sample):
    return _amalgam_group_dir(cohort_id, bed_id, sample_type) + "/quantification/" + sample + "_transcript_quantification.tsv"


rule _9N1_amalgam_merge_gtfs:
    # Step 2 of AMALGAM's own pipeline (see its README): GffCompare merges
    # every sample's StringTie GTF in this group (from rule 7's _7D) with
    # the reference annotation into one combined transcript set. The
    # reference annotation MUST be first in the input list (see AMALGAM's
    # README) -- gtf_list.tsv is built in that order below.
    #
    # This step -- discovering the cohort's shared transcript set from
    # every sample's own de-novo assembly -- is inherently cohort-level:
    # adding one sample to the group can change which transcripts get
    # discovered for everyone, so this (and every AMALGAM step after it)
    # always reruns on a cohort membership change, unlike _7D itself.
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


rule _9N2_amalgam_build_transcriptome:
    # Step 3: Build_Transcriptome.py identifies high-confidence,
    # full-length transcripts from the GffCompare merge, using AMALGAM's
    # bundled RefTSS/PolyASite reference BED files (config["amalgam_dir"]/
    # assets/, cloned alongside the tool itself -- see setup.sh). Keeps
    # BOTH the uncompressed filtered.gtf (needed as-is by _9N3 below,
    # matching AMALGAM's own README) and the sorted/bgzip/tabix-indexed
    # filtered.gtf.gz (needed by _9N4) as real tracked outputs.
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
        main_env       = workflow.basedir + "/envs/conda_env",  # bgzip/tabix (htslib) -- not part of AMALGAM's own env
        amalgam_dir    = config["amalgam_dir"],
        merge_prefix   = lambda wc: _amalgam_group_dir(wc.cohort_id, wc.bed_id, wc.sample_type) + "/annotation/merged",
    threads: 1
    resources:
        mem_mb  = lambda wc, attempt: attempt * 1024 * max(32, len(GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)])),
        runtime = config["time"],
    log:
        _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/logs/amalgam_build_transcriptome.log"
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


rule _9N3_amalgam_annotate_orf:
    # Step 4 (OPTIONAL per AMALGAM's own README): Annotate_ORF.py adds
    # open-reading-frame annotations to the filtered transcriptome.
    # NOTE: this step's output is currently NOT consumed by anything
    # downstream -- _9N4 (transcript quantification) runs against
    # _9N2's filtered.gtf.gz, not this rule's annotated.gtf.gz, matching
    # exactly how the group's own AMALGAM submission script was written
    # (Quantify_Transcripts.py -g pointed at step3's output, not step4's).
    # Kept as a real rule (rather than dropped) since the group's script
    # ran it unconditionally, but flagging here in case that was meant to
    # feed step 5 and didn't due to an oversight in the original script.
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


rule _9N4_amalgam_quantify_transcripts:
    # Step 5: Quantify_Transcripts.py, per sample, against the GROUP's
    # filtered transcriptome from _9N2 (every sample in a group shares the
    # same transcriptome; only the BAM being quantified differs). Unlike
    # rule 7's _7D (plain StringTie, no cohort dependency), this genuinely
    # needs the group's cohort-built transcriptome first, so it stays here
    # rather than moving to rule 7 -- bed_id/sample_type are real wildcards
    # for that reason, and output lives under the group's
    # by_amalgam/quantification/ directory, not the sample's own.
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


rule _9N5_amalgam_aggregate_matrices:
    # Step 6 (not part of AMALGAM itself -- the group's own aggregation
    # step from its submission script): combine every sample's
    # transcript-level quantification in this group into cohort-wide
    # transcript and gene matrices, via scripts/aggregate_amalgam_matrices.py.
    input:
        tsvs = lambda wc: [
            _amalgam_quantification_tsv(wc.cohort_id, wc.bed_id, wc.sample_type, s)
            for s in GROUPS[_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)]
        ],
    output:
        transcript_matrix = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_transcript_matrix.tsv",
        gene_matrix        = _cohort_outdir + "/{bed_id}/output/sample_types/{sample_type}/output/gene_quantification/by_amalgam/quantification/gene_amalgam_gene_matrix.tsv",
        # Always produced alongside `gene_matrix` (mirrors it verbatim when
        # no sample in this group has an alias configured); only actually
        # *requested* by rule all via all_outputs() when group_has_alias() is true.
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
        """
        mkdir -p $(dirname {output.gene_matrix}) $(dirname {log})
        python -u {params.script} \\
            --infiles   {input.tsvs} \\
            --samples   {params.samples} \\
            --outprefix {params.outprefix} \\
            --alias-map {params.alias_args} \\
        2>&1 | tee {log}
        """
