
import yaml
import os
import re
from pathlib import Path
from collections import defaultdict
from math import ceil
import shlex

shell.prefix("set -euo pipefail;")


# Load the run config
with open(config["run"]) as fh:
    run_cfg = yaml.safe_load(fh)

# Load sample names
SAMPLES = run_cfg["samples"]

# Set output directory keys
for _key, _val in run_cfg.items():
    if _key in ("samples", "output_dir"):
        continue
    config[_key] = _val

# Set top-level output directory
if not config.get("output_dir"):
    config["output_dir"] = run_cfg.get("output_dir", "")
if not config.get("output_dir"):
    raise ValueError(
        "output_dir must be set, either via --config output_dir=<path> or as a "
        "top-level 'output_dir:' key in the run config YAML (alongside 'samples:'), "
        "for the pipeline to know where to write its output layout."
    )

# Default each sample's output directory if it wasn't set explicitly.
for _sample, _entry in SAMPLES.items():
    if not _entry.get("outdir"):
        _entry["outdir"] = config["output_dir"] + "/samples/" + _sample + "/output"

# Config keys that point at reference/resource files bundled with the pipeline
_BUNDLED_PATH_KEYS = [
    "genome", "annotation",
    "conda_env", "conda_env_compile_variants", "gnomad_base", "clinvar_vcf", "annovar_dir",
    "amalgam_dir",
    "cadd_data_dir", "cadd_script", "cadd_local_prescored_snv", "cadd_local_prescored_indel",
    "gnomad_mito_vcf", "longcallr_bin",
    "nanots_model_unphased", "nanots_model_phased", "gtex_data_dir", "omim_file",
    "nanots_image", "clair3_rna_image", "deepvariant_image",
    "spliceai_prescored_snv_vcf", "spliceai_prescored_indel_vcf", "hpo_obo", "softclip_alu_annotation",
]
_URL_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*://")

# Set remote sentinels for gnomAD, ClinVar, etc. if necessary
def _is_remote_sentinel(value):
    return isinstance(value, str) and value.strip().lower() == "remote"
for _key in _BUNDLED_PATH_KEYS:
    _val = config.get(_key)
    if _val and not os.path.isabs(_val) and not _URL_RE.match(_val) and not _is_remote_sentinel(_val):
        config[_key] = os.path.join(workflow.basedir, _val)

# If a resource's config value is "remote", download it from these public URLs instead of using a local/bundled path
_REMOTE_GNOMAD_BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes"
_REMOTE_GNOMAD_MITO_VCF = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/"
                           "vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz")
_REMOTE_CLINVAR_VCF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
_REMOTE_CADD_PRESCORED_URL = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"

def _resolved_gnomad_base():
    return _REMOTE_GNOMAD_BASE if _is_remote_sentinel(config["gnomad_base"]) else config["gnomad_base"]

def _resolved_gnomad_mito_vcf():
    return _REMOTE_GNOMAD_MITO_VCF if _is_remote_sentinel(config["gnomad_base"]) else config["gnomad_mito_vcf"]

def _resolved_clinvar_vcf():
    return _REMOTE_CLINVAR_VCF if _is_remote_sentinel(config["clinvar_vcf"]) else config["clinvar_vcf"]

def _cadd_use_local():
    return not _is_remote_sentinel(config["cadd_script"])

# --- Sample/tissue/cohort/group bookkeeping ---------------------------------
# Samples are organized into cohorts (config-defined, plus an implicit "cohort_all"),
# and within a cohort into "groups" keyed by (cohort, bed file, sample_type) for group-level merge/analysis rules.

def _parse_tissues(raw):
    if isinstance(raw, list):
        return [str(t).strip() for t in raw]
    return [t.strip() for t in str(raw).strip("[]").split(",")]

def all_tissues():
    tissues = set()
    for s in SAMPLES.values():
        tissues.update(_parse_tissues(s["tissues"]))
    return sorted(tissues)

TISSUES = all_tissues()

def sample_tissues(sample):
    return _parse_tissues(SAMPLES[sample]["tissues"])

COHORTS = {"cohort_all": list(SAMPLES.keys())}
for _cid, _members in config.get("cohorts", {}).items():
    if _cid == "cohort_all":
        raise ValueError("'cohort_all' is a reserved cohort name (it always contains every sample "
                         "in the run config) and can't be redefined under config['cohorts'].")
    _unknown = [m for m in _members if m not in SAMPLES]
    if _unknown:
        raise ValueError("cohorts['" + str(_cid) + "'] references sample(s) not present in this run "
                         "config's 'samples:' block: " + str(_unknown))
    COHORTS[_cid] = list(_members)

def _bed_id(bed):
    return Path(bed).stem

def all_bed_files():
    return sorted(set(SAMPLES[s]["bed"] for s in SAMPLES))

def _group_id(cohort_id, bed, sample_type):
    return (str(cohort_id) + '_' + str(_bed_id(bed)) + '_' + str(sample_type))

def _group_id_from_ids(cohort_id, bed_id, sample_type):
    return (str(cohort_id) + '_' + str(bed_id) + '_' + str(sample_type))

def all_groups():
    groups = defaultdict(list)
    for cid, members in COHORTS.items():
        for s in members:
            if "sample_type" not in SAMPLES[s]:
                raise ValueError(
                    ("Sample '" + str(s) + "' is missing a 'sample_type' field in the run config, ")
                    + 'required for grouping samples during the cohort-level merge stages.'
                )
            bed = SAMPLES[s]["bed"]
            sample_type = SAMPLES[s]["sample_type"]
            gid = _group_id(cid, bed, sample_type)
            groups[gid].append(s)
            GROUP_COHORT_ID[gid] = cid
            GROUP_BED_ID[gid] = _bed_id(bed)
            GROUP_SAMPLE_TYPE[gid] = sample_type
    return dict(groups)

GROUP_COHORT_ID = {}
GROUP_BED_ID = {}
GROUP_SAMPLE_TYPE = {}
GROUPS = all_groups()

def group_tissues(group_id):
    tissues = set()
    for s in GROUPS[group_id]:
        tissues.update(sample_tissues(s))
    return sorted(tissues)

def group_outdir(group_id):
    return (str(config['output_dir']) + '/' + str(GROUP_COHORT_ID[group_id]) + '/' + str(GROUP_BED_ID[group_id]) + '/output/sample_types/' + str(GROUP_SAMPLE_TYPE[group_id]) + '/output')

# --- Cohort junction analysis settings ---------------
# 8A always runs to compute raw metrics. 8B runs only when n >= cohort_jxn_zscore_n_threshold (default 10).
# Small groups (10 <= n < 30) use modified z-score; larger groups (n >= 30) use beta-binomial
# (auto-selected by sample count unless overridden by cohort_jxn_method).

def _cja_method_for_group(group_id):
    method = config.get("cohort_jxn_method", "auto")
    if method in ("beta_binomial", "modified_zscore"):
        return method
    if method != "auto":
        raise ValueError(
            "config['cohort_jxn_method'] must be 'auto', 'beta_binomial', or "
            "'modified_zscore', got " + repr(method)
        )
    n     = len(GROUPS[group_id])
    min_n = config.get("cohort_jxn_beta_min_samples", 30)
    return "beta_binomial" if n >= min_n else "modified_zscore"

def _cja_thr_label(group_id):
    if _cja_method_for_group(group_id) == "beta_binomial":
        return "padj" + str(config["cohort_jxn_beta_padj_threshold"]) + "_delta" + str(config["cohort_jxn_beta_delta_threshold"])
    z = config.get("cohort_jxn_z_threshold", 3.5)
    d = config.get("cohort_jxn_z_delta_threshold", 0.1)
    return "z" + str(z) + "_delta" + str(d)

def _cja_outliers_filtered_path(group_id):
    thr_label = _cja_thr_label(group_id)
    return (group_outdir(group_id) + "/cohort_junction_analysis/" + group_id + "_" + thr_label
            + "/" + group_id + "_outliers_filtered.tsv")

def _cja_thr_flag(group_id):
    if _cja_method_for_group(group_id) == "beta_binomial":
        return "--bb-thresholds " + str(config["cohort_jxn_beta_padj_threshold"]) + ":" + str(config["cohort_jxn_beta_delta_threshold"])
    z = config.get("cohort_jxn_z_threshold", 3.5)
    d = config.get("cohort_jxn_z_delta_threshold", 0.1)
    return "--z-thresholds " + str(z) + ":" + str(d)

def cja_enabled(group_id):
    return flag("cohort_junction_analysis") and len(GROUPS[group_id]) >= config.get("cohort_jxn_zscore_n_threshold", 10)

def _cja_n_threshold(group_id):
    if _cja_method_for_group(group_id) == "beta_binomial":
        return config.get("cohort_jxn_beta_n_threshold", 30)
    return config.get("cohort_jxn_zscore_n_threshold", 10)

# --- Cohort transcript analysis settings ---------------
# Runs only for groups with at least cohort_jxn_beta_min_samples (default 30) samples, using beta-binomial/Dirichlet-multinomial tests.

def _ctx_dir(cohort_id, bed_id, sample_type):
    return config["output_dir"] + "/" + str(cohort_id) + "/" + str(bed_id) + "/output/sample_types/" + str(sample_type) + "/output/cohort_transcript_analysis"


def _ctx_method(wc):
    return _cja_method_for_group(_group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type))


def _ctx_thr_label(group_id):
    tvd = str(config.get("cohort_tx_tvd_threshold", 0.1))
    if _cja_method_for_group(group_id) == "beta_binomial":
        return "padj" + str(config["cohort_jxn_beta_padj_threshold"]) + "_delta" + str(config["cohort_jxn_beta_delta_threshold"]) + "_tvd" + tvd
    return "z" + str(config.get("cohort_jxn_z_threshold", 3.5)) + "_delta" + str(config.get("cohort_jxn_z_delta_threshold", 0.1)) + "_tvd" + tvd


def _ctx_filter_args(wc):
    gid = _group_id_from_ids(wc.cohort_id, wc.bed_id, wc.sample_type)
    tvd = " --tvd " + str(config.get("cohort_tx_tvd_threshold", 0.1))
    if _cja_method_for_group(gid) == "beta_binomial":
        return ("--padj " + str(config["cohort_jxn_beta_padj_threshold"])
                + " --delta " + str(config["cohort_jxn_beta_delta_threshold"]) + tvd)
    return ("--z " + str(config.get("cohort_jxn_z_threshold", 3.5))
            + " --delta " + str(config.get("cohort_jxn_z_delta_threshold", 0.1)) + tvd)


def ctx_enabled(group_id):
    return (flag("gene_quantification") and flag("cohort_transcript_analysis")
            and len(GROUPS[group_id]) >= config.get("cohort_jxn_beta_min_samples", 30))


def ctx_files(group_id):
    d = _ctx_dir(GROUP_COHORT_ID[group_id], GROUP_BED_ID[group_id], GROUP_SAMPLE_TYPE[group_id])
    st = str(GROUP_SAMPLE_TYPE[group_id])
    thr = _ctx_thr_label(group_id)
    return dict(
        usage_matrix=d + "/" + st + "_transcript_usage_matrix.tsv",
        gene_status=d + "/" + st + "_transcript_usage_gene_status.tsv",
        isoform_outliers=d + "/" + st + "_transcript_usage_outliers_" + thr + ".tsv",
        landscape_outliers=d + "/" + st + "_gene_isoform_landscape_outliers_" + thr + ".tsv",
        gtf=group_outdir(group_id) + "/gene_quantification/by_amalgam/annotation/filtered.gtf.gz",
        gtf_tbi=group_outdir(group_id) + "/gene_quantification/by_amalgam/annotation/filtered.gtf.gz.tbi",
    )


def ctx_outputs(group_id):
    f = ctx_files(group_id)
    return [f["isoform_outliers"], f["landscape_outliers"]]


def ctx_unavailable_reason(group_id):
    if not flag("gene_quantification"):
        return "gene quantification (AMALGAM) is turned off"
    if not flag("cohort_transcript_analysis"):
        return "the transcript analysis is turned off"
    n_min = config.get("cohort_jxn_beta_min_samples", 30)
    return ("the transcript analysis needs at least " + str(n_min) + " " + str(GROUP_SAMPLE_TYPE[group_id])
            + " samples (this cohort has " + str(len(GROUPS[group_id])) + ")")


def softclip_outputs(group_id):
    d = str(group_outdir(group_id)) + "/softclip_analysis/" + str(GROUP_SAMPLE_TYPE[group_id])
    return [d + "_softclip_outliers.tsv", d + "_softclip_clusters_all.tsv", d + "_softclip_outliers.pdf"]


def all_bed_groups():
    groups = defaultdict(list)
    for gid, members in GROUPS.items():
        bed = SAMPLES[members[0]]["bed"]
        groups[(GROUP_COHORT_ID[gid], _bed_id(bed))].append(gid)
    return dict(groups)

BED_GROUPS = all_bed_groups()

def bed_path(cohort_id, bed_id):
    gid = BED_GROUPS[(cohort_id, bed_id)][0]
    return SAMPLES[GROUPS[gid][0]]["bed"]

def bed_samples(cohort_id, bed_id):
    return [s for gid in BED_GROUPS[(cohort_id, bed_id)] for s in GROUPS[gid]]

def bed_outdir(cohort_id, bed_id):
    return (str(config['output_dir']) + '/' + str(cohort_id) + '/' + str(bed_id) + '/output')

def _quoted(items):
    return [shlex.quote(str(x)) for x in items]

# Nature palette in the order colours are handed out (same list as scripts/palette.py SELECT_ORDER)
_SELECT_ORDER = [
    "#c23637", "#009099", "#9d4588", "#c69528", "#3d892e", "#ea6302", "#0068a9",
    "#d95d5b", "#42b4b5", "#b271ab", "#e8c048", "#57aa3e", "#f48f3e", "#4c8fca",
    "#ea9a9c", "#8cccce", "#cca0ca", "#f4d77e", "#95c36e", "#fbb875", "#91c4e9",
    "#f8c9c6", "#c5e3ea", "#e8d1e6", "#ffefbb", "#d3e4be", "#fddbb4", "#c1e3fa",
    "#912321", "#005d6e", "#6e2769", "#926d17", "#1e662a", "#ae450b", "#002b58",
]

def _auto_color(name, names, configured):
    """Configured colour if any; otherwise the next SELECT_ORDER colour not used by a configured
    name, handed out to the unconfigured names in sorted order."""
    if name in configured and configured[name]:
        return configured[name]
    used = {str(c).lower() for c in configured.values() if c}
    free = [c for c in _SELECT_ORDER if c not in used] or _SELECT_ORDER
    unconfigured = sorted(n for n in set(names) if n not in configured)
    return free[unconfigured.index(name) % len(free)]

def sample_type_color(sample_type):
    return _auto_color(sample_type, {SAMPLES[s]["sample_type"] for s in SAMPLES},
                       config.get("sample_type_colors", {}) or {})

def validate_ref_color(tissue):
    return _auto_color(tissue, config.get("validate_ref_tissues", []),
                       config.get("validate_ref_colors", {}) or {})

def sample_alias(sample):
    return SAMPLES[sample].get("alias") or sample

def _has_alias(sample):
    return bool(SAMPLES[sample].get("alias"))

def group_has_alias(group_id):
    return any(_has_alias(s) for s in GROUPS[group_id])

def bed_has_alias(cohort_id, bed_id):
    return any(_has_alias(s) for s in bed_samples(cohort_id, bed_id))

def alias_map_args(samples):
    return [str(s) + '=' + str(SAMPLES[s]['alias']) for s in samples if _has_alias(s)]

def sample_fraction_threshold(group_id, fraction):
    return ceil(len(GROUPS[group_id]) * fraction)

# Pipeline stages are on by default; set config[<key>] = False to skip one.
def flag(key):
    return config.get(key, True)

def _rule_threads(wc, rule_key):
    return int(SAMPLES[wc.sample].get((str(rule_key) + '_threads'), config["threads"]))

def _group_threads(group_id, rule_key, default):
    return int(config.get("groups", {}).get(group_id, {}).get((str(rule_key) + '_threads'), default))

def sample_hpo(sample):
    terms = SAMPLES[sample].get("hpo") or []
    if isinstance(terms, str):
        terms = terms.split(";")
    return [str(t).strip() for t in terms if str(t).strip()]

def amelie_tsv(sample):
    return str(SAMPLES[sample]["outdir"]) + "/phenotype/" + str(sample) + "_amelie.tsv"

def hpo_terms_tsv(sample):
    return str(SAMPLES[sample]["outdir"]) + "/phenotype/" + str(sample) + "_hpo_terms.tsv"

def hpo_obo_path():
    return config.get("hpo_obo") or os.path.join(workflow.basedir, "resources", "hpo_data", "hp.obo")

# Builds the full list of target files for `rule all`,
# i.e. every output the pipeline should produce given the enabled stages
def all_outputs():
    outs = []
    for s in SAMPLES:
        od = SAMPLES[s]["outdir"]

        if flag("longcallr"):
            outs.append((str(od) + '/variant_calling/longcallR/' + str(s) + '_longcallR_norm.vcf.gz'))

        if flag("nanots"):
            outs.append((str(od) + '/variant_calling/nanoTS/' + str(s) + '_nanoTS_norm.vcf.gz'))

        if flag("clair3_rna"):
            outs.append((str(od) + '/variant_calling/clair3_rna/' + str(s) + '_clair3_rna_norm.vcf.gz'))

        if flag("deepvariant"):
            outs.append((str(od) + '/variant_calling/deepvariant/' + str(s) + '_deepvariant_norm.vcf.gz'))

        if SAMPLES[s].get("wgs_vcf"):
            outs.append((str(od) + '/variant_calling/wgs/' + str(s) + '_wgs_phased.vcf.gz'))

        if sample_hpo(s):
            outs.append(amelie_tsv(s))

        if flag("compile_variants"):
            outs.append((str(od) + '/variant_calling/compiled_variants/' + str(s) + '_filtered_variants.tsv'))

        if flag("phase_reads"):
            outs.append((str(od) + '/phased_reads/' + str(s) + '_phasing_summary.tsv'))

        if flag("ase_analysis"):
            outs.append((str(od) + '/ase_analysis/' + str(s) + '_binomial_ase_results.tsv'))

        if flag("junction_analysis"):
            for t in sample_tissues(s):
                outs.append(
                    (str(od) + '/junction_analysis/gtex_' + str(t) + '/' + str(s) + '_gtex_' + str(t) + '_outlier_junctions.tsv')
                )

        if flag("qc"):
            outs.append((str(od) + '/qc/' + str(s) + '_on_target.tsv'))
            outs.append((str(od) + '/qc/' + str(s) + '_read_attributes.tsv'))
            outs.append((str(od) + '/qc/' + str(s) + '_full_length_ratio.tsv'))
            outs.append((str(od) + '/qc/' + str(s) + '_avg_coverage.tsv'))

        if flag("gene_quantification"):
            outs.append((str(od) + '/gene_quantification/' + str(s) + '_gene_count.tsv'))
            outs.append((str(od) + '/gene_quantification/' + str(s) + '_gene_coverage.tsv'))
            outs.append((str(od) + '/gene_quantification/' + str(s) + '_gene_assignment.tsv'))

    if flag("merge_results"):
        # Junction-derived outputs run for each filter set ("default" and "stringent")
        _jxn_filter_sets = ["default", "stringent"]
        for (cid, bid) in BED_GROUPS:
            bod = bed_outdir(cid, bid)
            for fs in _jxn_filter_sets:
                outs.append((str(bod) + '/' + str(fs) + '/hits_upset_density.pdf'))
            # Per-sample final reports (QC page + IGV-style hit pages), default hit set only
            if flag("hit_reports"):
                for s in bed_samples(cid, bid):
                    outs.append((str(bod) + '/default/final_reports/' + str(s) + '_report.pdf'))
                    outs.append((str(bod) + '/default/final_reports/' + str(s) + '_review.html'))
        for gid in GROUPS:
            god = group_outdir(gid)
            for fs in _jxn_filter_sets:
                for fname in ('genes_with_pathogenic_variant_barplot.pdf', 'genes_with_ASE_barplot.pdf',
                              'genes_with_outlier_junction_barplot.pdf', 'genes_with_RNA_dysregulation_barplot.pdf'):
                    outs.append((str(god) + '/merged_hits/' + str(fs) + '/' + fname))

        if flag("junction_analysis"):
            for gid in GROUPS:
                if cja_enabled(gid):
                    outs.append(str(group_outdir(gid)) + '/merged_hits/default/method_comparison/'
                                + str(GROUP_SAMPLE_TYPE[gid]) + '_method_comparison.pdf')

    if flag("cohort_junction_analysis"):
        for gid in GROUPS:
            god = group_outdir(gid)
            if not cja_enabled(gid):
                continue  # 8A always runs; 8B (fitting) and 8C (outlier calling) require at least min_n samples
            outs.append(
                god + "/cohort_junction_analysis/" + gid + "_" + _cja_thr_label(gid)
                + "/" + gid + "_outliers.tsv"
            )
            if group_has_alias(gid):
                outs.append(
                    god + "/cohort_junction_analysis/" + gid + "_" + _cja_thr_label(gid)
                    + "/" + gid + "_outliers_alias.tsv"
                )

    for gid in GROUPS:
        if ctx_enabled(gid):
            outs.extend(ctx_outputs(gid))

    if flag("softclip_analysis"):
        for gid in GROUPS:
            outs.extend(softclip_outputs(gid))

    if flag("gene_quantification"):
        # cohort-level per-gene CPTM / MOTR box plots (one PDF per gene, by sample type)
        for (cid, bid) in BED_GROUPS:
            for _m in ("count", "coverage", "assignment", "amalgam"):
                outs.append(str(bed_outdir(cid, bid)) + "/cohort_qc/gene_quantification/by_" + _m)
        for gid in GROUPS:
            god = group_outdir(gid)
            outs.append(god + "/gene_quantification/by_count/gene_count_matrix_cptm.tsv")
            outs.append(god + "/gene_quantification/by_count/gene_count_matrix_motr.tsv")
            outs.append(god + "/gene_quantification/by_coverage/gene_coverage_matrix_cptm.tsv")
            outs.append(god + "/gene_quantification/by_coverage/gene_coverage_matrix_motr.tsv")
            outs.append(god + "/gene_quantification/by_assignment/gene_assignment_matrix_cptm.tsv")
            outs.append(god + "/gene_quantification/by_assignment/gene_assignment_matrix_motr.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/gene_amalgam_matrix_raw_all_genes.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/gene_amalgam_matrix_raw.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/transcript_amalgam_matrix_raw_all_genes.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/transcript_amalgam_matrix_raw.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/gene_amalgam_matrix_cptm.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/gene_amalgam_matrix_motr.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/annotation/annotated.gtf.gz")
            if group_has_alias(gid):
                outs.append(god + "/gene_quantification/by_count/gene_count_matrix_cptm_alias.tsv")
                outs.append(god + "/gene_quantification/by_count/gene_count_matrix_motr_alias.tsv")
                outs.append(god + "/gene_quantification/by_coverage/gene_coverage_matrix_cptm_alias.tsv")
                outs.append(god + "/gene_quantification/by_coverage/gene_coverage_matrix_motr_alias.tsv")
                outs.append(god + "/gene_quantification/by_assignment/gene_assignment_matrix_cptm_alias.tsv")
                outs.append(god + "/gene_quantification/by_assignment/gene_assignment_matrix_motr_alias.tsv")
                outs.append(god + "/gene_quantification/by_amalgam/gene_amalgam_matrix_cptm_alias.tsv")
                outs.append(god + "/gene_quantification/by_amalgam/gene_amalgam_matrix_motr_alias.tsv")

    for (cid, bid) in BED_GROUPS:
        if any(SAMPLES[s].get("wgs_vcf") for s in bed_samples(cid, bid)):
            wcd = str(bed_outdir(cid, bid)) + '/cohort_qc/wgs_concordance/' + str(bid)
            outs.append(wcd + '_wgs_concordance_summary.tsv')
            outs.append(wcd + '_wgs_concordance_heatmap.pdf')

    if flag("qc"):
        for (cid, bid) in BED_GROUPS:
            bod = bed_outdir(cid, bid)
            cqd = bod + "/cohort_qc"
            outs.append((str(bod) + '/cohort_qc/validate_sample_types/' + str(bid) + '_distance_heatmap.pdf'))
            outs.append((str(cqd) + '/on_target_rates/' + str(bid) + '_on_target_rates.pdf'))
            outs.append((str(cqd) + '/read_attributes/' + str(bid) + '_read_lengths.pdf'))
            for _m in ('avg_coverage', 'median_coverage', 'breadth_coverage'):
                outs.append((str(cqd) + '/avg_coverage/' + str(bid) + '_' + _m + '_matrix.tsv'))
                outs.append((str(cqd) + '/avg_coverage/' + str(bid) + '_' + _m + '_heatmap.pdf'))
            outs.append((str(cqd) + '/full_length_ratio/' + str(bid) + '_full_length_ratio_matrix.tsv'))
            outs.append((str(cqd) + '/full_length_ratio/' + str(bid) + '_full_length_ratio_heatmap.pdf'))
            if bed_has_alias(cid, bid):
                outs.append((str(bod) + '/cohort_qc/validate_sample_types/' + str(bid) + '_distance_heatmap_alias.pdf'))
                outs.append((str(cqd) + '/on_target_rates/' + str(bid) + '_on_target_rates_alias.pdf'))
                outs.append((str(cqd) + '/read_attributes/' + str(bid) + '_read_lengths_alias.pdf'))
                for _m in ('avg_coverage', 'median_coverage', 'breadth_coverage'):
                    outs.append((str(cqd) + '/avg_coverage/' + str(bid) + '_' + _m + '_matrix_alias.tsv'))
                    outs.append((str(cqd) + '/avg_coverage/' + str(bid) + '_' + _m + '_heatmap_alias.pdf'))
                outs.append((str(cqd) + '/full_length_ratio/' + str(bid) + '_full_length_ratio_matrix_alias.tsv'))

    return outs

rule all:
    input:
        all_outputs()


# Rule definitions, split by pipeline stage and applied in roughly this order
include: "rules/1_call_variants.smk"
include: "rules/2_compile_variants.smk"
include: "rules/3_phase_reads.smk"
include: "rules/4_ase_analysis.smk"
include: "rules/5_junction_analysis.smk"
include: "rules/6_sample_qc.smk"
include: "rules/7_sample_gene_quantification.smk"
include: "rules/8_cohort_junction_analysis.smk"
include: "rules/9_cohort_transcript_analysis.smk"
include: "rules/10_merge_results.smk"
include: "rules/11_softclip_analysis.smk"
