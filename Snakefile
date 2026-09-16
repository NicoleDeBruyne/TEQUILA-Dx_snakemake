"""
RNA-Dx Snakemake Pipeline
Runs variant calling, phasing, allele-specific expression (ASE), and splice junction analysis across multiple samples.
See docs/general.md for an overview of the pipeline structure.
Config is supplied via --config/--configfile on the command line (see README.md).
"""

import yaml
import os
import re
from pathlib import Path
from collections import defaultdict
from math import ceil
import shlex

shell.prefix("set -euo pipefail;")


# Load per-run sample config (passed via --config run=<path>)
with open(config["run"]) as fh:
    run_cfg = yaml.safe_load(fh)

SAMPLES = run_cfg["samples"]

# Override the defaults in config.yaml
for _key, _val in run_cfg.items():
    if _key in ("samples", "output_dir"):
        continue
    config[_key] = _val

# Get output_dir from the command line (--config output_dir=<path>, highest priority) or from the run config YAML
if not config.get("output_dir"):
    config["output_dir"] = run_cfg.get("output_dir", "")

if not config.get("output_dir"):
    raise ValueError(
        "output_dir must be set, either via --config output_dir=<path> or as a "
        "top-level 'output_dir:' key in the run config YAML (alongside 'samples:'), "
        "for the pipeline to know where to write its output layout."
    )

# By default each sample writes to {output_dir}/samples/{sample}/output.
# A sample can override this with its own 'outdir' key in the run config
for _sample, _entry in SAMPLES.items():
    if not _entry.get("outdir"):
        _entry["outdir"] = config["output_dir"] + "/samples/" + _sample + "/output"

# Turn relative reference-data/environment paths into absolute paths, anchored to the pipeline's own folder
_BUNDLED_PATH_KEYS = [
    "genome", "annotation",
    "conda_env", "conda_env_compile_variants", "gnomad_base", "clinvar_vcf", "annovar_dir",
    "amalgam_dir",
    "cadd_data_dir", "cadd_script", "cadd_local_prescored_snv", "cadd_local_prescored_indel",
    "gnomad_mito_vcf", "longcallr_bin",
    "nanots_model_unphased", "nanots_model_phased", "gtex_data_dir", "omim_file",
    "spliceai_prescored_snv_vcf", "spliceai_prescored_indel_vcf",
]
_URL_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*://")

def _is_remote_sentinel(value):
    """True if a config value is exactly the word 'remote', meaning
    "use the public online copy of this file instead of a local one"."""
    return isinstance(value, str) and value.strip().lower() == "remote"

for _key in _BUNDLED_PATH_KEYS:
    _val = config.get(_key)
    if _val and not os.path.isabs(_val) and not _URL_RE.match(_val) and not _is_remote_sentinel(_val):
        config[_key] = os.path.join(workflow.basedir, _val)

# Public download links for gnomAD/ClinVar/CADD, used when a config value is set to "remote" instead of a local file path
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
    """True unless config['cadd_script'] is set to "remote", in which case
    CADD scoring is done via the public pre-scored lookup instead of a local install."""
    return not _is_remote_sentinel(config["cadd_script"])

def _parse_tissues(raw):
    """Turn a sample's 'tissues' value into a plain list of tissue names."""
    if isinstance(raw, list):
        return [str(t).strip() for t in raw]
    return [t.strip() for t in str(raw).strip("[]").split(",")]

def all_tissues():
    """All tissue names used by any sample, sorted and without duplicates."""
    tissues = set()
    for s in SAMPLES.values():
        tissues.update(_parse_tissues(s["tissues"]))
    return sorted(tissues)

TISSUES = all_tissues()

def sample_tissues(sample):
    """Tissue names for one sample."""
    return _parse_tissues(SAMPLES[sample]["tissues"])

# Named cohorts: "cohort_all" always exists and includes every sample,
# Extra cohorts can be declared in the run config to analyze a chosen subset of samples together, e.g.:
#
#   cohorts:
#     cohort_subset:
#       - SAMPLE_A
#       - SAMPLE_C
#
# A sample can belong to as many cohorts as you like; it is always part of cohort_all regardless. "cohort_all" cannot be redefined.
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
                    + 'required for grouping samples during the merge_hits stage.'
                )
            bed = SAMPLES[s]["bed"]
            sample_type = SAMPLES[s]["sample_type"]
            gid = _group_id(cid, bed, sample_type)
            groups[gid].append(s)
            GROUP_COHORT_ID[gid] = cid
            GROUP_BED_ID[gid] = _bed_id(bed)
            GROUP_SAMPLE_TYPE[gid] = sample_type
    return dict(groups)

GROUP_COHORT_ID = {}    # {group_id: cohort_id}
GROUP_BED_ID = {}       # {group_id: bed_id}
GROUP_SAMPLE_TYPE = {}  # {group_id: sample_type}
GROUPS = all_groups()  # {group_id: [sample, ...]}

def group_tissues(group_id):
    tissues = set()
    for s in GROUPS[group_id]:
        tissues.update(sample_tissues(s))
    return sorted(tissues)

def group_outdir(group_id):
    return (str(config['output_dir']) + '/' + str(GROUP_COHORT_ID[group_id]) + '/' + str(GROUP_BED_ID[group_id]) + '/output/sample_types/' + str(GROUP_SAMPLE_TYPE[group_id]) + '/output')

def _cja_method_for_group(group_id):
    """Return "beta_binomial" or "modified_zscore" for a group, based on
    config["cohort_jxn_method"] or, if "auto", the group's sample count."""
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
    """Folder-name label for this group's outlier threshold, e.g.
    'padj0.05_delta0.1' or 'z3.5_delta0.1', depending on the method used."""
    if _cja_method_for_group(group_id) == "beta_binomial":
        return "padj" + str(config["cohort_jxn_beta_padj_threshold"]) + "_delta" + str(config["delta_psi_threshold"])
    z = config.get("cohort_jxn_z_threshold", 3.5)
    d = config.get("cohort_jxn_z_delta_threshold", 0.1)
    return "z" + str(z) + "_delta" + str(d)

def _cja_outliers_filtered_path(group_id):
    thr_label = _cja_thr_label(group_id)
    return (group_outdir(group_id) + "/cohort_junction_analysis/" + group_id + "_" + thr_label
            + "/" + group_id + "_outliers_filtered.tsv")

def _cja_thr_flag(group_id):
    if _cja_method_for_group(group_id) == "beta_binomial":
        return "--bb-thresholds " + str(config["cohort_jxn_beta_padj_threshold"]) + ":" + str(config["delta_psi_threshold"])
    z = config.get("cohort_jxn_z_threshold", 3.5)
    d = config.get("cohort_jxn_z_delta_threshold", 0.1)
    return "--z-thresholds " + str(z) + ":" + str(d)

def _cja_n_threshold(group_id):
    if _cja_method_for_group(group_id) == "beta_binomial":
        return config.get("cohort_jxn_beta_n_threshold", 30)
    return config.get("cohort_jxn_zscore_n_threshold", 10)

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
    """Shell-quote each item in a list, so values like hex colors ('#8BBF9F') aren't accidentally treated as a bash comment."""
    return [shlex.quote(str(x)) for x in items]

_DEFAULT_SAMPLE_TYPE_PALETTE = ["#8BBF9F", "#D27D7D", "#A78BC5", "#E8B04B", "#4A7C9B", "#C46B6B"]

def sample_type_color(sample_type):
    """Plot color for a sample_type: uses config['sample_type_colors'] if set, otherwise picks one from a default color palette."""
    configured = config.get("sample_type_colors", {})
    if sample_type in configured:
        return configured[sample_type]
    all_types = sorted({SAMPLES[s]["sample_type"] for s in SAMPLES})
    idx = all_types.index(sample_type) % len(_DEFAULT_SAMPLE_TYPE_PALETTE)
    return _DEFAULT_SAMPLE_TYPE_PALETTE[idx]

def sample_alias(sample):
    """A sample's configured alias, or its own ID if it has none."""
    return SAMPLES[sample].get("alias") or sample

def _has_alias(sample):
    return bool(SAMPLES[sample].get("alias"))

def group_has_alias(group_id):
    """True if any sample in this sample_type group has an alias configured."""
    return any(_has_alias(s) for s in GROUPS[group_id])

def bed_has_alias(cohort_id, bed_id):
    """True if any sample sharing this cohort+bed panel has an alias configured."""
    return any(_has_alias(s) for s in bed_samples(cohort_id, bed_id))

def alias_map_args(samples):
    """--alias-map CLI tokens ('sample=alias') for every sample in `samples`
    that actually has an alias configured (samples without one are simply
    omitted -- the receiving script's resolve() falls back to the sample ID)."""
    return [str(s) + '=' + str(SAMPLES[s]['alias']) for s in samples if _has_alias(s)]

def sample_fraction_threshold(group_id, fraction):
    return ceil(len(GROUPS[group_id]) * fraction)

def flag(key):
    return config.get(key, True)

def _rule_threads(wc, rule_key):
    return int(SAMPLES[wc.sample].get((str(rule_key) + '_threads'), config["threads"]))

def _group_threads(group_id, rule_key, default):
    return int(config.get("groups", {}).get(group_id, {}).get((str(rule_key) + '_threads'), default))

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

    if flag("merge_hits"):
        for (cid, bid) in BED_GROUPS:
            bod = bed_outdir(cid, bid)
            outs.append((str(bod) + '/hits_upset_density.pdf'))
            outs.append((str(bod) + '/merged_all_hits_simplified.tsv'))
            if bed_has_alias(cid, bid):
                outs.append((str(bod) + '/merged_all_hits_simplified_alias.tsv'))
        for gid in GROUPS:
            god = group_outdir(gid)
            for fname in ('genes_with_pathogenic_variant_boxplot.pdf', 'genes_with_ASE_boxplot.pdf',
                          'genes_with_outlier_junction_boxplot.pdf', 'genes_with_RNA_dysregulation_boxplot.pdf'):
                outs.append((str(god) + '/merged_hits/' + fname))

    if flag("cohort_junction_analysis"):
        for gid in GROUPS:
            god = group_outdir(gid)
            outs.append(
                god + "/cohort_junction_analysis/" + gid + "_" + _cja_thr_label(gid)
                + "/" + gid + "_outliers.tsv"
            )
            if group_has_alias(gid):
                outs.append(
                    god + "/cohort_junction_analysis/" + gid + "_" + _cja_thr_label(gid)
                    + "/" + gid + "_outliers_alias.tsv"
                )

    if flag("quantify_genes"):
        for gid in GROUPS:
            god = group_outdir(gid)
            outs.append(god + "/gene_quantification/by_count/gene_count_matrix.tsv")
            outs.append(god + "/gene_quantification/by_coverage/gene_coverage_matrix.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/quantification/gene_amalgam_gene_matrix.tsv")
            outs.append(god + "/gene_quantification/by_amalgam/annotation/annotated.gtf.gz")
            outs.append(god + "/gene_quantification/by_assignment/gene_assignment_matrix.tsv")
            if group_has_alias(gid):
                outs.append(god + "/gene_quantification/by_count/gene_count_matrix_alias.tsv")
                outs.append(god + "/gene_quantification/by_coverage/gene_coverage_matrix_alias.tsv")
                outs.append(god + "/gene_quantification/by_amalgam/quantification/gene_amalgam_gene_matrix_alias.tsv")
                outs.append(god + "/gene_quantification/by_assignment/gene_assignment_matrix_alias.tsv")

    if flag("cohort_qc"):
        for (cid, bid) in BED_GROUPS:
            bod = bed_outdir(cid, bid)
            cqd = bod + "/cohort_qc"
            outs.append((str(bod) + '/cohort_qc/validate_sample_types/' + str(bid) + '_distance_heatmap.pdf'))
            outs.append((str(cqd) + '/on_target_rates/' + str(bid) + '_on_target_rates.pdf'))
            outs.append((str(cqd) + '/read_attributes/' + str(bid) + '_read_lengths.pdf'))
            outs.append((str(cqd) + '/full_length_ratio/' + str(bid) + '_full_length_ratio_matrix.tsv'))
            outs.append((str(cqd) + '/full_length_ratio/' + str(bid) + '_full_length_ratio_heatmap.pdf'))
            if bed_has_alias(cid, bid):
                outs.append((str(bod) + '/cohort_qc/validate_sample_types/' + str(bid) + '_distance_heatmap_alias.pdf'))
                outs.append((str(cqd) + '/on_target_rates/' + str(bid) + '_on_target_rates_alias.pdf'))
                outs.append((str(cqd) + '/read_attributes/' + str(bid) + '_read_lengths_alias.pdf'))
                outs.append((str(cqd) + '/full_length_ratio/' + str(bid) + '_full_length_ratio_matrix_alias.tsv'))

    return outs

rule all:
    input:
        all_outputs()


# Load the rule definitions
include: "rules/1_call_variants.smk"
include: "rules/2_compile_variants.smk"
include: "rules/3_phase_reads.smk"
include: "rules/4_ase_analysis.smk"
include: "rules/5_junction_analysis.smk"
include: "rules/6_sample_qc.smk"
include: "rules/7_sample_gene_quantification.smk"
include: "rules/8_cohort_junction_analysis.smk"
include: "rules/9_merge_results.smk"
