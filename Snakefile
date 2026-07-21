"""
RNA-Dx Snakemake Pipeline
Multi-sample wrapper for variant calling, phasing, ASE, and splice junction analysis.

Note: `configfile:` is not declared here, since that path is resolved
relative to Snakemake's working directory, not this file's location.
The config is instead supplied via `--configfile` on the command line
(see submit_snakemake.sh).
"""

import yaml
import os
import re
from pathlib import Path
from collections import defaultdict
from math import ceil
import shlex

# Ensure pipe failures aren't masked (e.g. `cmd 2>&1 | tee {log}` would
# otherwise report tee's exit status instead of cmd's).
shell.prefix("set -euo pipefail;")


# ---------------------------------------------------------------------------
# Load per-run sample config (passed via --config run=<path>)
# ---------------------------------------------------------------------------
with open(config["run"]) as fh:
    run_cfg = yaml.safe_load(fh)

SAMPLES = run_cfg["samples"]

# Merge every other top-level run-config key over config.yaml's defaults
# (e.g. the "Pipeline stage booleans" block, or any run-specific override).
# samples/merged_outdir keep their own special-cased handling below.
for _key, _val in run_cfg.items():
    if _key in ("samples", "merged_outdir"):
        continue
    config[_key] = _val

# merged_outdir is cohort-specific, so it can be set either on the command
# line (--config merged_outdir=<path>, takes priority) or as a top-level
# key in the run config YAML alongside "samples:".
if not config.get("merged_outdir"):
    config["merged_outdir"] = run_cfg.get("merged_outdir", "")

if not config.get("merged_outdir"):
    raise ValueError(
        "merged_outdir must be set, either via --config merged_outdir=<path> or as a "
        "top-level 'merged_outdir:' key in the run config YAML (alongside 'samples:'), "
        "for the merge_hits stage to know where to write cross-sample merged results."
    )

# ---------------------------------------------------------------------------
# Resolve pipeline-relative reference-data / environment paths against the
# pipeline directory itself (workflow.basedir), so config.yaml can use
# relative paths (e.g. "resources/gnomad_data") instead of machine-specific
# absolute ones -- keeping the whole folder self-contained and relocatable.
# Absolute paths and remote URLs in config.yaml are left untouched.
# ---------------------------------------------------------------------------
_BUNDLED_PATH_KEYS = [
    "genome", "annotation",
    "conda_env", "conda_env_compile_variants", "gnomad_base", "clinvar_vcf", "annovar_dir",
    "cadd_data_dir", "cadd_script", "cadd_local_prescored_snv", "cadd_local_prescored_indel",
    "gnomad_mito_vcf", "longcallr_bin",
    "nanots_model_unphased", "nanots_model_phased", "gtex_data_dir", "omim_file",
    "spliceai_prescored_snv_vcf", "spliceai_prescored_indel_vcf",
]
_URL_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*://")

def _is_remote_sentinel(value):
    """True if a config value is the literal keyword 'remote' (case-insensitive,
    whitespace-insensitive) -- an explicit request to use the public HTTPS
    resource instead of a local copy. See _resolved_*() helpers below."""
    return isinstance(value, str) and value.strip().lower() == "remote"

for _key in _BUNDLED_PATH_KEYS:
    _val = config.get(_key)
    if _val and not os.path.isabs(_val) and not _URL_RE.match(_val) and not _is_remote_sentinel(_val):
        config[_key] = os.path.join(workflow.basedir, _val)

# ---------------------------------------------------------------------------
# Canonical public HTTPS locations for gnomAD/ClinVar/CADD, used when a
# config value is the "remote" sentinel above, and as the fallback
# compile_variants.py retries against if a configured local copy fails.
# ---------------------------------------------------------------------------
_REMOTE_GNOMAD_BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes"
_REMOTE_GNOMAD_MITO_VCF = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/"
                           "vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz")
_REMOTE_CLINVAR_VCF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
_REMOTE_CADD_PRESCORED_URL = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"

def _resolved_gnomad_base():
    """config['gnomad_base'], or the canonical public gnomAD base URL if the
    config value is the 'remote' sentinel."""
    return _REMOTE_GNOMAD_BASE if _is_remote_sentinel(config["gnomad_base"]) else config["gnomad_base"]

def _resolved_gnomad_mito_vcf():
    """config['gnomad_mito_vcf'], or the canonical public gnomAD v3.1 mito URL
    if gnomad_base is the 'remote' sentinel -- mito follows the same
    local/remote switch as the autosomal/X/Y chromosomes."""
    return _REMOTE_GNOMAD_MITO_VCF if _is_remote_sentinel(config["gnomad_base"]) else config["gnomad_mito_vcf"]

def _resolved_clinvar_vcf():
    return _REMOTE_CLINVAR_VCF if _is_remote_sentinel(config["clinvar_vcf"]) else config["clinvar_vcf"]

def _cadd_use_local():
    """False if config['cadd_script'] is the 'remote' sentinel -- i.e. the
    user wants to skip the local CADD-scripts install and rely on the
    remote pre-scored SNV lookup instead."""
    return not _is_remote_sentinel(config["cadd_script"])

# ---------------------------------------------------------------------------
# Derive the set of tissues across all samples (for junction outlier rules)
# ---------------------------------------------------------------------------
def _parse_tissues(raw):
    """Normalize a sample's 'tissues' value to a list of tissue names.
    Accepts either a YAML list (parsed directly by PyYAML) or a plain
    comma-separated string, for backward compatibility."""
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


# ---------------------------------------------------------------------------
# Group samples by (bed, sample_type) for cross-sample merging (rules/merge_hits.smk).
# Each sample's run-config entry must include a "sample_type" field (e.g. "fibroblasts",
# "PBMC") alongside the existing "bed" field. Samples sharing both the same BED panel
# and the same sample_type are merged together, mirroring how merge_candidate_hits.sh
# was previously run once per (panel, sample_type) combination.
# ---------------------------------------------------------------------------
def _bed_id(bed):
    """Filesystem-safe identifier for a BED panel, e.g. 'IEI422_gene_symbols'."""
    return Path(bed).stem


def _group_id(bed, sample_type):
    """Filesystem-safe identifier for a (bed, sample_type) group, e.g. 'IEI422_fibroblasts'."""
    return f"{_bed_id(bed)}_{sample_type}"


def _group_id_from_ids(bed_id, sample_type):
    """Reconstruct a group_id from its already-split bed_id/sample_type wildcards
    (used by rules whose output path has {bed_id}/{sample_type} as separate
    wildcards, to look back up into GROUPS/GROUP_* dicts keyed by group_id)."""
    return f"{bed_id}_{sample_type}"


def all_groups():
    """Return {group_id: [sample, sample, ...]} for every unique (bed, sample_type)
    combination present in the run config. Also populates GROUP_BED_ID and
    GROUP_SAMPLE_TYPE, since group_id is just a concatenation and shouldn't be
    re-split (bed stems or sample_types could themselves contain underscores)."""
    groups = defaultdict(list)
    for s in SAMPLES:
        if "sample_type" not in SAMPLES[s]:
            raise ValueError(
                f"Sample '{s}' is missing a 'sample_type' field in the run config, "
                f"required for grouping samples during the merge_hits stage."
            )
        bed = SAMPLES[s]["bed"]
        sample_type = SAMPLES[s]["sample_type"]
        gid = _group_id(bed, sample_type)
        groups[gid].append(s)
        GROUP_BED_ID[gid] = _bed_id(bed)
        GROUP_SAMPLE_TYPE[gid] = sample_type
    return dict(groups)


GROUP_BED_ID = {}       # {group_id: bed_id}
GROUP_SAMPLE_TYPE = {}  # {group_id: sample_type}
GROUPS = all_groups()  # {group_id: [sample, ...]}


def group_tissues(group_id):
    """Union of tissues across all samples in a group."""
    tissues = set()
    for s in GROUPS[group_id]:
        tissues.update(sample_tissues(s))
    return sorted(tissues)


def group_outdir(group_id):
    """Shared output directory for a group's merged results, nested under its BED
    panel: config['merged_outdir']/{bed_id}/{sample_type}."""
    return f"{config['merged_outdir']}/{GROUP_BED_ID[group_id]}/{GROUP_SAMPLE_TYPE[group_id]}"


# ---------------------------------------------------------------------------
# Group (bed, sample_type) groups by BED panel alone, for stages that operate
# across all sample types on the same panel: validating sample types against
# GTEx reference tissues, and the final cross-sample-type merge of all_hits.
# ---------------------------------------------------------------------------
def all_bed_groups():
    """Return {bed_id: [group_id, ...]} for every BED panel present in the run,
    where each group_id is one of the (bed, sample_type) keys in GROUPS."""
    groups = defaultdict(list)
    for gid, members in GROUPS.items():
        bed = SAMPLES[members[0]]["bed"]
        groups[_bed_id(bed)].append(gid)
    return dict(groups)


BED_GROUPS = all_bed_groups()  # {bed_id: [group_id, ...]}


def bed_path(bed_id):
    """Actual BED file path for a given bed_id."""
    gid = BED_GROUPS[bed_id][0]
    return SAMPLES[GROUPS[gid][0]]["bed"]


def bed_outdir(bed_id):
    """Shared output directory for a BED panel's cross-sample-type results."""
    return f"{config['merged_outdir']}/{bed_id}"


def _quoted(items):
    """Shell-quote each item in a list (e.g. hex colors like '#8BBF9F' would
    otherwise be treated as a bash comment when passed unquoted in a shell: block)."""
    return [shlex.quote(str(x)) for x in items]


_DEFAULT_SAMPLE_TYPE_PALETTE = ["#8BBF9F", "#D27D7D", "#A78BC5", "#E8B04B", "#4A7C9B", "#C46B6B"]


def sample_type_color(sample_type):
    """Color for a sample_type in validate_sample_types plots. Looks up
    config['sample_type_colors'][sample_type] first; otherwise assigns a
    deterministic color from a default palette based on sorted sample_type order."""
    configured = config.get("sample_type_colors", {})
    if sample_type in configured:
        return configured[sample_type]
    all_types = sorted({SAMPLES[s]["sample_type"] for s in SAMPLES})
    idx = all_types.index(sample_type) % len(_DEFAULT_SAMPLE_TYPE_PALETTE)
    return _DEFAULT_SAMPLE_TYPE_PALETTE[idx]


def sample_fraction_threshold(group_id, fraction):
    """Round-up sample-count threshold for a given fraction of a group's sample size
    (mirrors the bash script's `echo "(n * frac + 0.9999)/1" | bc` rounding)."""
    return ceil(len(GROUPS[group_id]) * fraction)


# ---------------------------------------------------------------------------
# Helper: booleans from config (default True)
# ---------------------------------------------------------------------------
def flag(key):
    return config.get(key, True)


# ---------------------------------------------------------------------------
# Helper: per-sample, per-rule thread count.
#
# Rules look up their thread count with:
#   threads: lambda wc: _rule_threads(wc, "rule_key")
#
# Resolution order:
#   1. run_config.yaml sample entry  e.g. longcallr_threads: 4
#   2. config/config.yaml            threads: 8   (global default)
#
# The following rules are always single-threaded and do NOT use this helper:
#   detect_ase_outliers, get_junction_counts
# ---------------------------------------------------------------------------
def _rule_threads(wc, rule_key):
    """Return the thread count for a given rule and sample.
    Looks for <rule_key>_threads in the sample's run-config entry first,
    then falls back to the global config["threads"] default."""
    return int(SAMPLES[wc.sample].get(f"{rule_key}_threads", config["threads"]))


# ---------------------------------------------------------------------------
# Per-group / per-bed-panel thread & memory overrides.
#
# The group- and bed-level rules (build_group_junction_matrix,
# run_cohort_junction_outlier_analysis, the merge_hits.smk stages --
# merge_group_variants, merge_group_ase, merge_group_junctions,
# split_group_hits_by_sample, merge_sample_hits, concat_group_hits,
# plot_group_hits -- final_merge, validate_sample_types) don't have a
# single sample to key off
# of the way _rule_threads() does above -- they run once per (bed,
# sample_type) group, or once per bed panel. Overrides for these are looked
# up by group_id (e.g. "IEI422_gene_symbols_fibroblasts") or bare bed_id
# (for the two rules that operate per-BED-panel rather than per-group), via
# a top-level `groups:` block in the run config YAML:
#
#   groups:
#     IEI422_gene_symbols_fibroblasts:            # group-level: bed_id_sample_type
#       build_group_junction_matrix_threads: 1
#       build_group_junction_matrix_mem_gb: 40
#       run_cohort_junction_outlier_analysis_threads: 16
#       run_cohort_junction_outlier_analysis_mem_gb: 96
#     IEI422_gene_symbols:                        # bed-level: bare bed_id
#       validate_sample_types_mem_gb: 300
#
# Resolution order: groups.<id>.<rule_key>_<field> in the run config, else
# whatever default the rule itself passes in.
# ---------------------------------------------------------------------------
def _group_threads(group_id, rule_key, default):
    return int(config.get("groups", {}).get(group_id, {}).get(f"{rule_key}_threads", default))

def _group_mem_gb(group_id, rule_key, default_gb):
    return float(config.get("groups", {}).get(group_id, {}).get(f"{rule_key}_mem_gb", default_gb))


# ---------------------------------------------------------------------------
# Collect all final outputs across samples
# ---------------------------------------------------------------------------
def all_outputs():
    outs = []
    for s in SAMPLES:
        od = SAMPLES[s]["outdir"]

        if flag("longcallr"):
            outs.append(f"{od}/variant_calling/longcallR/{s}_longcallR_norm.vcf.gz")

        if flag("nanots"):
            outs.append(f"{od}/variant_calling/nanoTS/{s}_nanoTS_norm.vcf.gz")

        if flag("clair3_rna"):
            outs.append(f"{od}/variant_calling/clair3_rna/{s}_clair3_rna_norm.vcf.gz")

        if flag("deepvariant"):
            outs.append(f"{od}/variant_calling/deepvariant/{s}_deepvariant_norm.vcf.gz")

        if flag("compile_variants"):
            outs.append(f"{od}/variant_calling/compiled_variants/{s}_compiled_variants.tsv")

        if flag("phase_reads"):
            outs.append(f"{od}/phased_reads/{s}_phasing_summary.tsv")

        if flag("ase_analysis"):
            outs.append(f"{od}/ase_analysis/{s}_binomial_ase_results.tsv")

        if flag("junction_analysis"):
            # Requesting the chain's final output (identify_junction_outliers)
            # pulls get_junction_counts -> perform_binomial_tests along with it,
            # since this flag now gates the whole junction_analysis.smk chain
            # as a single unit.
            for t in sample_tissues(s):
                outs.append(
                    f"{od}/junction_analysis/gtex_{t}/{s}_gtex_{t}_outlier_junctions.tsv"
                )

    if flag("merge_hits"):
        # Requesting final_merge's output pulls the whole merge_hits.smk chain
        # along with it (merge_group_variants -> merge_group_ase ->
        # merge_group_junctions -> split_group_hits_by_sample ->
        # merge_sample_hits -> concat_group_hits -> plot_group_hits), since
        # this flag gates that whole chain as a unit.
        for bid in BED_GROUPS:
            bod = bed_outdir(bid)
            outs.append(f"{bod}/merged_all_hits.tsv")

    if flag("cohort_junction_analysis"):
        for gid in GROUPS:
            god = group_outdir(gid)
            outs.append(
                f"{god}/cohort_junction_analysis/{gid}_padj{config['padj_threshold']}"
                f"_delta{config['delta_psi_threshold']}/{gid}_outliers.tsv"
            )

    if flag("validate_sample_types"):
        # Requesting validate_sample_types' output pulls build_group_junction_matrix
        # along with it, since this flag gates the whole validate_sample_types.smk chain.
        for bid in BED_GROUPS:
            bod = bed_outdir(bid)
            outs.append(f"{bod}/validate_sample_types/{bid}_distance_heatmap.pdf")

    return outs


rule all:
    input:
        all_outputs()


# ---------------------------------------------------------------------------
# Include modular rule files
# ---------------------------------------------------------------------------
include: "rules/longcallr.smk"
include: "rules/nanots.smk"
include: "rules/clair3_rna.smk"
include: "rules/deepvariant.smk"
include: "rules/compile_variants.smk"
include: "rules/phase_reads.smk"
include: "rules/ase_analysis.smk"
include: "rules/junction_analysis.smk"
include: "rules/cohort_junction_analysis.smk"
include: "rules/merge_hits.smk"
include: "rules/validate_sample_types.smk"
