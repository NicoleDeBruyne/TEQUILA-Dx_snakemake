"""
RNA-Dx Snakemake Pipeline
Multi-sample wrapper for variant calling, phasing, ASE, and splice junction analysis.
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

# Ensure pipe failures aren't masked (e.g. tee's exit status doesn't
# shadow the real command's).
shell.prefix("set -euo pipefail;")


# ---------------------------------------------------------------------------
# Load per-run sample config (passed via --config run=<path>)
# ---------------------------------------------------------------------------
with open(config["run"]) as fh:
    run_cfg = yaml.safe_load(fh)

SAMPLES = run_cfg["samples"]

# Merge every other top-level run-config key over config.yaml's defaults.
# samples/merged_outdir keep their own special-cased handling below.
for _key, _val in run_cfg.items():
    if _key in ("samples", "merged_outdir"):
        continue
    config[_key] = _val

# merged_outdir can be set via --config merged_outdir=<path> (takes
# priority) or as a top-level key in the run config YAML.
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
# pipeline directory itself (workflow.basedir). Absolute paths, remote
# URLs, and the "remote" sentinel are left untouched. See docs/general.md.
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
    """True if a config value is the literal keyword 'remote' -- an explicit
    request to use the public HTTPS resource instead of a local copy."""
    return isinstance(value, str) and value.strip().lower() == "remote"

for _key in _BUNDLED_PATH_KEYS:
    _val = config.get(_key)
    if _val and not os.path.isabs(_val) and not _URL_RE.match(_val) and not _is_remote_sentinel(_val):
        config[_key] = os.path.join(workflow.basedir, _val)

# ---------------------------------------------------------------------------
# Canonical public HTTPS locations for gnomAD/ClinVar/CADD -- used for the
# "remote" sentinel, and as compile_variants.py's fallback if a local copy fails.
# ---------------------------------------------------------------------------
_REMOTE_GNOMAD_BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes"
_REMOTE_GNOMAD_MITO_VCF = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/"
                           "vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz")
_REMOTE_CLINVAR_VCF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
_REMOTE_CADD_PRESCORED_URL = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"

def _resolved_gnomad_base():
    """config['gnomad_base'], or the canonical public gnomAD base URL if
    the config value is the 'remote' sentinel."""
    return _REMOTE_GNOMAD_BASE if _is_remote_sentinel(config["gnomad_base"]) else config["gnomad_base"]

def _resolved_gnomad_mito_vcf():
    """config['gnomad_mito_vcf'], or the canonical public gnomAD v3.1 mito
    URL if gnomad_base is the 'remote' sentinel."""
    return _REMOTE_GNOMAD_MITO_VCF if _is_remote_sentinel(config["gnomad_base"]) else config["gnomad_mito_vcf"]

def _resolved_clinvar_vcf():
    return _REMOTE_CLINVAR_VCF if _is_remote_sentinel(config["clinvar_vcf"]) else config["clinvar_vcf"]

def _cadd_use_local():
    """False if config['cadd_script'] is the 'remote' sentinel (skip the
    local CADD-scripts install, use the remote pre-scored lookup instead)."""
    return not _is_remote_sentinel(config["cadd_script"])

# ---------------------------------------------------------------------------
# Derive the set of tissues across all samples (for junction outlier rules)
# ---------------------------------------------------------------------------
def _parse_tissues(raw):
    """Normalize a sample's 'tissues' value to a list of tissue names.
    Accepts either a YAML list or a plain comma-separated string."""
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
# Group samples by (bed, sample_type) for cross-sample merging (rules/6_merge_hits.smk).
# Each sample's run-config entry must include a "sample_type" field. Samples
# sharing both the same BED panel and sample_type are merged together.
# ---------------------------------------------------------------------------
def _bed_id(bed):
    """Filesystem-safe identifier for a BED panel, e.g. 'IEI422_gene_symbols'."""
    return Path(bed).stem


def _group_id(bed, sample_type):
    """Filesystem-safe identifier for a (bed, sample_type) group, e.g. 'IEI422_fibroblasts'."""
    return f"{_bed_id(bed)}_{sample_type}"


def _group_id_from_ids(bed_id, sample_type):
    """Reconstruct a group_id from its already-split bed_id/sample_type wildcards."""
    return f"{bed_id}_{sample_type}"


def all_groups():
    """Return {group_id: [sample, sample, ...]} for every unique (bed, sample_type)
    combination present in the run config. Also populates GROUP_BED_ID and
    GROUP_SAMPLE_TYPE (group_id is a concatenation and shouldn't be re-split,
    since bed stems or sample_types could themselves contain underscores)."""
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
# across all sample types on the same panel: validating sample types, and
# the final cross-sample-type merge of all_hits.
# ---------------------------------------------------------------------------
def all_bed_groups():
    """Return {bed_id: [group_id, ...]} for every BED panel present in the run."""
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
    """Shell-quote each item in a list (hex colors like '#8BBF9F' would
    otherwise be treated as a bash comment when unquoted)."""
    return [shlex.quote(str(x)) for x in items]


_DEFAULT_SAMPLE_TYPE_PALETTE = ["#8BBF9F", "#D27D7D", "#A78BC5", "#E8B04B", "#4A7C9B", "#C46B6B"]


def sample_type_color(sample_type):
    """Color for a sample_type in validate_sample_types plots. Uses
    config['sample_type_colors'][sample_type] if set, else a deterministic
    color from a default palette."""
    configured = config.get("sample_type_colors", {})
    if sample_type in configured:
        return configured[sample_type]
    all_types = sorted({SAMPLES[s]["sample_type"] for s in SAMPLES})
    idx = all_types.index(sample_type) % len(_DEFAULT_SAMPLE_TYPE_PALETTE)
    return _DEFAULT_SAMPLE_TYPE_PALETTE[idx]


def sample_fraction_threshold(group_id, fraction):
    """Round-up sample-count threshold for a given fraction of a group's sample size."""
    return ceil(len(GROUPS[group_id]) * fraction)


# ---------------------------------------------------------------------------
# Helper: booleans from config (default True)
# ---------------------------------------------------------------------------
def flag(key):
    return config.get(key, True)


# ---------------------------------------------------------------------------
# Helper: per-sample, per-rule thread count. Rules use:
#   threads: lambda wc: _rule_threads(wc, "rule_key")
# Resolution order: run_config.yaml sample entry (e.g. longcallr_threads: 4),
# else config/config.yaml's global "threads:" default.
# (detect_ase_outliers and get_junction_counts are always single-threaded
# and don't use this helper.)
# ---------------------------------------------------------------------------
def _rule_threads(wc, rule_key):
    return int(SAMPLES[wc.sample].get(f"{rule_key}_threads", config["threads"]))


# ---------------------------------------------------------------------------
# Per-group / per-bed-panel thread & memory overrides, for rules that run
# once per (bed, sample_type) group or once per bed panel rather than per
# sample. Set via a top-level `groups:` block in the run config YAML:
#
#   groups:
#     IEI422_gene_symbols_fibroblasts:            # group-level: bed_id_sample_type
#       build_group_junction_matrix_threads: 1
#       build_group_junction_matrix_mem_gb: 40
#       identify_cohort_junction_outliers_threads: 16
#       identify_cohort_junction_outliers_mem_gb: 96
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
            # Requesting the chain's final output pulls the rest of the
            # junction_analysis.smk chain along with it.
            for t in sample_tissues(s):
                outs.append(
                    f"{od}/junction_analysis/gtex_{t}/{s}_gtex_{t}_outlier_junctions.tsv"
                )

    if flag("merge_hits"):
        # Requesting final_merge's output pulls the whole merge_hits.smk
        # chain along with it (see docs/rules/6_merge_hits.md).
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
        # Requesting validate_sample_types' output pulls build_group_junction_matrix along with it.
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
include: "rules/1_call_variants.smk"
include: "rules/2_compile_variants.smk"
include: "rules/3_phase_reads.smk"
include: "rules/4_ase_analysis.smk"
include: "rules/5_junction_analysis.smk"
include: "rules/6_merge_hits.smk"
include: "rules/7_cohort_junction_analysis.smk"
include: "rules/8_validate_sample_types.smk"
