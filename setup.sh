#!/usr/bin/env bash
#
# setup_resources.sh
#
# Builds this pipeline's conda environments and populates snakemake/resources/
# with everything config/config.yaml expects by default. See docs/setup.md
# for what each step does and why.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # snakemake/
RESOURCES_DIR="$SCRIPT_DIR/resources"

mkdir -p "$RESOURCES_DIR"

log()  { echo -e "\n=== $* ==="; }
skip() { echo "  [skip] $1 already exists"; }
error() {
    # Reports a step failure without aborting the whole script -- set -e
    # still applies within each step, but the rest of the script proceeds.
    echo -e "\n  ERROR: $* -- see above (and any log file mentioned) for" >&2
    echo -e "  details. Continuing to the next step.\n" >&2
}

# Prefer mamba over conda for this script's own env-creation calls (faster
# solver, same envs/channels/behavior). Doesn't affect CADD-scripts' own
# install.sh, which has its own conda/mamba selection logic.
if command -v mamba >/dev/null 2>&1; then
    CONDA_BIN=mamba
else
    CONDA_BIN=conda
fi
echo "Using '$CONDA_BIN' for this script's own environment creation."

# Works around a libmamba-solver bug on some systems where strict channel
# priority causes unrelated packages to be reported as mutually
# unsatisfiable. See docs/setup.md.
export CONDA_CHANNEL_PRIORITY=flexible

fetch() {
    # fetch <url> <target_path>
    local url="$1" target="$2"
    if [ -s "$target" ]; then
        skip "$target"
        return
    fi
    echo "  Downloading $target"
    mkdir -p "$(dirname "$target")"
    # --progress=dot:giga: periodic lines instead of a carriage-return bar,
    # readable in a log file and safe with concurrent fetch() calls.
    wget -q --progress=dot:giga -O "$target.partial" "$url"
    mv "$target.partial" "$target"
}
export -f fetch skip

##############################################################################
log "Conda environments"
##############################################################################
# Builds both conda environments this pipeline needs -- see docs/setup.md
# for what each one is for and why they're separate. Both idempotent
# (skipped if already present), same as every other step in this script.
CONDA_ENV_LOG="$RESOURCES_DIR/.setup_logs/conda_envs.log"
mkdir -p "$RESOURCES_DIR/.setup_logs"
echo "  Progress: tail -f $CONDA_ENV_LOG"
(
    env_create_or_hint() {
        # Wraps `$CONDA_BIN env create ...`; on failure, prints a pointer
        # to the channel-priority workaround above.
        if ! "$CONDA_BIN" env create "$@"; then
            cat <<'EOF'

  Environment creation failed. If the error above mentions
  "SOLVER_RULE_STRICT_REPO_PRIORITY" or reports unrelated packages
  (htslib/pysam/snakemake/...) as mutually unsatisfiable, this is a known
  libmamba-solver bug under strict channel priority. This script already
  sets CONDA_CHANNEL_PRIORITY=flexible to work around it, but if that
  didn't take effect for some reason, try setting it globally and re-run:
      conda config --set channel_priority flexible
EOF
            return 1
        fi
    }

    CONDA_ENV_DIR="$SCRIPT_DIR/envs/conda_env"
    if [ -x "$CONDA_ENV_DIR/bin/python" ]; then
        skip "$CONDA_ENV_DIR (main env)"
    else
        echo "  Creating the main environment at $CONDA_ENV_DIR from environment.yaml..."
        env_create_or_hint -p "$CONDA_ENV_DIR" -f "$SCRIPT_DIR/environment.yaml"
    fi

    COMPILE_VARIANTS_ENV_DIR="$SCRIPT_DIR/envs/conda_env_compile_variants"
    if [ -x "$COMPILE_VARIANTS_ENV_DIR/bin/snakemake" ]; then
        skip "$COMPILE_VARIANTS_ENV_DIR (compile_variants env)"
    else
        echo "  Creating the compile_variants environment at $COMPILE_VARIANTS_ENV_DIR"
        echo "  from environment_compile_variants.yaml..."
        env_create_or_hint -p "$COMPILE_VARIANTS_ENV_DIR" -f "$SCRIPT_DIR/environment_compile_variants.yaml"
    fi
) > "$CONDA_ENV_LOG" 2>&1 || error "Conda environments step failed (see $CONDA_ENV_LOG)"

##############################################################################
log "Genome (GENCODE release 44)"
##############################################################################
(
    GENOME_FA="$RESOURCES_DIR/gencode_data/GRCh38.primary_assembly.genome.fa"
    if [ -s "$GENOME_FA" ]; then
        skip "$GENOME_FA"
    else
        fetch "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz" \
            "$GENOME_FA.gz"
        echo "  Decompressing..."
        gunzip "$GENOME_FA.gz"
    fi
    if [ ! -s "$GENOME_FA.fai" ]; then
        if command -v samtools >/dev/null 2>&1; then
            echo "  Indexing genome (samtools faidx)..."
            samtools faidx "$GENOME_FA"
        else
            echo "  WARNING: samtools not found on PATH -- run 'samtools faidx $GENOME_FA' manually."
        fi
    fi
) || error "Genome step failed"

##############################################################################
log "Annotation (GENCODE v44 GTF)"
##############################################################################
(
    ANNOTATION_GTF="$RESOURCES_DIR/gencode_data/gencode.v44.primary_assembly.annotation.gtf"
    if [ -s "$ANNOTATION_GTF" ]; then
        skip "$ANNOTATION_GTF"
    else
        fetch "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" \
            "$ANNOTATION_GTF.gz"
        echo "  Decompressing..."
        gunzip "$ANNOTATION_GTF.gz"
    fi
) || error "Annotation step failed"

##############################################################################
log "NanoTS"
##############################################################################
(
    NANOTS_DIR="$RESOURCES_DIR/NanoTS"
    if [ -d "$NANOTS_DIR/.git" ] || [ -d "$NANOTS_DIR/model" ]; then
        skip "$NANOTS_DIR"
    else
        echo "  Cloning NanoTS..."
        git clone https://github.com/Xinglab/NanoTS.git "$NANOTS_DIR" \
            || git clone git@github.com:Xinglab/NanoTS.git "$NANOTS_DIR"
    fi
) || error "NanoTS step failed"

##############################################################################
log "longcallR v1.12.0"
##############################################################################
# Installed via bioconda as part of the main conda_env -- this just checks
# it's really there after the "Conda environments" step above.
(
    CONDA_ENV_DIR="$SCRIPT_DIR/envs/conda_env"
    LONGCALLR_BIN="$CONDA_ENV_DIR/bin/longcallR"
    if [ -x "$LONGCALLR_BIN" ]; then
        skip "$LONGCALLR_BIN"
    else
        cat <<EOF
  MISSING: $LONGCALLR_BIN

  longcallR should have been installed via bioconda as part of the main
  conda_env build (see environment.yaml) in the "Conda environments" step
  above. If that step failed or hasn't been run yet, re-run this script.
  If conda_env built successfully but this binary still isn't here, check
  "$CONDA_ENV_DIR/bin/" directly for what conda actually installed it as --
  update config.yaml's longcallr_bin to match if the name differs from
  "longcallR".
EOF
    fi
) || error "longcallR check failed"

CADD_DIR="$RESOURCES_DIR/CADD-scripts-1.7.1"
CLINVAR_DIR="$RESOURCES_DIR/clinvar_data"
mkdir -p "$RESOURCES_DIR/.setup_logs"

##############################################################################
log "GTEx junction count matrices (v11, filtered per sample_type)"
##############################################################################
GTEX_LOG="$RESOURCES_DIR/.setup_logs/gtex.log"
echo "  Progress: tail -f $GTEX_LOG"
(
# Edit this map to add/change which GTEx SMTSD tissue label(s) each
# sample_type corresponds to. Multiple SMTSD values are comma-separated.
declare -A GTEX_TISSUE_MAP=(
    ["brain"]="Brain - Amygdala,Brain - Anterior cingulate cortex (BA24),Brain - Caudate (basal ganglia),Brain - Cerebellar Hemisphere,Brain - Cerebellum,Brain - Cortex,Brain - Frontal Cortex (BA9),Brain - Hippocampus,Brain - Hypothalamus,Brain - Nucleus accumbens (basal ganglia),Brain - Putamen (basal ganglia),Brain - Spinal cord (cervical c-1),Brain - Substantia Nigra"
    ["fibroblasts"]="Cells - Cultured fibroblasts"
    ["wholeblood"]="Whole Blood"
    ["lymphocytes"]="Cells - EBV-transformed lymphocytes"
)

GTEX_DIR="$RESOURCES_DIR/gtex_data"
GTEX_RAW_DIR="$GTEX_DIR/raw"
mkdir -p "$GTEX_RAW_DIR"

JUNCTIONS_GZ="$GTEX_RAW_DIR/GTEx_Analysis_2025-08-22_v11_STARv2.7.11b_junctions.gct.gz"
SAMPLE_ATTRIBUTES="$GTEX_RAW_DIR/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"

fetch "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_STARv2.7.11b_junctions.gct.gz" "$JUNCTIONS_GZ"
fetch "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt" "$SAMPLE_ATTRIBUTES"

# Phase 1: identify sample IDs per tissue (cheap -- sample attributes file
# is tiny compared to the junctions matrix). Builds one combined
# sample->tissue map so phase 2 reads the huge junctions file only once.
COMBINED_MAP="$GTEX_RAW_DIR/combined_sample_tissue_map.txt"
: > "$COMBINED_MAP"
tissues_to_filter=()
for tissue in "${!GTEX_TISSUE_MAP[@]}"; do
    smtsd_list="${GTEX_TISSUE_MAP[$tissue]}"
    outfile="$GTEX_DIR/gtex_${tissue}_jxn_counts.txt"

    if [ -s "$outfile" ]; then
        skip "$outfile"
        continue
    fi

    echo "  Identifying $tissue sample IDs (SMTSD: $smtsd_list)..."
    sampids_file="$GTEX_RAW_DIR/${tissue}_SAMPIDs.txt"
    : > "$sampids_file"
    echo "$smtsd_list" | tr ',' '\n' | while read -r smtsd; do
        awk -F'\t' -v smtsd="$smtsd" '
            NR==1 {
                for (i=1; i<=NF; i++) {
                    if ($i == "SAMPID")  sampid_col  = i
                    if ($i == "SMTSD")   smtsd_col   = i
                    if ($i == "SMAFRZE") smafrze_col = i
                }
                if (!sampid_col || !smtsd_col || !smafrze_col) {
                    print "ERROR: could not find SAMPID/SMTSD/SMAFRZE columns in sample attributes header" > "/dev/stderr"
                    exit 1
                }
                next
            }
            $smtsd_col == smtsd && $smafrze_col == "RNASEQ" { print $sampid_col }
        ' "$SAMPLE_ATTRIBUTES"
    done | sort -u > "$sampids_file"

    n_samples=$(wc -l < "$sampids_file")
    if [ "$n_samples" -eq 0 ]; then
        echo "  WARNING: no samples found for $tissue (SMTSD: $smtsd_list) -- skipping." >&2
        continue
    fi
    echo "  Found $n_samples samples for $tissue."

    awk -v t="$tissue" '{print $0 "\t" t}' "$sampids_file" >> "$COMBINED_MAP"
    tissues_to_filter+=("$tissue")
done

if [ "${#tissues_to_filter[@]}" -eq 0 ]; then
    echo "  Nothing to filter -- all tissue matrices already present."
else
    echo "  Filtering junction count matrix for ${#tissues_to_filter[@]} tissue(s) in a single pass: ${tissues_to_filter[*]}..."

    # Phase 2: one streaming pass over the (decompressed) junctions file,
    # streaming the index column + each tissue's columns to that tissue's
    # temp output file. O(1) memory regardless of how many tissues are filtered.
    zcat "$JUNCTIONS_GZ" | tail -n +3 | awk -F'\t' -v OFS='\t' -v idx_name='Name' -v outdir="$GTEX_RAW_DIR" '
        NR==FNR { tissue_of[$1] = $2; next }
        FNR==1 {
            idx_col = 0
            ntissues = 0
            for (i = 1; i <= NF; i++) {
                if ($i == idx_name) idx_col = i
                if ($i in tissue_of) {
                    t = tissue_of[$i]
                    if (!(t in tfile)) {
                        ntissues++
                        tissue_list[ntissues] = t
                        tfile[t] = outdir "/gtex_" t "_jxn_counts.txt.tmp"
                    }
                    n[t]++
                    cols[t, n[t]] = i
                }
            }
            if (!idx_col) {
                print "ERROR: index column \"" idx_name "\" not found in junctions file header" > "/dev/stderr"
                exit 1
            }
            for (k = 1; k <= ntissues; k++) {
                t = tissue_list[k]
                line = $idx_col
                for (j = 1; j <= n[t]; j++) line = line OFS $(cols[t, j])
                print line > tfile[t]
            }
            next
        }
        {
            for (k = 1; k <= ntissues; k++) {
                t = tissue_list[k]
                line = $idx_col
                for (j = 1; j <= n[t]; j++) line = line OFS $(cols[t, j])
                print line > tfile[t]
            }
        }
    ' "$COMBINED_MAP" -

    for tissue in "${tissues_to_filter[@]}"; do
        outfile="$GTEX_DIR/gtex_${tissue}_jxn_counts.txt"
        tmpfile="$GTEX_RAW_DIR/gtex_${tissue}_jxn_counts.txt.tmp"

        # The GTEx junction count matrix contains duplicate rows; dedupe,
        # then check for leftover duplicate junction IDs.
        num_dup_rows=$(sort "$tmpfile" | uniq -d | wc -l)
        if [ "$num_dup_rows" -gt 0 ]; then
            echo "  Removing $num_dup_rows duplicate rows from $tissue junction count matrix..."
            head -n 1 "$tmpfile" > "$outfile"
            tail -n +2 "$tmpfile" | sort -u >> "$outfile"
        else
            mv "$tmpfile" "$outfile"
        fi
        num_dup_idx=$(awk -F'\t' 'NR>1{print $1}' "$outfile" | sort | uniq -d | wc -l)
        if [ "$num_dup_idx" -gt 0 ]; then
            echo "  WARNING: $tissue still has $num_dup_idx duplicate junction IDs after dedup." >&2
        fi
        rm -f "$tmpfile"
    done
fi
echo "  Done. (Raw downloads kept in $GTEX_RAW_DIR for re-filtering later, e.g. if"
echo "  you add a sample_type to GTEX_TISSUE_MAP -- safe to delete if you don't"
echo "  expect to add more and want to reclaim the disk space.)"
) > "$GTEX_LOG" 2>&1 || error "GTEx step failed (see $GTEX_LOG)"

##############################################################################
log "gnomAD"
##############################################################################
# Downloaded locally by default -- see docs/setup.md for why (remote
# querying needs reliable outbound HTTPS from compute nodes, which has been
# unreliable in practice).
#
#   gnomAD v4.1 genomes sites, chroms below + chrM from v3.1 -> ~300GB+
#
# Fetched GNOMAD_PARALLEL files at a time (default 6, override via env var).
# Safe to re-run/resume. To query the remote HTTPS URL instead, set
# SKIP_GNOMAD=y and point gnomad_base/gnomad_mito_vcf at their https:// values.
#
# Keep in sync with config.yaml's gnomad_chroms list.
GNOMAD_LOG="$RESOURCES_DIR/.setup_logs/gnomad.log"
echo "  Progress: tail -f $GNOMAD_LOG"
GNOMAD_CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
               chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
GNOMAD_DIR="$RESOURCES_DIR/gnomad_data"
GNOMAD_PARALLEL="${GNOMAD_PARALLEL:-6}"

(
    if [ "${SKIP_GNOMAD:-n}" = "y" ]; then
        echo "  SKIP_GNOMAD=y set -- skipping. Remember to point gnomad_base/"
        echo "  gnomad_mito_vcf back at their https:// values in config.yaml."
    else
        echo "  Downloading gnomAD v4.1 (${#GNOMAD_CHROMS[@]} autosome/X/Y VCFs) + chrM"
        echo "  (v3.1), $GNOMAD_PARALLEL files at a time..."
        {
            for chrom in "${GNOMAD_CHROMS[@]}"; do
                fname="gnomad.genomes.v4.1.sites.${chrom}.vcf.bgz"
                printf '%s\t%s\n' "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/$fname" "$GNOMAD_DIR/$fname"
                printf '%s\t%s\n' "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/$fname.tbi" "$GNOMAD_DIR/$fname.tbi"
            done
            printf '%s\t%s\n' "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz" "$GNOMAD_DIR/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
            printf '%s\t%s\n' "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi" "$GNOMAD_DIR/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi"
        } | xargs -P "$GNOMAD_PARALLEL" -n 2 bash -c 'fetch "$1" "$2"' _
        echo "  gnomAD done."
    fi
) > "$GNOMAD_LOG" 2>&1 || error "gnomAD step failed (see $GNOMAD_LOG)"

##############################################################################
log "ClinVar"
##############################################################################
# Downloaded locally by default, same reasoning as gnomAD above. Set
# SKIP_CLINVAR=y to skip and keep clinvar_vcf's https:// value in config.yaml.
#
#   ClinVar (GRCh38 VCF) -> ~200MB
CLINVAR_LOG="$RESOURCES_DIR/.setup_logs/clinvar.log"
echo "  Progress: tail -f $CLINVAR_LOG"
(
    if [ "${SKIP_CLINVAR:-n}" = "y" ]; then
        echo "  SKIP_CLINVAR=y set -- skipping. Remember to point clinvar_vcf back"
        echo "  at its https:// value in config.yaml."
    else
        echo "  Downloading ClinVar..."
        fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
            "$CLINVAR_DIR/clinvar.vcf.gz"
        fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" \
            "$CLINVAR_DIR/clinvar.vcf.gz.tbi"
        echo "  ClinVar done."
    fi
) > "$CLINVAR_LOG" 2>&1 || error "ClinVar step failed (see $CLINVAR_LOG)"

##############################################################################
log "CADD-scripts v1.7.1"
##############################################################################
# CADD.sh -m needs a `snakemake` >=8.25.2 binary on PATH (this pipeline's
# main env runs Snakemake 7.x). A shim directory is built below containing
# only a symlink to conda_env_compile_variants's snakemake, kept off the
# rest of that env's PATH to avoid leaking its `perl` into CADD.sh's
# per-rule conda environments. See docs/setup.md for the full explanation,
# and why config.yaml's cadd_script should point at the generated
# CADD_wrapper.sh rather than CADD.sh directly.
COMPILE_VARIANTS_ENV_DIR="$SCRIPT_DIR/envs/conda_env_compile_variants"
CADD_SNAKEMAKE_SHIM_DIR="$RESOURCES_DIR/.cadd_snakemake_shim"
mkdir -p "$CADD_SNAKEMAKE_SHIM_DIR"
ln -sf "$COMPILE_VARIANTS_ENV_DIR/bin/snakemake" "$CADD_SNAKEMAKE_SHIM_DIR/snakemake"
CADD_LOG="$RESOURCES_DIR/.setup_logs/cadd_install.log"
echo "  Progress: tail -f $CADD_LOG"
(
    export PATH="$CADD_SNAKEMAKE_SHIM_DIR:$PATH"

    # Same libmamba-solver workaround as the envs/ build above -- also
    # needed for CADD's own internal conda env builds.
    export CONDA_CHANNEL_PRIORITY=flexible

    # Count complete conda envs: dir + .yaml + .env_setup_done all present.
    count_complete_cadd_envs() {
        n=0
        for envdir in "$CADD_DIR"/envs/conda/*/; do
            [ -d "$envdir" ] || continue
            h="${envdir%/}"
            [ -f "${h}.yaml" ] && [ -f "${h}.env_setup_done" ] && n=$((n + 1))
        done
        echo "$n"
    }
    dir_nonempty() { [ -d "$1" ] && [ -n "$(ls -A "$1" 2>/dev/null)" ]; }

    if [ ! -f "$CADD_DIR/CADD.sh" ]; then
        echo "  Cloning CADD-scripts v1.7.1..."
        git clone --branch v1.7.1 --depth 1 https://github.com/kircherlab/CADD-scripts.git "$CADD_DIR"
    else
        skip "$CADD_DIR/CADD.sh (repo already cloned)"
    fi

    # CADD-scripts' own regulatory-sequence env pins an old TensorFlow
    # without a compatible protobuf pin, which breaks at runtime on a
    # fresh solve. Patched here idempotently -- see docs/setup.md for the
    # full explanation, including why this is keyed off the known filename
    # rather than a content grep (envs/mmsplice.yml also mentions tensorflow).
    REGSEQ_YML="$CADD_DIR/envs/regulatorySequence.yml"
    if [ ! -f "$REGSEQ_YML" ]; then
        echo "  WARNING: could not find CADD-scripts' regulatory-sequence env yaml"
        echo "  (expected $REGSEQ_YML) -- skipping the protobuf pin below. If"
        echo "  annotate_regseq later fails with a protobuf error (\"Descriptors"
        echo "  cannot be created directly\" or similar), pin protobuf<3.20 in"
        echo "  that env's yaml manually and rebuild it."
    elif grep -q "protobuf" "$REGSEQ_YML"; then
        skip "protobuf pin in $REGSEQ_YML (already present)"
    else
        echo "  Pinning protobuf<3.20 in $REGSEQ_YML (works around a TensorFlow/"
        echo "  protobuf version mismatch in CADD-scripts' own env spec -- see"
        echo "  comment above)."
        cp "$REGSEQ_YML" "$REGSEQ_YML.orig"
        if python3 -c "import yaml" 2>/dev/null; then
            python3 - "$REGSEQ_YML" <<'PYEOF'
import sys, yaml
path = sys.argv[1]
with open(path) as f:
    spec = yaml.safe_load(f)
deps = spec.setdefault('dependencies', [])
# Prefer an existing pip: sub-list if present, else add a plain conda dependency.
pip_list = None
for d in deps:
    if isinstance(d, dict) and 'pip' in d:
        pip_list = d['pip']
        break
if pip_list is not None:
    pip_list.append('protobuf<3.20')
else:
    deps.append('protobuf<3.20')
with open(path, 'w') as f:
    yaml.safe_dump(spec, f, default_flow_style=False, sort_keys=False)
PYEOF
        else
            echo "  WARNING: python3's yaml module isn't available to safely edit"
            echo "  $REGSEQ_YML -- falling back to a plain text append (assumes"
            echo "  standard 2-space YAML list indentation; verify $REGSEQ_YML"
            echo "  by hand if CADD-scripts' own formatting differs)."
            printf '\n# Pinned by setup_resources.sh -- works around a TensorFlow/protobuf\n# version mismatch (see setup_resources.sh CADD-scripts section comment).\n  - protobuf<3.20\n' >> "$REGSEQ_YML"
        fi
    fi

    if dir_nonempty "$CADD_DIR/data/annotations/GRCh38_v1.7"; then
        skip "CADD annotations (data/annotations/GRCh38_v1.7 already present)"
        ann_answer=n
    else
        ann_answer=y
    fi

    if dir_nonempty "$CADD_DIR/data/prescored/GRCh38_v1.7/no_anno"; then
        skip "CADD prescored variants (data/prescored/GRCh38_v1.7/no_anno already present)"
        pre_answer=n
    else
        pre_answer=y
    fi

    echo "  Running CADD's own installer, answering its prompts automatically"
    echo "  (safe to rerun -- annotations/prescored are skipped above if already"
    echo "  present on disk, since install.sh has no such check of its own;"
    echo "  conda env creation is always requested since Snakemake's conda"
    echo "  integration skips any env that's already complete):"
    echo "    1. Install conda/mamba environments?  -> y (idempotent -- Snakemake"
    echo "       skips any env under $CADD_DIR/envs/conda already complete)"
    echo "    2. Install CADD for GRCh37/hg19?       -> n (not used by this pipeline;"
    echo "       NOTE this prompt defaults to YES if left unanswered -- a 261GB+ download)"
    echo "    3. Install CADD for GRCh38/hg38?       -> y (must always be y -- the"
    echo "       annotations/prescored downloads below are nested under this answer)"
    echo "    4. Load annotations?                   -> \$ann_answer (~336GB if needed)"
    echo "    5. Load prescored variants?             -> \$pre_answer (~81GB if needed;"
    echo "       speeds up scoring -- the prescored SNV file covers all possible"
    echo "       genome-wide SNVs, so any SNV this pipeline queries hits the cache"
    echo "       instead of being computed from scratch)"
    echo "    6-8. (only asked if #5 is y) with-anno / without-anno / InDels"
    echo "       prescored -> n / \$pre_answer / \$pre_answer (skipped entirely if #5 is n --"
    echo "       install.sh does NOT ask these when prescored=n, so sending answers for"
    echo "       them anyway would misalign with the final 'Ready to continue?' prompt"
    echo "       and cancel the install)"
    echo "    9. Ready to continue?                   -> y"
    if [ "$pre_answer" = "y" ]; then
        # Prompt 5=y -> install.sh also asks 6, 7, 8 -> 9 prompts total.
        printf 'y\nn\ny\n%s\ny\nn\ny\ny\ny\n' "$ann_answer" | \
            ( cd "$CADD_DIR" && bash install.sh )
    else
        # Prompt 5=n -> install.sh skips 6, 7, 8 entirely -> only 6 prompts total.
        printf 'y\nn\ny\n%s\nn\ny\n' "$ann_answer" | \
            ( cd "$CADD_DIR" && bash install.sh )
    fi
    echo "  CADD-scripts install.sh step done."
    chmod +x "$CADD_DIR/CADD.sh"

    # If the regulatory-sequence env was already built before the yaml pin
    # above existed, the patch doesn't retroactively fix it (Snakemake only
    # rebuilds an env if its yaml hash changes) -- so also patch any
    # already-built env directly. Matched by content against $REGSEQ_YML
    # (or its .orig backup, if present) rather than by grepping for
    # "tensorflow", to avoid mismatching mmsplice.yml's built copy. See
    # docs/setup.md for details.
    REGSEQ_YML_FOR_MATCH="$REGSEQ_YML"
    [ -f "$REGSEQ_YML.orig" ] && REGSEQ_YML_FOR_MATCH="$REGSEQ_YML.orig"

    BUILT_REGSEQ_YAML=""
    if [ -f "$REGSEQ_YML_FOR_MATCH" ]; then
        for candidate in "$CADD_DIR"/envs/conda/*.yaml; do
            [ -f "$candidate" ] || continue
            if diff -q <(grep -v '^prefix:' "$REGSEQ_YML_FOR_MATCH") <(grep -v '^prefix:' "$candidate") >/dev/null 2>&1; then
                BUILT_REGSEQ_YAML="$candidate"
                break
            fi
        done
    fi
    if [ -n "$BUILT_REGSEQ_YAML" ]; then
        BUILT_REGSEQ_ENV_DIR="${BUILT_REGSEQ_YAML%.yaml}"
        if [ -x "$BUILT_REGSEQ_ENV_DIR/bin/pip" ]; then
            CURRENT_PROTOBUF="$("$BUILT_REGSEQ_ENV_DIR/bin/pip" show protobuf 2>/dev/null | awk '/^Version:/{print $2}')"
            # Only the major.minor matters here (3.20 is the first
            # incompatible release), so a simple string comparison suffices.
            case "$CURRENT_PROTOBUF" in
                ""|3.19.*|3.1[0-8].*|3.[0-9].*|2.*)
                    ;;  # already <3.20 (or pip show failed) -- nothing to do
                *)
                    echo "  Downgrading protobuf ($CURRENT_PROTOBUF -> <3.20) in the already-built"
                    echo "  regulatory-sequence env ($BUILT_REGSEQ_ENV_DIR) -- see the protobuf"
                    echo "  pin comment above for why."
                    "$BUILT_REGSEQ_ENV_DIR/bin/pip" install "protobuf<3.20" --quiet
                    ;;
            esac
        fi
    fi

    # Generate CADD_wrapper.sh -- this, not CADD.sh directly, is what
    # config.yaml's cadd_script should point to. Self-contained (builds its
    # own clean PATH, stripping conda_env_compile_variants's bin/ to avoid
    # leaking its `perl` into CADD.sh's per-rule conda environments) so
    # compile_variants.py can invoke it directly. See docs/setup.md.
    CADD_WRAPPER="$CADD_DIR/CADD_wrapper.sh"
    cat > "$CADD_WRAPPER" <<WRAPPER_EOF
#!/bin/bash
# Auto-generated by setup_resources.sh -- do not edit directly; rerun
# setup_resources.sh to regenerate. See docs/setup.md.
set -euo pipefail
SHIM_DIR="$CADD_SNAKEMAKE_SHIM_DIR"
STRIP_DIR="$COMPILE_VARIANTS_ENV_DIR/bin"
CLEAN_PATH=""
IFS=':' read -ra _PARTS <<< "\$PATH"
for _p in "\${_PARTS[@]}"; do
    if [ "\$_p" != "\$STRIP_DIR" ]; then
        CLEAN_PATH="\${CLEAN_PATH:+\$CLEAN_PATH:}\$_p"
    fi
done
export PATH="\$SHIM_DIR:\$CLEAN_PATH"
exec "$CADD_DIR/CADD.sh" "\$@"
WRAPPER_EOF
    chmod +x "$CADD_WRAPPER"
    echo "  Generated $CADD_WRAPPER (point config.yaml's cadd_script here, not at CADD.sh)."

    # Force a real scoring pass here, once, serially, against CADD-scripts'
    # bundled test VCF -- builds every conda env CADD.sh actually needs
    # (beyond install.sh's own narrower env-build target) and confirms
    # scoring works before any real sample touches it, avoiding a race if
    # multiple concurrent compile_variants.py runs each try to build the
    # same missing env. Goes through the wrapper, not CADD.sh directly, so
    # this exercises exactly what compile_variants.py invokes at runtime.
    # See docs/setup.md.
    n_envs="$(count_complete_cadd_envs)"
    echo "  Running a real scoring pass (via CADD_wrapper.sh) against the bundled"
    echo "  test VCF ($n_envs/5 expected conda envs currently complete) to"
    echo "  force-build the full env set (beyond install.sh's own narrower"
    echo "  env-build target) and confirm scoring works. Already-complete envs"
    echo "  are skipped, so this is cheap when everything's already built, and"
    echo "  slow (conda solves + several envs, several GB each) the first time."
    ( cd "$CADD_DIR" && ./CADD_wrapper.sh -m -c 1 -o /tmp/cadd_setup_test_output.tsv.gz -g GRCh38 test/input.vcf.gz )
    rm -f /tmp/cadd_setup_test_output.tsv.gz
    n_envs="$(count_complete_cadd_envs)"
    echo "  CADD test-scoring pass done -- $n_envs/5 envs now complete."
) > "$CADD_LOG" 2>&1 || error "CADD-scripts step failed (see $CADD_LOG)"

##############################################################################
log "ANNOVAR (requires free registration -- cannot be auto-downloaded)"
##############################################################################
(
    ANNOVAR_DIR="$RESOURCES_DIR/annovar"
    if [ -d "$ANNOVAR_DIR" ] && [ -n "$(ls -A "$ANNOVAR_DIR" 2>/dev/null)" ]; then
        skip "$ANNOVAR_DIR"
    else
        cat <<EOF
  MISSING: $ANNOVAR_DIR

  ANNOVAR requires a free academic registration before download:
    1. Register at https://www.openbioinformatics.org/annovar/annovar_download_form.php
    2. Download the annovar.latest.tar.gz link emailed to you
    3. Extract it so that its contents land directly in:
         $ANNOVAR_DIR
       (i.e. $ANNOVAR_DIR/annotate_variation.pl etc., not a nested subfolder)
    4. Download the annotation databases this pipeline uses (refGene, etc.)
       via ANNOVAR's own annotate_variation.pl -downdb -webfrom annovar ...
EOF
    fi
) || error "ANNOVAR check failed"

##############################################################################
log "SpliceAI precomputed scores (requires a free BaseSpace account -- cannot be auto-downloaded)"
##############################################################################
# Masked (not raw) is the right choice -- matches this pipeline's live
# `spliceai -M 1` invocation. See docs/setup.md for why.
(
    SPLICEAI_DIR="$RESOURCES_DIR/spliceai_data"
    SPLICEAI_SNV="$SPLICEAI_DIR/spliceai_scores.masked.snv.hg38.vcf.gz"
    SPLICEAI_INDEL="$SPLICEAI_DIR/spliceai_scores.masked.indel.hg38.vcf.gz"
    if [ -s "$SPLICEAI_SNV" ] && [ -s "$SPLICEAI_SNV.tbi" ] && \
       [ -s "$SPLICEAI_INDEL" ] && [ -s "$SPLICEAI_INDEL.tbi" ]; then
        skip "$SPLICEAI_DIR"
    else
        cat <<EOF
  MISSING: $SPLICEAI_SNV (+.tbi)
           $SPLICEAI_INDEL (+.tbi)

  SpliceAI's precomputed scores are only distributed via Illumina BaseSpace,
  which requires a free account -- they can't be fetched with a plain
  wget/curl the way gnomAD/ClinVar are above:
    1. Create a free BaseSpace account (if you don't already have one):
         https://basespace.illumina.com
    2. Go to the SpliceAI precomputed scores project:
         https://basespace.illumina.com/s/otSPW8hnhaZR
    3. Download the GRCh38, MASKED SNV and INDEL files:
         spliceai_scores.masked.snv.hg38.vcf.gz (+ .tbi)
         spliceai_scores.masked.indel.hg38.vcf.gz (+ .tbi)
       (NOT the raw.* files, and NOT the hg37/hg19 build -- see the note
       above on why masked is the right choice for this pipeline)
    4. Place all four files directly in:
         $SPLICEAI_DIR
       (i.e. $SPLICEAI_SNV etc., not a nested subfolder)

  These files are free for academic/non-profit use; other use requires a
  commercial license from Illumina, Inc. (same terms as the spliceai
  package/models this pipeline already runs live via run_SpliceAI_chunk).
EOF
    fi
) || error "SpliceAI precomputed scores check failed"

##############################################################################
log "OMIM (bundled with this repo -- see below for how to refresh it)"
##############################################################################
# OMIM.tsv is small enough to ship as part of this repo directly rather
# than being fetched by this script -- this section just checks it's
# present and reminds how to refresh it later. See docs/setup.md.
(
    OMIM_FILE="$RESOURCES_DIR/omim_data/OMIM.tsv"
    if [ -s "$OMIM_FILE" ]; then
        skip "$OMIM_FILE (bundled with this repo)"
    else
        cat <<EOF
  MISSING: $OMIM_FILE

  This is unexpected -- OMIM.tsv is normally bundled with this repo, so if
  it's missing, something didn't copy over correctly when this repo was
  cloned/copied (check whether $RESOURCES_DIR/omim_data/ itself exists
  and what it actually contains).

  If you want to update to a newer OMIM release instead (OMIM's own data is
  periodically updated; this repo's bundled copy is a point-in-time
  snapshot), OMIM data requires an institutional license
  (https://omim.org/downloads). Download an updated copy and replace:
    $OMIM_FILE

  Whatever you put there must be a tab-separated file containing (at least)
  these three columns -- scripts/merge_hits.py reads only these, any others
  are ignored:
    approved_gene_symbol   Gene symbol -- joined against this pipeline's own
                            ANNOVAR-derived gene symbols (ANNOVAR_Gene.refGene),
                            so naming convention/casing needs to match those.
    phenotypes              Associated disease phenotype(s), passed through
                            as-is into the final merged output.
    inheritance_patterns    Associated inheritance pattern(s) (e.g. autosomal
                            recessive), passed through as-is.
EOF
    fi
) || error "OMIM check failed"

##############################################################################
log "Done. Review any MISSING sections above before running the pipeline."
##############################################################################