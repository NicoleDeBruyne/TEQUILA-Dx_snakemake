# RNA-Dx Snakemake Pipeline

A Snakemake pipeline for diagnostic analysis of long-read RNA-seq data. It calls variants with multiple tools, 
phases reads, finds allele-specific expression (ASE) outliers, and finds splice-junction outliers. 
Results from all samples are then merged into one candidate diagnostic hit list.

## What it does

For each sample (a long-read RNA-seq BAM file), the pipeline:

1. **Calls variants** with four independent tools: [longcallR](https://github.com/huangnengCSU/longcallR), [NanoTS](https://github.com/Xinglab/NanoTS), [Clair3-RNA](https://github.com/HKU-BAL/Clair3-RNA), and [DeepVariant](https://github.com/google/deepvariant).
2. **Compiles variants** (`compile_variants`) — combines the four VCFs into one and adds annotations from gnomAD, ClinVar,  ANNOVAR, CADD, and SpliceAI.
3. **Phases reads** (`phase_reads`) — Produces one BAM file per haplotype (where possible) and a table mapping each gene to its phased BAM files, used in later steps.
4. **Finds ASE outliers** (`ase_analysis`) — runs a binomial test per gene on the phased haplotype read counts to detect allele-specific expression.
5. **Finds splice-junction outliers vs. GTEx** (`junction_analysis`) — compares each sample's junction usage to **GTEx reference tissue** data, using a beta-binomial test.
6. **Computes per-sample QC** (`sample_qc`) — mapping/on-target read counts, a read-length five-number summary, and a per-gene "full-length ratio", all computed directly from that sample's own BAM (no cohort dependency).
7. **Quantifies gene expression per sample** (`sample_gene_quantification`) — approximate per-sample relative gene expression by read count, by peak coverage, and by splice-site-sharing assignment, plus a per-sample StringTie assembly for isoform-aware quantification (AMALGAM) later.
8. **Finds splice-junction outliers vs. cohort** (`cohort_junction_analysis`) — compares each sample's junction usage to **the rest of the cohort's own samples**, instead of GTEx.
9. **Merges results across samples** (`merge_hits` / `cohort_qc` / `quantify_genes`) — combines every sample's step 1-7 results (plus step 8's cohort-junction results) into cohort-level outputs: a ranked `merged_all_hits.tsv` candidate-hit list per BED panel, cohort-wide QC matrices/plots, and gene-expression matrices (by count, coverage, assignment, and AMALGAM). This is the only stage that reads results from more than one sample at once, which is what keeps rerunning it (e.g. after adding a sample, or configuring a sample alias) cheap -- it never touches a BAM directly.

**Low-expression outlier score:** each of step 9's four gene-expression matrices (by count, coverage, assignment, and AMALGAM) also gets an `..._outlier_zscores.tsv` (full genes x samples matrix) and an `..._outliers.tsv` (just the flagged rows) alongside it -- a robust, nonparametric per-gene z-score (leave-one-out median/MAD in log2 space, each gene's MAD shrunk toward a cohort-wide MAD-vs-expression trend so genes with too few samples don't produce spurious extreme scores) flagging samples with unusually **low** expression of a gene relative to the rest of that (bed, sample_type) group. See `scripts/expression_outliers.py`'s module docstring for the full algorithm and `gene_outlier_*` in `config/config.yaml` for its tunable parameters. The splice-site-assignment method's score is additionally annotated onto `merged_all_hits.tsv` as `gene_expression_zscore`/`gene_expression_outlier` columns -- informational only for now, it does not currently affect a gene's tier/ranking.


## Repository structure

```
TEQUILA-Dx_snakemake/
├── Snakefile 
├── config/
│   └── config.yaml                    # Default paths, resource databases, thresholds, and on/off switches for each stage
├── profile/                           # SLURM cluster profile for Snakemake 7
│   ├── config.yaml 
│   ├── slurm-submit.py                # Submits jobs to SLURM
│   ├── slurm-status.py                # Checks SLURM job status
│   ├── slurm-jobscript.sh             # Wrapper script run for each job (activates conda env, etc.)
│   └── slurm_utils.py
├── rules/   
│   ├── 1_call_variants.smk 
│   ├── 2_compile_variants.smk
│   ├── 3_phase_reads.smk
│   ├── 4_ase_analysis.smk
│   ├── 5_junction_analysis.smk
│   ├── 6_sample_qc.smk                 # Per-sample: on-target rate, read-length summary, full-length ratio
│   ├── 7_sample_gene_quantification.smk # Per-sample: gene expression by count, coverage, assignment; StringTie assembly
│   ├── 8_cohort_junction_analysis.smk  # Splice-junction outliers vs. the rest of the cohort
│   └── 9_merge_results.smk             # Merges everything above into cohort-wide hit rankings, QC, and gene-expression matrices
├── scripts/                            # Python scripts called by the rules above
├── resources/
│   └── omim_data/OMIM.tsv              # Included OMIM gene -> phenotype/inheritance table
│                                       # (other reference data is downloaded by setup.sh)
├── environment.yaml                    # Main conda environment (`RNA-Dx`), Snakemake 7.x,
├── environment_compile_variants.yaml   # Separate conda environment (`RNA-Dx-compile-variants`),
│                                       # Snakemake >= 8.25.2, used only by the compile_variants
│                                       # rule (because CADD runs Snakemake 8.x internally)
└── setup.sh                            # Builds both conda environments and downloads/prepares all reference data into resources/
```

## Requirements

- Linux, SLURM cluster (the included `profile/` is written for SLURM specifically)
- conda or mamba (mamba is preferred — `setup.sh` uses it automatically if it's installed)
- Singularity/Apptainer and/or Docker, for the containerized tools (NanoTS, Clair3-RNA, DeepVariant)
- A lot of disk space for reference data — gnomAD 596GB, CADD 462GB, SpliceAI 91GB
  data can add several hundred GB more (see `setup.sh` for what gets downloaded and how to skip parts of it)
- Free registration is required for two resources that `setup.sh` **cannot** download automatically:
  - **ANNOVAR** (requires academic registration)
  - **SpliceAI precomputed scores** (requires an Illumina BaseSpace account)

## Setup

```bash
git clone <this-repo>
cd TEQUILA-Dx_snakemake
./setup.sh
```

`setup.sh` can be run more than once safely — it skips anything that's already set up. It will:
- Create the `envs/conda_env` and `envs/conda_env_compile_variants` conda environments
- Download the GENCODE v44 GRCh38 genome and annotation
- Clone NanoTS and check that longcallR is installed via conda
- Build per-sample-type GTEx junction count matrices (v11)
- Download gnomAD v4.1 genomes and ClinVar variant annotations
- Clone and install CADD-scripts v1.7.1 (and generate `CADD_wrapper.sh`)
- Check whether ANNOVAR and the SpliceAI precomputed scores are present, and print manual setup instructions if not
- Check that the included `resources/omim_data/OMIM.tsv` file is present

By default, reference-data and environment paths are relative, pointing at `resources/` and `envs/`
You can set any of these paths to an absolute path instead. `gnomad_base`,
`clinvar_vcf`, and `cadd_script` can each also be set to the literal value `"remote"`, which tells
the pipeline to query the public HTTPS source directly instead of using a local copy.

## Configuration

`config/config.yaml` holds the pipeline's global defaults: reference/database paths, on/off switches for each stage, and filtering thresholds

**Per-run sample list:** each run also needs its own YAML file listing the samples for that run
(its path is passed via `--config run=<path>`):

```yaml
output_dir: "/path/to/output"   # or pass via --config output_dir=<path>

samples:
  sample1:
    bam: "/path/to/sample1.bam"
    bed: "/path/to/panel.bed"
    tissues: ["fibroblasts", "wholeblood"]   # GTEx reference tissue(s) to compare this sample against
    sample_type: "fibroblasts" 
    # outdir: "/path/to/sample1_output"      # optional: use this instead of the default {output_dir}/samples/sample1/output for this sample
    # alias: "PT01"                          # optional: de-identified label to substitute for this sample's real ID in _alias outputs (see below)
```

Samples that share the same `bed` panel and `sample_type` are grouped together for cohort-level analyses.

**Sample aliases:** any sample can be given an optional `alias:` in the run config. For every
cohort-level output (steps 6-9) that shows real sample IDs -- as a TSV column/index or as a plot
label -- the pipeline also writes a companion `..._alias` file (e.g. `merged_all_hits_simplified.tsv`
→ `merged_all_hits_simplified_alias.tsv`) with those IDs replaced by their alias. Samples with no
`alias` configured fall back to their own real ID in the alias file. An `_alias` file is only
produced/requested for a given BED panel or sample-type group if at least one of its samples has
an alias configured; outputs that never show a sample identifier (e.g. `annotated.gtf.gz`, or the
per-sample-unlabeled boxplot/upset plots in `merged_hits/`) don't get an alias companion at all.

## Usage

Dry run (shows what would happen without actually running anything):
```bash
conda activate envs/conda_env
snakemake -n --config run=/path/to/run_config.yaml
```

Run on a SLURM cluster using the included profile:
```bash
snakemake --profile profile/ --use-conda --config run=/path/to/run_config.yaml
```

Any `config.yaml` value can be overridden on the command line, for example to turn a stage off:
```bash
snakemake --profile profile/ --config run=/path/to/run_config.yaml merge_hits=False
```

## Output

Everything is written under the single `output_dir` set in the run config:

```
{output_dir}/
  samples/{sample}/output/...       -- each sample's own results:
                                        variant_calling, phased_reads, ase_analysis, junction_analysis,
                                        qc/                          -- per-sample QC (on-target rate, read-length
                                                                        summary, full-length ratio)
                                        gene_quantification/         -- per-sample gene expression (count, coverage,
                                                                        assignment) + StringTie assembly for AMALGAM
  samples/{sample}/logs/...         -- each sample's own logs
  cohort/{bed_id}/
    output/
      validate_sample_types/        -- sample-identity check plots (covers all sample types on this panel)
      cohort_qc/                    -- cohort-wide QC plots merged from every sample's own qc/ TSVs above
                                        (on-target rate, read lengths, full-length ratio; covers all sample
                                        types on this panel)
      merged_all_hits.tsv           -- the final diagnostic result for this BED panel
      sample_types/{sample_type}/
        output/...                  -- results for this group of samples: merged_variant_calling,
                                        merged_ase_analysis, merged_junction_analysis, merged_hits,
                                        cohort_junction_analysis, gene_quantification (merged from every
                                        sample's own gene_quantification/ TSVs above)
        logs/...                    -- logs for this group
    logs/...                        -- logs for the whole BED panel
```

Every output folder of results has a matching `logs/` folder next to it, one level up:
this applies at the sample level (`samples/{sample}/`), the BED-panel level (`cohort/{bed_id}/`),
and the sample-type group level (`cohort/{bed_id}/output/sample_types/{sample_type}/`).