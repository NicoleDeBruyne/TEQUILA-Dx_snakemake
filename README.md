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
6. **Merges hits across samples** (`merge_hits`) — for each (BED panel, sample type) group, filters and combines the variant, ASE, and junction results into one final `merged_all_hits.tsv` file per panel.
7. **Finds splice-junction outliers vs. cohort** (`cohort_junction_analysis`) — compares each sample's junction usage to **the rest of the cohort's own samples**, instead of GTEx.
8. **Checks cohort-wide QC** (`cohort_qc`) — for each BED panel: mapping/on-target rates, read-length distributions, a "full-length ratio" per gene, and a sample-identity check
9. **Quantifies gene expression** (`quantify_genes`) — for each (BED panel, sample type) group: approximate relative gene expression by read count, by peak coverage, by splice-site-sharing assignment, and by isoform-aware quantification (AMALGAM).


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
│   ├── 6_merge_hits.smk
│   ├── 7_cohort_junction_analysis.smk
│   ├── 8_cohort_qc.smk                # On-target rates, read attributes, full-length ratio, sample-identity check
│   └── 9_quantify_genes.smk           # Gene expression matrices (by count, coverage, assignment, and AMALGAM)
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
```

Samples that share the same `bed` panel and `sample_type` are grouped together for cohort-level analyses.

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
  samples/{sample}/output/...       -- each sample's own results (variant_calling, phased_reads, ase_analysis, junction_analysis)
  samples/{sample}/logs/...         -- each sample's own logs
  cohort/{bed_id}/
    output/
      validate_sample_types/        -- sample-identity check plots (covers all sample types on this panel)
      cohort_qc/                    -- QC plots for on-target rate / read attributes (covers all sample types on this panel)
      merged_all_hits.tsv           -- the final diagnostic result for this BED panel
      sample_types/{sample_type}/
        output/...                  -- results for this group of samples: merged_variant_calling,
                                        merged_ase_analysis, merged_junction_analysis, merged_hits,
                                        cohort_junction_analysis, gene_quantification
        logs/...                    -- logs for this group
    logs/...                        -- logs for the whole BED panel
```

Every output folder of results has a matching `logs/` folder next to it, one level up:
this applies at the sample level (`samples/{sample}/`), the BED-panel level (`cohort/{bed_id}/`),
and the sample-type group level (`cohort/{bed_id}/output/sample_types/{sample_type}/`).