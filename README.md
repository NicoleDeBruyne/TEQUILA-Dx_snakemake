# RNA-Dx Snakemake Pipeline

*(repo name: `TEQUILA-Dx_snakemake`)*

A Snakemake pipeline for diagnostic analysis of long-read RNA-seq data: multi-caller variant
calling, read phasing, allele-specific expression (ASE) outlier detection, and splice-junction
outlier analysis, with cross-sample merging into a candidate diagnostic hit list annotated
against OMIM.

## What it does

For each sample (a long-read RNA-seq BAM), the pipeline:

1. **Calls variants** with four independent callers: [longcallR](https://github.com/huangnengCSU/longcallR), [NanoTS](https://github.com/Xinglab/NanoTS), [Clair3-RNA](https://github.com/HKU-BAL/Clair3-RNA), and [DeepVariant](https://github.com/google/deepvariant).
2. **Compiles variants** (`compile_variants`) — merges the four VCFs and annotates with ANNOVAR, gnomAD allele frequency, ClinVar significance, CADD, and SpliceAI.
3. **Phases reads** (`phase_reads`) — builds a merged VCF from all callers and runs WhatsHap-based per-gene phasing, producing per-haplotype BAMs and a gene→BAM mapping file used downstream.
4. **Detects ASE outliers** (`ase_analysis`) — binomial test for allele-specific expression per gene from the phased haplotype counts.
5. **Detects splice-junction outliers**, two ways:
   - `junction_analysis` — per-sample junction usage vs. **GTEx reference tissue** matrices (beta-binomial test).
   - `cohort_junction_analysis` — per-sample junction usage vs. the **rest of the cohort's own bulk BAMs** (rather than GTEx).
6. **Validates sample types** (`validate_sample_types`) — PSI-based check of each sample's tissue identity against GTEx reference tissues, to flag potential sample mislabeling.
7. **Merges hits across samples** (`merge_hits`) — for each (BED panel, sample type) group, filters and merges variant / ASE / junction hits across samples, annotates against a bundled **OMIM** gene→phenotype table, and concatenates everything into a final `merged_all_hits.tsv` per panel.

## Repository structure

```
TEQUILA-Dx_snakemake/
├── Snakefile                          # Multi-sample orchestration: loads run config, builds
│                                       # sample/group wildcards, defines rule `all`
├── config/
│   └── config.yaml                    # Default paths, resource DBs, thresholds, stage on/off flags
├── profile/                            # Snakemake 7 SLURM cluster profile
│   ├── config.yaml                    # --profile settings (jobs, latency-wait, etc.)
│   ├── slurm-submit.py                # SLURM job submission
│   ├── slurm-status.py                # SLURM job status polling
│   ├── slurm-jobscript.sh             # Per-job wrapper (activates conda env, etc.)
│   └── slurm_utils.py
├── rules/                              # One .smk file per pipeline stage (included by Snakefile)
│   ├── 1_call_variants.smk            # 4 independent callers: longcallR, NanoTS, Clair3-RNA, DeepVariant
│   ├── 2_compile_variants.smk
│   ├── 3_phase_reads.smk
│   ├── 4_ase_analysis.smk
│   ├── 5_junction_analysis.smk
│   ├── 6_merge_hits.smk
│   ├── 7_cohort_junction_analysis.smk
│   └── 8_validate_sample_types.smk
├── scripts/                            # Python scripts invoked by the rules above
├── resources/
│   └── omim_data/OMIM.tsv             # Bundled OMIM gene → phenotype/inheritance table
│                                       # (other resources — genome, gnomAD, ClinVar, CADD,
│                                       # GTEx, etc. — are downloaded by setup.sh, not bundled)
├── environment.yaml                    # Main conda env (`RNA-Dx`) — Snakemake 7.x, used by
│                                       # every rule except compile_variants
├── environment_compile_variants.yaml   # Dedicated conda env (`RNA-Dx-compile-variants`) —
│                                       # Snakemake ≥8.25.2, used only by the compile_variants
│                                       # rule (needed for CADD's own internal snakemake call)
└── setup.sh                            # Builds both conda envs and downloads/prepares all
                                        # reference data under resources/
```

## Requirements

- Linux, SLURM cluster (the bundled `profile/` targets SLURM specifically)
- conda or mamba (mamba preferred — `setup.sh` uses it automatically if present)
- Singularity/Apptainer and/or Docker, for the containerized callers (NanoTS, Clair3-RNA, DeepVariant)
- Substantial disk space for reference data — gnomAD alone is ~300GB+; CADD annotations/prescored
  data can add several hundred GB more (see `setup.sh` for what's downloaded and how to skip pieces)
- Free registrations required for two resources that `setup.sh` **cannot** auto-download:
  - **ANNOVAR** (academic registration)
  - **SpliceAI precomputed scores** (Illumina BaseSpace account) — optional; the pipeline runs
    SpliceAI live via the bioconda `spliceai` package by default and only falls back to these
    precomputed files if that fails

## Setup

```bash
git clone <this-repo>
cd TEQUILA-Dx_snakemake
./setup.sh
```

`setup.sh` is idempotent (safe to re-run — it skips anything already present) and will:
- Create the `envs/conda_env` and `envs/conda_env_compile_variants` conda environments
- Download the GENCODE v44 GRCh38 genome + annotation
- Clone NanoTS and confirm longcallR is installed
- Build per-sample-type GTEx junction count matrices (v11)
- Download gnomAD v4.1 genomes + ClinVar
- Clone and install CADD-scripts v1.7.1 (and generate `CADD_wrapper.sh`)
- Check for ANNOVAR and SpliceAI precomputed scores, printing manual setup instructions if missing
- Verify the bundled `resources/omim_data/OMIM.tsv` is present

Check the script's final output for any `MISSING:` sections before running the pipeline, and see
`resources/.setup_logs/` for per-step logs.

Reference-data and environment paths default to relative paths under `resources/` and `envs/`
inside the pipeline directory (see `config/config.yaml`), so the whole folder is self-contained
and relocatable. Point any of them at an absolute path instead if you want to share a copy across
multiple pipeline checkouts. `gnomad_base`, `clinvar_vcf`, and `cadd_script` can each also be set
to the literal value `"remote"` to query the public HTTPS source directly instead of a local copy.

## Configuration

`config/config.yaml` holds pipeline-wide defaults: reference/database paths, stage on/off flags
(`longcallr`, `nanots`, `clair3_rna`, `deepvariant`, `compile_variants`, `phase_reads`,
`ase_analysis`, `junction_analysis`, `cohort_junction_analysis`, `merge_hits`,
`validate_sample_types`), and filtering thresholds (gnomAD AF, CADD, SpliceAI, ASE p-adj,
splice-junction padj/delta-PSI, etc.).

**Per-run sample manifest:** each run additionally needs a YAML file (path passed via
`--config run=<path>`) defining the samples for that run:

```yaml
output_dir: "/path/to/output"   # or pass via --config output_dir=<path>

samples:
  sample1:
    bam: "/path/to/sample1.bam"
    bed: "/path/to/panel.bed"
    tissues: ["fibroblasts", "wholeblood"]   # GTEx reference tissue(s) to compare against
    sample_type: "fibroblasts"               # used for grouping in the merge_hits stage
    # outdir: "/path/to/sample1_output"      # optional: overrides the default
                                              # {output_dir}/samples/sample1/output for this sample
```

Samples sharing the same `bed` panel and `sample_type` are grouped together for cross-sample
merging (`rules/6_merge_hits.smk`) and cohort-level analyses.

## Usage

Dry run:
```bash
conda activate envs/conda_env
snakemake -n --config run=/path/to/run_config.yaml
```

Run on a SLURM cluster via the bundled profile:
```bash
snakemake --profile profile/ --use-conda --config run=/path/to/run_config.yaml
```

Override any `config.yaml` value at the command line, e.g. to disable a stage:
```bash
snakemake --profile profile/ --config run=/path/to/run_config.yaml merge_hits=False
```

> Note: several comments in the Snakefile/config reference a `submit_snakemake.sh` wrapper
> (which would generate the per-run YAML and inject `output_dir` automatically) — that script
> isn't included in this copy of the repo, so run configs currently need to be written by hand
> as shown above.

## Output

Everything lands under the single `output_dir` given in the run config:

```
{output_dir}/
  samples/{sample}/output/...       -- per-sample stage outputs (variant_calling, phased_reads,
                                        ase_analysis, junction_analysis)
  samples/{sample}/logs/...         -- per-sample stage logs
  cohort/{bed_id}/
    output/
      validate_sample_types/        -- sample-type validation plots (across all sample_types on this panel)
      cohort_qc/                    -- on-target-rate / read-attribute QC plots (across all sample_types on this panel)
      merged_all_hits.tsv           -- the final diagnostic output for this BED panel, annotated
                                        with OMIM phenotype and inheritance information
      sample_types/{sample_type}/
        output/...                  -- per-group merged_variant_calling, merged_ase_analysis,
                                        merged_junction_analysis, merged_hits, cohort_junction_analysis,
                                        gene_quantification (the "merged_" prefix distinguishes these
                                        from the per-sample folders of the same base name)
        logs/...                    -- per-group logs
    logs/...                        -- bed-panel-level logs (validate_sample_types, cohort_qc, final merge) --
                                        a sibling of output/ above, not nested inside it
```

Every output-producing directory at every level has a sibling `logs/` directory one level up from it
(samples/{sample}/, cohort/{bed_id}/, and cohort/{bed_id}/output/sample_types/{sample_type}/ each follow
this same output/ + logs/ sibling pattern).

A sample's `outdir` can be overridden per-sample in the run config (see Configuration above);
its `logs/` directory always sits next to whatever that resolves to.

## License

[Add license information.]
