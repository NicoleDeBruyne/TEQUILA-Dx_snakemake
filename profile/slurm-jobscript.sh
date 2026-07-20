#!/bin/bash
# profile/slurm-jobscript.sh
# Standard Snakemake jobscript template. Snakemake fills in {properties} as
# a JSON comment and {exec_job} with the actual rule command before handing
# this off to slurm-submit.py for sbatch submission.
#
# Every rule in this pipeline uses the exact same conda environment (no rule
# ever uses a different one), so rules/*.smk no longer declare a per-rule
# `conda:` environment for Snakemake to manage via --use-conda. Declaring one
# used to trigger a Snakemake bug: after a job finished successfully,
# Snakemake tries to record provenance by running `conda env export --name
# '<path>'`, which is invalid syntax for a path-based (non-registered) conda
# env and crashes the whole run -- see
# https://github.com/snakemake/snakemake/issues/1674
#
# Instead, we activate the one shared environment directly here, since each
# SLURM job is a separate sbatch submission that doesn't otherwise inherit
# conda activation state. CONDA_ENV_DIR arrives via the environment (set by
# slurm-submit.py/slurm_utils.py's default_conda_env_dir(), passed through
# via `sbatch --export`) rather than computed here from this script's own
# location: SLURM copies submitted scripts into its own spool directory
# before running them, so $0/$BASH_SOURCE inside a running job reflects the
# spool copy, not this file's real path.
if [ -z "$CONDA_ENV_DIR" ]; then
    echo "WARNING: CONDA_ENV_DIR was not set in the job environment --" >&2
    echo "expected slurm-submit.py to export it via sbatch --export. Conda" >&2
    echo "env won't be activated for this job." >&2
elif [ -d "$CONDA_ENV_DIR/bin" ]; then
    export PATH="$CONDA_ENV_DIR/bin:$PATH"
else
    echo "WARNING: $CONDA_ENV_DIR/bin not found -- check conda_env in" >&2
    echo "config.yaml matches where the env actually is." >&2
fi

# properties = {properties}
{exec_job}