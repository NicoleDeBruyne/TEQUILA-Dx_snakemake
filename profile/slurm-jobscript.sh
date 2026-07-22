#!/bin/bash
# profile/slurm-jobscript.sh
# Standard Snakemake jobscript template. Snakemake fills in {properties} as
# a JSON comment and {exec_job} with the actual rule command before handing
# this off to slurm-submit.py for sbatch submission.
#
# Activates the one shared conda env directly here (each SLURM job is a
# separate sbatch submission that doesn't inherit conda activation state).
# See docs/slurm.md for why rules don't declare a per-rule `conda:` env,
# and why CONDA_ENV_DIR is passed in via the environment rather than
# computed from this script's own location.
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