#!/usr/bin/env python3
"""
profile/slurm-submit.py
Snakemake's --cluster script. Called once per job as:
    slurm-submit.py <jobscript>
Must print ONLY the numeric SLURM job ID to stdout — Snakemake captures this
to correlate with --cluster-status later.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slurm_utils import (
    parse_jobscript,
    get_job_properties,
    build_sbatch_command,
    submit_job,
)

# Directory where ALL per-job SLURM stdout/stderr logs are written, flat
# (no subdirectories) -- e.g. logs/slurm/slurm-<jobid>_<rule>_<wildcards>.out.
# Normally set by the orchestrator (submit_snakemake.sh) to the same
# timestamped run directory Snakemake's own run log lives under, so every
# log from one invocation lands in one place; falls back to logs/slurm/
# under the current directory if run without that orchestrator.
LOG_DIR = os.environ.get(
    "SNAKEMAKE_SLURM_LOG_DIR",
    os.path.join(os.getcwd(), "logs", "slurm"),
)


def main():
    jobscript = parse_jobscript(sys.argv)
    job_properties = get_job_properties(jobscript)
    cmd = build_sbatch_command(job_properties, jobscript, LOG_DIR)

    try:
        job_id = submit_job(cmd)
    except Exception as e:
        sys.stderr.write(f"slurm-submit.py: failed to submit job: {e}\n")
        sys.exit(1)

    # Snakemake requires ONLY the job ID on stdout.
    print(job_id)


if __name__ == "__main__":
    main()