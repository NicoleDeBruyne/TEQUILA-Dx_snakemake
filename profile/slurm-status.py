#!/usr/bin/env python3
"""
profile/slurm-status.py
Snakemake's --cluster-status script. Called periodically as:
    slurm-status.py <jobid>
Must print exactly one of: running, success, failed
so Snakemake knows whether to keep waiting, mark the job done, or retry/abort.

Uses `sacct` (SLURM accounting) rather than `squeue`, since squeue only shows
currently-queued/running jobs and loses the job the moment it finishes —
sacct retains job history and reliably reports the final exit state.
"""

import subprocess
import sys

# SLURM job states that map to each Snakemake-reported status.
RUNNING_STATES = {"PENDING", "RUNNING", "SUSPENDED", "COMPLETING", "CONFIGURING"}
SUCCESS_STATES = {"COMPLETED"}
FAILED_STATES = {
    "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL",
    "OUT_OF_MEMORY", "BOOT_FAIL", "DEADLINE", "PREEMPTED",
}


def get_job_state(job_id):
    """Query sacct for the job's current state. Returns the raw SLURM state
    string (e.g. 'COMPLETED', 'RUNNING', 'FAILED')."""
    cmd = [
        "sacct",
        "-j", job_id,
        "--format=JobID,State",
        "--noheader",
        "--parsable2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 or not result.stdout.strip():
        # sacct unavailable or job not found yet — treat as still running
        # rather than failing the whole pipeline on a transient query issue.
        return "RUNNING"

    # sacct lists the main job plus sub-steps (e.g. "12345.batch"); the
    # first line (bare job ID, no suffix) reflects the overall job state.
    for line in result.stdout.strip().splitlines():
        parts = line.split("|")
        if len(parts) >= 2 and parts[0] == job_id:
            # State can include extra info like "CANCELLED by 12345"; take
            # just the leading word.
            return parts[1].split()[0]

    return "RUNNING"


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: slurm-status.py <jobid>\n")
        sys.exit(1)

    job_id = sys.argv[1]
    state = get_job_state(job_id)

    if state in SUCCESS_STATES:
        print("success")
    elif state in FAILED_STATES:
        print("failed")
    elif state in RUNNING_STATES:
        print("running")
    else:
        # Unknown state — don't prematurely fail the pipeline; let Snakemake
        # keep polling.
        print("running")


if __name__ == "__main__":
    main()
