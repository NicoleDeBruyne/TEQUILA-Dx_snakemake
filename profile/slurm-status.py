#!/usr/bin/env python3
"""
profile/slurm-status.py
Snakemake's --cluster-status script. Called periodically as:
    slurm-status.py <jobid>

Uses `sacct` (SLURM accounting)
"""

import subprocess
import sys

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
        return "RUNNING"

    for line in result.stdout.strip().splitlines():
        parts = line.split("|")
        if len(parts) >= 2 and parts[0] == job_id:
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
        print("running")


if __name__ == "__main__":
    main()
