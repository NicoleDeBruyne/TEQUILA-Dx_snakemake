#!/usr/bin/env python3
import os
import re
import subprocess
import sys

RUNNING_STATES = {
    "PENDING", "RUNNING", "SUSPENDED", "COMPLETING", "CONFIGURING",
    "REQUEUED", "REQUEUE_FED", "REQUEUE_HOLD", "RESIZING", "SIGNALING",
    "STAGE_OUT", "STOPPED", "RESV_DEL_HOLD",
}
SUCCESS_STATES = {"COMPLETED"}
FAILED_STATES = {
    "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY",
    "BOOT_FAIL", "DEADLINE", "PREEMPTED", "REVOKED", "SPECIAL_EXIT",
}

EXIT_CODE_DIR = os.environ.get(
    "SNAKEMAKE_SLURM_EXIT_CODE_DIR",
    os.path.join(os.getcwd(), ".snakemake", "slurm_exit_codes"),
)


def _run(cmd):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    return result.stdout


def state_from_sacct(job_id):
    out = _run(["sacct", "-j", job_id, "--format=JobID,State", "--noheader", "--parsable2"])
    if not out:
        return None
    for line in out.strip().splitlines():
        parts = line.split("|")
        if len(parts) >= 2 and parts[0] == job_id and parts[1].strip():
            return parts[1].split()[0].rstrip("+")
    return None


def state_from_scontrol(job_id):
    out = _run(["scontrol", "show", "job", "-o", job_id])
    if not out:
        return None
    match = re.search(r"\bJobState=(\S+)", out)
    return match.group(1) if match else None


def state_from_exit_code_file(job_id):
    path = os.path.join(EXIT_CODE_DIR, job_id)
    try:
        with open(path) as fh:
            code = fh.read().strip()
    except OSError:
        return None
    return "COMPLETED" if code == "0" else "FAILED"


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: slurm-status.py <jobid>\n")
        sys.exit(1)

    job_id = sys.argv[1]
    state = (
        state_from_sacct(job_id)
        or state_from_scontrol(job_id)
        or state_from_exit_code_file(job_id)
    )

    if state is None:
        print("running")
        return

    if state in SUCCESS_STATES:
        print("success")
    elif state in FAILED_STATES:
        print("failed")
    else:
        if state not in RUNNING_STATES:
            sys.stderr.write(f"slurm-status.py: unrecognised state {state!r} for job {job_id}; treating as running\n")
        print("running")


if __name__ == "__main__":
    main()
