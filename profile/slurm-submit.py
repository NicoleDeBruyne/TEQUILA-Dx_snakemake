#!/usr/bin/env python3
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slurm_utils import (
    parse_jobscript,
    get_job_properties,
    build_sbatch_command,
    submit_job,
)


def main():
    jobscript = parse_jobscript(sys.argv)
    job_properties = get_job_properties(jobscript)

    try:
        cmd = build_sbatch_command(job_properties, jobscript)
        job_id = submit_job(cmd)
    except Exception as e:
        sys.stderr.write(f"slurm-submit.py: failed to submit job: {e}\n")
        sys.exit(1)

    print(job_id)


if __name__ == "__main__":
    main()
