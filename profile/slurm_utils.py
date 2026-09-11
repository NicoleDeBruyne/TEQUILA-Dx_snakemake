"""
profile/slurm_utils.py
Shared helpers for slurm-submit.py and slurm-status.py.
"""

import os
import re
import subprocess
from snakemake.utils import read_job_properties

# Nodes excluded from every SLURM job on CHOP HPC, applied uniformly across the pipeline.
EXCLUDE_NODES = (
    "m-09-01,m-09-02,m-09-03,m-09-04,m-09-05,m-09-06,"
    "m-09-07,m-09-09,m-09-10,m-12-08"
)


def parse_jobscript(argv):
    return argv[-1]

def get_job_properties(jobscript):
    return read_job_properties(jobscript)

def runtime_to_hms(runtime_minutes):
    runtime_minutes = int(runtime_minutes)
    hours, minutes = divmod(runtime_minutes, 60)
    return f"{hours:02d}:{minutes:02d}:00"

def default_conda_env_dir():
    profile_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(profile_dir, "..", "envs", "conda_env")

def build_sbatch_command(job_properties, jobscript, log_dir):
    resources = job_properties.get("resources", {})
    wildcards = job_properties.get("wildcards", {})
    rule_name = job_properties.get("rule", "job")
    threads   = job_properties.get("threads", 1)

    mem_mb  = resources.get("mem_mb", 4000)
    runtime = resources.get("runtime", 720)

    SAFE_WILDCARD_KEYS = {"bed_id", "sample", "sample_type", "tissue"}
    safe_values = [
        str(v) for k, v in wildcards.items()
        if k in SAFE_WILDCARD_KEYS
    ]
    wildcard_str = "_".join(safe_values)
    job_label = f"{rule_name}_{wildcard_str}" if wildcard_str else rule_name

    os.makedirs(log_dir, exist_ok=True)

    cmd = [
        "sbatch",
        "--parsable",
        f"--job-name=smk_{job_label}",
        f"--cpus-per-task={threads}",
        f"--mem={mem_mb}M",
        f"--time={runtime_to_hms(runtime)}",
        f"--exclude={EXCLUDE_NODES}",
        f"--output={log_dir}/slurm-%j_{job_label}.out",
        f"--export=ALL,CONDA_ENV_DIR={default_conda_env_dir()}",
        jobscript,
    ]
    return cmd

def submit_job(cmd):
    """Run sbatch and return the numeric job ID"""
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = result.stdout.strip()
    match = re.search(r"(\d+)", job_id)
    if not match:
        raise RuntimeError(f"Could not parse job ID from sbatch output: {job_id!r}")
    return match.group(1)