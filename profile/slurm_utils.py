"""
profile/slurm_utils.py
Shared helpers for slurm-submit.py and slurm-status.py.
"""

import os
import re
import subprocess
from snakemake.utils import read_job_properties


# Nodes excluded from every SLURM job, applied uniformly across the pipeline.
EXCLUDE_NODES = (
    "m-09-01,m-09-02,m-09-03,m-09-04,m-09-05,m-09-06,"
    "m-09-07,m-09-09,m-09-10,m-12-08"
)


def parse_jobscript(argv):
    """The last argument passed by Snakemake to the --cluster script is the
    generated jobscript path. Everything else is informational."""
    return argv[-1]


def get_job_properties(jobscript):
    """Read the job_properties JSON Snakemake embeds as a comment at the
    top of the generated jobscript (rule name, threads, resources, etc.)."""
    return read_job_properties(jobscript)


def runtime_to_hms(runtime_minutes):
    """Convert an integer-minutes runtime resource into SLURM's HH:MM:SS."""
    runtime_minutes = int(runtime_minutes)
    hours, minutes = divmod(runtime_minutes, 60)
    return f"{hours:02d}:{minutes:02d}:00"


def default_conda_env_dir():
    """Absolute path to envs/conda_env/, a sibling of this profile/ directory.
    Computed here (not in the jobscript) since __file__ reliably reflects
    this file's real location -- see docs/slurm.md for why the jobscript
    itself can't self-locate the same way."""
    profile_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(profile_dir, "..", "envs", "conda_env")


def build_sbatch_command(job_properties, jobscript, log_dir):
    """Construct the full sbatch command line from job_properties."""
    resources = job_properties.get("resources", {})
    wildcards = job_properties.get("wildcards", {})
    rule_name = job_properties.get("rule", "job")
    threads   = job_properties.get("threads", 1)

    mem_mb  = resources.get("mem_mb", 4000)
    runtime = resources.get("runtime", 720)  # minutes

    # Build a readable job/log name: rule + every short, filename-safe
    # wildcard used anywhere in this pipeline (bed_id, sample, sample_type,
    # tissue -- none of these are full filesystem paths, unlike e.g.
    # `outdir`, which would break --output and is deliberately excluded).
    # Including bed_id/sample_type (not just sample/tissue) matters for
    # group-level rules like _10C3_amalgam_build_transcriptome, which have
    # no `sample` wildcard at all -- without this, every group's job for
    # such a rule shared the exact same job_label.
    SAFE_WILDCARD_KEYS = {"bed_id", "sample", "sample_type", "tissue"}
    safe_values = [
        str(v) for k, v in wildcards.items()
        if k in SAFE_WILDCARD_KEYS
    ]
    wildcard_str = "_".join(safe_values)
    job_label = f"{rule_name}_{wildcard_str}" if wildcard_str else rule_name

    # All SLURM stdout/stderr logs go directly into log_dir, flat -- no
    # per-sample or other subdirectories. `%j` (SLURM job ID) is globally
    # unique per submission, so slurm-<jobid>_<job_label>.out can never
    # collide with another job's file even without one, and past attempts'
    # logs are deliberately left in place (not deleted) for later
    # debugging -- see get_job_properties's docstring / docs/slurm.md.
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
    """Run sbatch and return the numeric job ID (sbatch --parsable prints
    just the job ID to stdout)."""
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = result.stdout.strip()
    match = re.search(r"(\d+)", job_id)
    if not match:
        raise RuntimeError(f"Could not parse job ID from sbatch output: {job_id!r}")
    return match.group(1)