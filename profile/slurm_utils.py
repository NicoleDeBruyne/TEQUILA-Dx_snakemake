"""
profile/slurm_utils.py
Shared helpers for slurm-submit.py and slurm-status.py.
"""

import glob
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

    Computed here rather than inside the generated jobscript: this file is
    invoked directly by Python, so __file__ reliably reflects its real
    location. The jobscript, by contrast, is submitted via sbatch, and
    SLURM copies job scripts into its own spool directory before running
    them -- so self-locating via $0/$BASH_SOURCE inside a running job
    doesn't give the script's real path, only the spool copy's."""
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

    # Build a readable job/log name: rule + sample wildcard only.
    # NOTE: wildcards can include path-like values (e.g. `outdir`, which is
    # a full filesystem path used as a wildcard in this pipeline's rules).
    # Joining ALL wildcard values blindly previously produced filenames
    # containing "/" characters, which broke --output (sbatch silently
    # created directories instead of log files). Only use short,
    # filename-safe wildcards here — primarily `sample`.
    SAFE_WILDCARD_KEYS = {"sample", "tissue"}
    safe_values = [
        str(v) for k, v in wildcards.items()
        if k in SAFE_WILDCARD_KEYS
    ]
    wildcard_str = "_".join(safe_values)
    job_label = f"{rule_name}_{wildcard_str}" if wildcard_str else rule_name

    # Route logs into a per-sample subdirectory when a sample wildcard
    # exists, so each run's directory is organized by sample:
    #   ${LOG_DIR}/<sample>/slurm-<jobid>_<rule>_<sample>.out
    # Falls back to LOG_DIR directly for jobs with no sample wildcard.
    sample = wildcards.get("sample")
    target_dir = os.path.join(log_dir, sample) if sample else log_dir
    os.makedirs(target_dir, exist_ok=True)

    # `%j` in --output below is the SLURM job ID, which is different every
    # submission -- so a rerun (retry, or a person re-launching a failed/updated
    # job) leaves the previous attempt's log file behind under a different name.
    # Remove any prior log(s) for this exact rule+sample before submitting the
    # new one, so only the most recent attempt's log file remains.
    stale_pattern = os.path.join(target_dir, f"slurm-*_{job_label}.out")
    for stale_log in glob.glob(stale_pattern):
        try:
            os.remove(stale_log)
        except OSError:
            pass

    cmd = [
        "sbatch",
        "--parsable",
        f"--job-name=smk_{job_label}",
        f"--cpus-per-task={threads}",
        f"--mem={mem_mb}M",
        f"--time={runtime_to_hms(runtime)}",
        f"--exclude={EXCLUDE_NODES}",
        f"--output={target_dir}/slurm-%j_{job_label}.out",
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