import os
import re
import shlex
import subprocess

import yaml
from snakemake.utils import read_job_properties

PROFILE_DIR = os.path.dirname(os.path.abspath(__file__))
SETTINGS_FILE = os.environ.get(
    "SNAKEMAKE_SLURM_SETTINGS",
    os.path.join(PROFILE_DIR, "slurm_settings.yaml"),
)

DEFAULTS = {
    "partition": "",
    "account": "",
    "qos": "",
    "exclude_nodes": "",
    "extra_sbatch_args": [],
    "default_mem_mb": 4000,
    "default_runtime": 720,
    "log_dir": "",
    "conda_env_dir": "",
}

ENV_OVERRIDES = {
    "partition": "SNAKEMAKE_SLURM_PARTITION",
    "account": "SNAKEMAKE_SLURM_ACCOUNT",
    "qos": "SNAKEMAKE_SLURM_QOS",
    "exclude_nodes": "SNAKEMAKE_SLURM_EXCLUDE_NODES",
    "extra_sbatch_args": "SNAKEMAKE_SLURM_EXTRA_ARGS",
    "default_mem_mb": "SNAKEMAKE_SLURM_DEFAULT_MEM_MB",
    "default_runtime": "SNAKEMAKE_SLURM_DEFAULT_RUNTIME",
    "log_dir": "SNAKEMAKE_SLURM_LOG_DIR",
    "conda_env_dir": "CONDA_ENV_DIR",
}

SAFE_WILDCARD_KEYS = ("bed_id", "sample", "sample_type", "tissue")


def load_settings():
    settings = dict(DEFAULTS)
    if os.path.isfile(SETTINGS_FILE):
        with open(SETTINGS_FILE) as fh:
            loaded = yaml.safe_load(fh) or {}
        unknown = sorted(set(loaded) - set(DEFAULTS))
        if unknown:
            raise ValueError(f"Unknown key(s) in {SETTINGS_FILE}: {', '.join(unknown)}")
        settings.update({k: v for k, v in loaded.items() if v is not None})
    for key, env_var in ENV_OVERRIDES.items():
        value = os.environ.get(env_var)
        if value:
            settings[key] = value
    if isinstance(settings["extra_sbatch_args"], str):
        settings["extra_sbatch_args"] = shlex.split(settings["extra_sbatch_args"])
    if isinstance(settings["exclude_nodes"], (list, tuple)):
        settings["exclude_nodes"] = ",".join(str(n) for n in settings["exclude_nodes"])
    if not settings["log_dir"]:
        settings["log_dir"] = os.path.join(os.getcwd(), "logs", "slurm")
    if not settings["conda_env_dir"]:
        settings["conda_env_dir"] = os.path.normpath(
            os.path.join(PROFILE_DIR, "..", "envs", "conda_env")
        )
    return settings


def exit_code_dir():
    return os.environ.get(
        "SNAKEMAKE_SLURM_EXIT_CODE_DIR",
        os.path.join(os.getcwd(), ".snakemake", "slurm_exit_codes"),
    )


def parse_jobscript(argv):
    return argv[-1]


def get_job_properties(jobscript):
    return read_job_properties(jobscript)


def runtime_to_hms(runtime_minutes):
    runtime_minutes = int(float(runtime_minutes))
    hours, minutes = divmod(runtime_minutes, 60)
    return f"{hours:02d}:{minutes:02d}:00"


def _sanitize(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(text)).strip("-")


def job_label(job_properties):
    if job_properties.get("type") == "group":
        name = job_properties.get("groupid", "group")
    else:
        name = job_properties.get("rule", "job")
    wildcards = job_properties.get("wildcards", {}) or {}
    parts = [str(wildcards[k]) for k in SAFE_WILDCARD_KEYS if k in wildcards]
    label = "_".join([name] + parts)
    return _sanitize(label)[:200] or "job"


def build_sbatch_command(job_properties, jobscript, settings=None):
    settings = settings or load_settings()
    resources = job_properties.get("resources", {}) or {}
    threads = job_properties.get("threads", 1)

    mem_mb = resources.get("mem_mb", settings["default_mem_mb"])
    runtime = resources.get("runtime", settings["default_runtime"])
    partition = resources.get("slurm_partition") or settings["partition"]
    account = resources.get("slurm_account") or settings["account"]
    qos = resources.get("slurm_qos") or settings["qos"]

    label = job_label(job_properties)
    log_dir = settings["log_dir"]
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(exit_code_dir(), exist_ok=True)

    cmd = [
        "sbatch",
        "--parsable",
        f"--job-name=smk_{label}",
        f"--cpus-per-task={threads}",
        f"--mem={int(float(mem_mb))}M",
        f"--time={runtime_to_hms(runtime)}",
        f"--output={log_dir}/slurm-%j_{label}.out",
        f"--export=ALL,CONDA_ENV_DIR={settings['conda_env_dir']},"
        f"SNAKEMAKE_SLURM_EXIT_CODE_DIR={exit_code_dir()}",
    ]
    if partition:
        cmd.append(f"--partition={partition}")
    if account:
        cmd.append(f"--account={account}")
    if qos:
        cmd.append(f"--qos={qos}")
    if settings["exclude_nodes"]:
        cmd.append(f"--exclude={settings['exclude_nodes']}")
    cmd.extend(str(a) for a in settings["extra_sbatch_args"])
    cmd.append(jobscript)
    return cmd


def submit_job(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"sbatch exited with code {result.returncode}: {result.stderr.strip()}"
        )
    output = result.stdout.strip()
    match = re.match(r"(\d+)", output)
    if not match:
        raise RuntimeError(f"Could not parse job ID from sbatch output: {output!r}")
    return match.group(1)
