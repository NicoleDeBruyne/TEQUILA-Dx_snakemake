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

# Snakemake's own automatic cleanup-on-failure only touches a rule's
# `output:` files, never its `log:` file(s) -- so a failed job otherwise
# leaves behind a `.log` file that looks like a completed run on the next
# glance, right next to output/ (see rules/*.smk's `log:` paths, all
# living in a `logs/` dir that's a sibling of `output/`). Extract this
# job's `log:` path(s) from the properties JSON embedded above so they can
# be removed below if the job fails. This is unrelated to (and doesn't
# touch) the SLURM stdout/stderr .out file that slurm-submit.py routes to
# $SNAKEMAKE_SLURM_LOG_DIR -- those are deliberately never removed.
JOB_LOG_FILES=$(python3 -c "
import json
with open('$0') as fh:
    for line in fh:
        if line.startswith('# properties ='):
            props = json.loads(line.split('=', 1)[1])
            print('\n'.join(props.get('log', [])))
            break
")

set +e
{exec_job}
EXIT_CODE=$?
set -e

if [ "$EXIT_CODE" -ne 0 ] && [ -n "$JOB_LOG_FILES" ]; then
    echo "Job failed (exit $EXIT_CODE) -- removing this rule's log file(s) so a rerun doesn't leave a stale one behind:" >&2
    echo "$JOB_LOG_FILES" | while IFS= read -r f; do
        [ -n "$f" ] && rm -fv -- "$f" >&2
    done
fi

exit "$EXIT_CODE"