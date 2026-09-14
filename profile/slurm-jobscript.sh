#!/bin/bash
# properties = {properties}
# profile/slurm-jobscript.sh
# Standard Snakemake jobscript template

# Activates the one shared conda env directly here
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

# Remove a failed job's `.log` (which lives in a `logs/` dir that's a sibling of `output/`)
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