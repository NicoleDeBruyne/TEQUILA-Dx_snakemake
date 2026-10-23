#!/bin/bash
# properties = {properties}

if [ -z "$CONDA_ENV_DIR" ]; then
    echo "WARNING: CONDA_ENV_DIR is not set in the job environment; the conda env will not be activated for this job." >&2
elif [ -d "$CONDA_ENV_DIR/bin" ]; then
    export PATH="$CONDA_ENV_DIR/bin:$PATH"
else
    echo "WARNING: $CONDA_ENV_DIR/bin not found; set conda_env_dir in profile/slurm_settings.yaml or CONDA_ENV_DIR to the pipeline's conda env." >&2
fi

JOB_LOG_FILES=$(python3 - "$0" <<'PYEOF'
import json, sys
with open(sys.argv[1]) as fh:
    for line in fh:
        if line.startswith('# properties ='):
            props = json.loads(line.split('=', 1)[1])
            print('\n'.join(props.get('log', [])))
            break
PYEOF
) || JOB_LOG_FILES=""

set +e
{exec_job}
EXIT_CODE=$?
set -e

if [ -n "$SNAKEMAKE_SLURM_EXIT_CODE_DIR" ] && [ -n "$SLURM_JOB_ID" ]; then
    mkdir -p "$SNAKEMAKE_SLURM_EXIT_CODE_DIR" 2>/dev/null \
        && echo "$EXIT_CODE" > "$SNAKEMAKE_SLURM_EXIT_CODE_DIR/$SLURM_JOB_ID" 2>/dev/null || true
fi

if [ "$EXIT_CODE" -ne 0 ] && [ -n "$JOB_LOG_FILES" ]; then
    echo "Job failed (exit $EXIT_CODE); removing this rule's log file(s) so a rerun doesn't leave a stale one behind:" >&2
    echo "$JOB_LOG_FILES" | while IFS= read -r f; do
        [ -n "$f" ] && rm -fv -- "$f" >&2
    done
fi

exit "$EXIT_CODE"
