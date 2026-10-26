#!/bin/bash
# Runs one full snakemake run inside a single SLURM allocation with one GPU.
#
# Unlike the slurm_mix profile, no jobs are submitted per sample: snakemake
# executes all rules locally on the allocated node. GPU jobs run one at a time
# (--resources gpu=1) and rule threads are capped at the allocated cores.
#
# Usage (from the snakeMSA root, so that slurm/ exists for the job log):
#   mkdir -p slurm
#   sbatch tools/slurm_run.sh CONFIGFILE [NUM_JOBS] [dry] [snakemake args...]
#
# "dry" adds -n to snakemake (print the plan, run nothing) and disables job
# chaining. Run it directly on the login node before submitting:
#   bash tools/slurm_run.sh CONFIGFILE dry
#
# NUM_JOBS (default 1) is the total number of chained jobs. If a run does not
# finish within the time limit, the next job continues it (snakemake skips
# finished rules). Each job submits its successor at start with
# --dependency=afternotok, so the successor only runs after a timeout, node
# failure or preemption. Once snakemake returns (success or error), the
# successor is cancelled, so a failing rule does not consume the remaining jobs.
#
# Examples:
#   sbatch tools/slurm_run.sh configs/learnMSA_homfam.json
#   sbatch tools/slurm_run.sh configs/learnMSA_homfam.json 3
#   sbatch tools/slurm_run.sh configs/learnMSA.json --config cleanup=true
#   sbatch tools/slurm_run.sh configs/learnMSA.json 4 --config cleanup=true
#   sbatch --time=24:00:00 --mem=128G --partition=vision \
#       tools/slurm_run.sh configs/learnMSA.json
#
# The #SBATCH defaults below can be overridden on the sbatch command line.
# Environment variables (pass e.g. with sbatch --export=ALL,VAR=value):
#   SNAKEMSA_DIR    snakeMSA checkout (default: $HOME/src/snakeMSA)
#   SNAKEMAKE_ENV   conda env with snakemake (default: snakeMSA)
#   CONDA_EXE       conda binary (default: $HOME/miniforge3/bin/conda)
#
#SBATCH --job-name=snakeMSA
#SBATCH --partition=storm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/%x-%j.out
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: sbatch $0 CONFIGFILE [NUM_JOBS] [dry] [snakemake args...]" >&2
    exit 1
fi
configfile="$1"
shift
num_jobs=1
dry=false
while [ $# -ge 1 ]; do
    if [[ "$1" =~ ^[0-9]+$ ]]; then
        num_jobs=$((10#$1))
    elif [ "$1" = dry ]; then
        dry=true
    else
        break
    fi
    shift
done
if $dry; then
    # Dry run: snakemake -n only prints the plan. No lock is touched and no
    # follow-up jobs are submitted.
    num_jobs=1
    set -- -n "$@"
fi
if [ "$num_jobs" -lt 1 ]; then
    echo "NUM_JOBS must be at least 1" >&2
    exit 1
fi

snakemsa_dir="${SNAKEMSA_DIR:-$HOME/src/snakeMSA}"
snakemake_env="${SNAKEMAKE_ENV:-snakeMSA}"
conda="${CONDA_EXE:-$HOME/miniforge3/bin/conda}"
cores="${SLURM_CPUS_PER_TASK:-4}"

cd "$snakemsa_dir"

echo "Job ${SLURM_JOB_ID:-local} on $(hostname), $(date)"
echo "Config: $configfile, cores: $cores, jobs left: $num_jobs, extra args: $*"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi -L || true

# Submit the successor first: at the time limit the script is killed and can
# no longer do it. Resource overrides given to sbatch are passed on.
next_job=""
if [ "$num_jobs" -gt 1 ] && [ -n "${SLURM_JOB_ID:-}" ]; then
    sbatch_opts=(
        --parsable
        --dependency="afternotok:$SLURM_JOB_ID"
        --kill-on-invalid-dep=yes
        --chdir "${SLURM_SUBMIT_DIR:-$PWD}"
        --job-name "${SLURM_JOB_NAME:-snakeMSA}"
        --cpus-per-task "$cores"
    )
    time_limit="$(squeue -h -j "$SLURM_JOB_ID" -o %l || true)"
    if [ -n "${SLURM_JOB_PARTITION:-}" ]; then
        sbatch_opts+=(--partition "$SLURM_JOB_PARTITION")
    fi
    if [ -n "$time_limit" ]; then
        sbatch_opts+=(--time "$time_limit")
    fi
    if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
        sbatch_opts+=(--mem "${SLURM_MEM_PER_NODE}M")
    elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
        sbatch_opts+=(--mem-per-cpu "${SLURM_MEM_PER_CPU}M")
    fi
    # $0 is sbatch's spool copy, so submit the script from the checkout.
    next_job="$(sbatch "${sbatch_opts[@]}" "$snakemsa_dir/tools/slurm_run.sh" \
        "$configfile" "$((num_jobs - 1))" "$@")"
    next_job="${next_job%%;*}"
    echo "Submitted follow-up job $next_job ($((num_jobs - 1)) left)"
fi

# A hard kill at the time limit can leave the snakemake lock behind.
if ! $dry; then
    "$conda" run --live-stream -n "$snakemake_env" snakemake \
        --configfile "$configfile" --cores 1 --unlock "$@" || true
fi

rc=0
"$conda" run --live-stream -n "$snakemake_env" snakemake \
    --configfile "$configfile" \
    --cores "$cores" \
    --resources gpu=1 \
    --rerun-incomplete \
    --keep-going \
    "$@" || rc=$?

# Reaching this point means no timeout: either everything is done or a rule
# failed. Either way the follow-up job would not help.
if [ -n "$next_job" ]; then
    scancel "$next_job" || true
    echo "Cancelled follow-up job $next_job"
fi

echo "Finished $(date), snakemake exit code $rc"
exit "$rc"
