#!/bin/bash
# Runs one full snakemake run inside a single SLURM allocation with one GPU.
#
# Unlike the slurm_mix profile, no jobs are submitted per sample: snakemake
# executes all rules locally on the allocated node. GPU jobs run one at a time
# (--resources gpu=1) and rule threads are capped at the allocated cores.
#
# Usage (from the snakeMSA root, so that slurm/ exists for the job log):
#   mkdir -p slurm
#   sbatch tools/slurm_run.sh CONFIGFILE [snakemake args...]
#
# Examples:
#   sbatch tools/slurm_run.sh configs/learnMSA_homfam.json
#   sbatch tools/slurm_run.sh configs/learnMSA.json --config cleanup=true
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
    echo "Usage: sbatch $0 CONFIGFILE [snakemake args...]" >&2
    exit 1
fi
configfile="$1"
shift

snakemsa_dir="${SNAKEMSA_DIR:-$HOME/src/snakeMSA}"
snakemake_env="${SNAKEMAKE_ENV:-snakeMSA}"
conda="${CONDA_EXE:-$HOME/miniforge3/bin/conda}"
cores="${SLURM_CPUS_PER_TASK:-4}"

cd "$snakemsa_dir"

echo "Job ${SLURM_JOB_ID:-local} on $(hostname), $(date)"
echo "Config: $configfile, cores: $cores, extra args: $*"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi -L || true

"$conda" run --live-stream -n "$snakemake_env" snakemake \
    --configfile "$configfile" \
    --cores "$cores" \
    --resources gpu=1 \
    --rerun-incomplete \
    --keep-going \
    "$@"

echo "Finished $(date)"
