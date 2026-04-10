#!/bin/bash
# opt_perf.sh
# Optimize Gromacs mdrun settings with mdtbx opt_perf.
#
# This example assumes that run_trial.sh prepares trial-local inputs and
# launches one short benchmark run.
#
# Usage:
#   bash opt_perf.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TRIAL_HELPER="${SCRIPT_DIR}/run_trial.sh"
TRIAL_ROOT="$(pwd)/opt_perf_trials"

export INPUT_TPR="$(pwd)/prd.tpr"

mdtbx opt_perf \
    --mdrun-command "bash ${TRIAL_HELPER} {trial_dir} {log_path} {n_gpu} {n_core} {ntomp} {ntmpi} {gpu_ids}" \
    --workdir "${TRIAL_ROOT}" \
    --gpu 1 2 \
    --core 8 16 \
    --omp 2 4 8 \
    --mpi 1 2 4 \
    --sampler grid \
    --history-output opt_perf_history.csv \
    -o opt_perf_best.json

# External MPI examples:
# mdtbx opt_perf \
#     --mpi-launcher "mpirun -np {ntmpi}" \
#     --mdrun-command "gmx_mpi mdrun -deffnm prd -gpu_id {gpu_ids}" \
#     --gpu 1 2 \
#     --core 8 16 \
#     --omp 2 4 8 \
#     --mpi 1 2 4
#
# mdtbx opt_perf \
#     --mpi-launcher "srun --ntasks {ntmpi} --cpus-per-task {ntomp}" \
#     --mdrun-command "gmx_mpi mdrun -deffnm prd -gpu_id {gpu_ids}" \
#     --gpu 1 2 \
#     --core 8 16 \
#     --omp 2 4 8 \
#     --mpi 1 2 4

echo "opt_perf done -> opt_perf_best.json, opt_perf_history.csv"
