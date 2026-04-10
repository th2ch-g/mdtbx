#!/bin/bash
# run_trial.sh
# Helper script for opt_perf.sh.
#
# Arguments:
#   1. trial directory
#   2. log path
#   3. number of GPUs
#   4. number of CPU cores
#   5. number of OpenMP threads
#   6. number of MPI ranks
#   7. GPU id string
#
# Required files in the parent directory:
#   ${INPUT_TPR}
#
# Adjust the environment variables and mdrun options to match your cluster.

set -e

TRIAL_DIR="$1"
LOG_PATH="$2"
N_GPU="$3"
N_CORE="$4"
NTOMP="$5"
NTMPI="$6"
GPU_IDS="$7"

export OMP_NUM_THREADS="${NTOMP}"

# Example GPU-related environment variables.
export GMX_CUDA_GRAPH=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_GPU_DD_COMMS=1
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_GPU_PME_PP_COMMS=1

# Prepare trial-local inputs.
if [ -z "${INPUT_TPR}" ]; then
    echo "INPUT_TPR is not set" >&2
    exit 1
fi

ln -sf "${INPUT_TPR}" "${TRIAL_DIR}/prd.tpr"

cd "${TRIAL_DIR}"

# Edit this command for your environment.
# Examples:
#   gmx mdrun ...
#   mpirun -np "${NTMPI}" gmx_mpi mdrun ...
gmx mdrun \
    -deffnm prd \
    -s prd.tpr \
    -g "${LOG_PATH}" \
    -ntomp "${NTOMP}" \
    -ntmpi "${NTMPI}" \
    -gpu_id "${GPU_IDS}" \
    -pin on \
    -dlb no \
    -nb gpu \
    -pme gpu \
    -pmefft gpu \
    -bonded gpu \
    -update gpu

echo "trial done: n_gpu=${N_GPU} n_core=${N_CORE} ntomp=${NTOMP} ntmpi=${NTMPI}"
