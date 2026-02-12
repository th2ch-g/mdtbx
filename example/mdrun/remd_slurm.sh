#!/bin/bash
#SBATCH -p all_q
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J test_md
#SBATCH -o test_md.out
#SBATCH -e test_md.err
#SBATCH --gres=gpu:1
set -e

source $MODULESHOME/init/bash
module purge
module load gcc/13.3.0 cuda/12.9 cmake/3.31.6 openmpi/5.0.7
module load /home/hori/works/tools/hpc_sdk/modulefiles/nvhpc/25.7
export PATH="$TOOLS/gromacs/2025.1-mpi/gromacs-2025.1/bin:$PATH"

OMP=1
MPI=16
REPLEX=500
N_REPLICA=16
DEFFNM="reus" # or tremd
SIMULATION_CONTINUE=true
SIMULATION_OVERWRITE=false
MAXWARN=10

# for single GPU
export GMX_CUDA_GRAPH=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_GPU_DD_COMMS=1
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_GPU_PME_PP_COMMS=1
GMX_CMD="gmx_mpi"
MDRUN_OPTION="-dlb no -pin on -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu"

# for CPU
# GMX_CMD="srun -np $SLURM_NTASKS gmx_mpi"
# MDRUN_OPTION="-dlb no -nomp $SLURM_CPUS_PER_TASK"

echo "hostname: $(hostname)"
echo "jobid: $SLURM_JOB_ID"
echo "node list: $SLURM_JOB_NODELIST"
echo "nodes: $SLURM_NNODES"
echo "cpus: $SLURM_CPUS_ON_NODE"
echo "gpus: $SLURM_GPUS_ON_NODE"
echo "omp: $SLURM_CPUS_PER_TASK"
echo "mpi: $SLURM_NTASKS"

restart_gro="gmx.gro"
eval mpirun -np $MPI \
    gmx_mpi mdrun -deffnm ${DEFFNM} -ntomp $OMP \
    -multidir rep{1..${N_REPLICA}} -replex ${REPLEX} \
    ${MDRUN_OPTION}

# remove temporary files
rm -f *cpt
rm -f \#*

echo done
