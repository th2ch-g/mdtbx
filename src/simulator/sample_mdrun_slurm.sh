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
module load gromacs/2025.1

# continue simulation => true
# start from scratch => false
SIMULATION_CONTINUE=true

MAXWARN=10
PRODUCTION_STEPS=10

# for GPU
export GMX_CUDA_GRAPH=1
# export GMX_ENABLE_DIRECT_GPU_COMM=1
# export GMX_FORCE_GPU_AWARE_MPI=1
GMX_CMD="gmx"
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

# minimization
i=0
if [ ! -f "step${i}_minimization.gro" ]; then
    if [[ $SIMULATION_CONTINUE == "false" || ! -f "step${i}_minimization.tpr" ]]; then
        $GMX_CMD grompp \
            -f step${i}_minimization.mdp \
            -c ${restart_gro} \
            -r gmx.gro \
            -n index.ndx \
            -p gmx.top \
            -maxwarn ${MAXWARN} \
            -o step${i}_minimization.tpr
    fi
    $GMX_CMD mdrun -v -deffnm step${i}_minimization
fi
restart_gro="step${i}_minimization.gro"

# nvt npt
for i in {1..6};
do
    if [ ! -f "step${i}_equilibration.gro" ]; then
        if [[ $SIMULATION_CONTINUE == "false" || ! -f "step${i}_equilibration.tpr" ]]; then
            $GMX_CMD grompp \
                -f step${i}_equilibration.mdp \
                -c ${restart_gro} \
                -r gmx.gro \
                -n index.ndx \
                -p gmx.top \
                -maxwarn ${MAXWARN} \
                -o step${i}_equilibration.tpr
        fi
        if $SIMULATION_CONTINUE; then
            $GMX_CMD mdrun -v -deffnm step${i}_equilibration ${MDRUN_OPTION} -cpi step${i}_equilibration.cpt
        else
            $GMX_CMD mdrun -v -deffnm step${i}_equilibration ${MDRUN_OPTION}
        fi
    fi
    restart_gro="step${i}_equilibration.gro"
done

# production
for i in $(seq 1 ${PRODUCTION_STEPS});
do
    if [ ! -f "step7_${i}.gro" ]; then
        # Use restraints option with force=0 for system converted by acpype
        #   restraints force is 0 during production MD
        #   restraints is enabled to prevent segmentation fault
        if [[ $SIMULATION_CONTINUE == "false" || ! -f "step7_${i}.tpr" ]]; then
            $GMX_CMD grompp \
                -f step7_production.mdp \
                -c ${restart_gro} \
                -n index.ndx \
                -p gmx.top \
                -maxwarn ${MAXWARN} \
                -o step7_${i}.tpr
                # -r gmx.gro \
        fi
        if $SIMULATION_CONTINUE; then
            $GMX_CMD mdrun -v -deffnm step7_${i} ${MDRUN_OPTION} -cpi step7_${i}.cpt
        else
            $GMX_CMD mdrun -v -deffnm step7_${i} ${MDRUN_OPTION}
        fi
    fi
    restart_gro="step7_${i}.gro"
done

# remove temporary files
rm -f *cpt
rm -f \#*
# remove intermediate files
rm -f step0_minimization.xtc step0_minimization.trr
rm -f step{1..6}_equilibration.xtc step{1..6}_equilibration.trr

echo done
