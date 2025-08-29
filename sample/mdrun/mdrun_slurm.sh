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
module load gromacs/2025.2

SIMULATION_CONTINUE=true
SIMULATION_OVERWRITE=false
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
if [[ $SIMULATION_OVERWRITE == "true" || ! -f "mini.gro" ]]; then
    if [[ $SIMULATION_OVERWRITE == "true" || $SIMULATION_CONTINUE == "false" || ! -f "mini.tpr" ]]; then
        $GMX_CMD grompp \
            -f mini.mdp \
            -c ${restart_gro} \
            -r gmx.gro \
            -n index.ndx \
            -p gmx.top \
            -maxwarn ${MAXWARN} \
            -o mini.tpr
    fi
    $GMX_CMD mdrun -v -deffnm mini
fi
restart_gro="mini.gro"

# nvt npt
for i in {1..6};
do
    if [[ $SIMULATION_OVERWRITE == "true" || ! -f "eq${i}.gro" ]]; then
        if [[ $SIMULATION_OVERWRITE == "true" || $SIMULATION_CONTINUE == "false" || ! -f "eq${i}.tpr" ]]; then
            $GMX_CMD grompp \
                -f eq${i}.mdp \
                -c ${restart_gro} \
                -r gmx.gro \
                -n index.ndx \
                -p gmx.top \
                -maxwarn ${MAXWARN} \
                -o eq${i}.tpr
        fi
        if [[ $SIMULATION_CONTINUE == "true" && $SIMULATION_OVERWRITE == "false" ]]; then
            $GMX_CMD mdrun -v -deffnm eq${i} ${MDRUN_OPTION} -cpi eq${i}.cpt
        else
            $GMX_CMD mdrun -v -deffnm eq${i} ${MDRUN_OPTION}
        fi
    fi
    restart_gro="eq${i}.gro"
done

# production
for i in $(seq 1 ${PRODUCTION_STEPS});
do
    if [[ $SIMULATION_OVERWRITE == "true" || ! -f "prd${i}.gro" ]]; then
        # Use restraints option with force=0 for system converted by acpype
        #   restraints force is 0 during production MD
        #   restraints is enabled to prevent segmentation fault
        if [[ $SIMULATION_OVERWRITE == "true" || $SIMULATION_CONTINUE == "false" || ! -f "prd${i}.tpr" ]]; then
            $GMX_CMD grompp \
                -f prd.mdp \
                -c ${restart_gro} \
                -n index.ndx \
                -p gmx.top \
                -maxwarn ${MAXWARN} \
                -o prd${i}.tpr
                # -r gmx.gro \
        fi
        if [[ $SIMULATION_CONTINUE == "true" && $SIMULATION_OVERWRITE == "false" ]]; then
            $GMX_CMD mdrun -v -deffnm prd${i} ${MDRUN_OPTION} -cpi prd${i}.cpt
        else
            $GMX_CMD mdrun -v -deffnm prd${i} ${MDRUN_OPTION}
        fi
    fi
    restart_gro="prd${i}.gro"
done

# remove temporary files
rm -f *cpt
rm -f \#*
# remove intermediate files
rm -f mini.xtc mini.trr
rm -f eq{1..6}.xtc eq{1..6}.trr

echo done
