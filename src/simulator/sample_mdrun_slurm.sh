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

GMX_CMD="gmx"
MAXWARN=10
PRODUCTION_STEPS=10
MDRUN_OPTION="-dlb no -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu"
SIMULATION_CONTINUE=true
# continue simulation => true
# start from scratch => false

echo "hostname: $(hostname)"
echo "cpus: $SLURM_CPUS_ON_NODE"
echo "gpus: $SLURM_GPUS_ON_NODE"

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

echo done
