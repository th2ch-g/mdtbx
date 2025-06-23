#!/bin/bash
#SBATCH -p all_q
#SBATCH -n 16
#SBATCH -J test_md
#SBATCH -o test_md.out
#SBATCH -e test_md.err
#SBATCH --gres=gpu:1
set -e

source $MODULESHOME/init/bash
module purge
module load gromacs/2025.1

restart_gro="gmx.gro"
MAXWARN=10
gmx_cmd="gmx"

# minimization
i=0
if [ ! -f "step6.${i}_minimization.gro" ]; then
    $gmx_cmd grompp \
        -f step6.${i}_minimization.mdp \
        -c ${restart_gro} \
        -r gmx.gro \
        -n index.ndx \
        -p gmx.top \
        -maxwarn ${MAXWARN} \
        -o step6.${i}_minimization.tpr
    $gmx_cmd run -v -deffnm step6.${i}_minimization
fi
restart_gro="step6.${i}_minimization.gro"

# nvt npt
for i in {1..6};
do
    if [ ! -f "step6.${i}_equilibration.gro" ]; then
        $gmx_cmd grompp \
            -f step6.${i}_equilibration.mdp \
            -c ${restart_gro} \
            -r gmx.gro \
            -n index.ndx \
            -p gmx.top \
            -maxwarn ${MAXWARN} \
            -o step6.${i}_equilibration.tpr
        $gmx_cmd run -v -deffnm step6.${i}_equilibration
    fi
    restart_gro="step6.${i}_equilibration.gro"
done

# production
for i in {1..10};
do
    if [ ! -f "step7_${i}.gro" ]; then
        # Use restraints option with force=0 for system converted by acpype
        #   restraints force is 0 during production MD
        #   restraints is enabled to prevent segmentation fault
        $gmx_cmd grompp \
            -f step7_production.mdp \
            -c ${restart_gro} \
            -n index.ndx \
            -p gmx.top \
            -maxwarn ${MAXWARN} \
            -o step7_${i}.tpr
            # -r gmx.gro \
        $gmx_cmd mdrun -v -deffnm step7_${i} -dlb no
    fi
    restart_gro="step7_${i}.gro"
done

echo done
