#!/bin/bash
set -e

TEMPLATE_MDP="$TOOLS/mdtbx/example/mdp/solution/us_dist.mdp"
SUBMIT_SCRIPT="$TOOLS/mdtbx/example/mdrun/remd_slurm.sh"
TOPOLOGY="gmx.top"
INDEX="index.ndx"
ITP="*.itp"

source /home/apps/Modules/init/bash
module purge
module load gcc/13.3.0 cuda/12.9 cmake/3.31.6 openmpi/5.0.7
module load /home/hori/works/tools/hpc_sdk/modulefiles/nvhpc/25.7
export PATH="$TOOLS/gromacs/2025.1-mpi/gromacs-2025.1/bin:$PATH"

cat <<EOF > target_structures_distances.txt
init1.gro     15
init2.gro     18
init3.gro     21
init4.gro     24
init5.gro     27
init6.gro     30
init7.gro     33
init8.gro     36
init9.gro     39
init10.gro    42
init11.gro    45
init12.gro    48
init13.gro    51
init14.gro    54
init15.gro    57
init16.gro    60
EOF

N_REPLICA=$(wc -l target_structures_distances.txt | awk '{print $1}')
echo "Number of replicas: ${N_REPLICA}"

for rep in $(seq 1 $N_REPLICA);
do
    STRUCTURE=$(head -n $rep target_structures_distances.txt | tail -n 1 | awk '{print $1}')
    TARGET_DISTANCE=$(head -n $rep target_structures_distances.txt | tail -n 1 | awk '{print $2}')
    mkdir replica${rep}
    cp $STRUCTURE replica${rep}/gmx.gro
    cp $TOPOLOGY replica${rep}/gmx.top
    cp $INDEX replica${rep}/index.ndx
    cp $ITP replica${rep}/
    cp $TEMPLATE_MDP replica${rep}/reus.mdp
    sed -i -e "s/TARGET_DISTANCE/${TARGET_DISTANCE}/g" replica${rep}/reus.mdp
    cp $SUBMIT_SCRIPT replica${rep}/
    mdtbx cmd gmx grompp \
        -f reus.mdp \
        -c gmx.gro \
        -n index.ndx \
        -p gmx.top \
        -maxwarn ${MAXWARN} \
        -o reus.tpr
done

rm -f target_structures_distances.txt

echo done
