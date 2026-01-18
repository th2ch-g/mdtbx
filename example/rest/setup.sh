#!/bin/bash
set -e

TEMPLATE_MDP="$TOOLS/mdtbx/example/mdp/solution/prd.mdp"
SUBMIT_SCRIPT="$TOOLS/mdtbx/example/mdrun/remd_plumed_slurm.sh"
N_REPLICA=16
TMIN=310
TMAX=500
TOPOLOGY="gmx.top"
INDEX="index.ndx"
ITP="*.itp"

source /home/apps/Modules/init/bash
modle purge
module load $TOOLS/plumed-2.10.0/build/lib/plumed/modulefile
module load gcc/13.3.0 cuda/12.9 cmake/3.31.6 openmpi/5.0.7
export PATH=$TOOLS/gromacs/2022.5-mpi-plumed/gromacs-2022.5/bin:$PATH

# ref: https://www.ag.kagawa-u.ac.jp/charlesy/2020/02/17/plumed-patched-gromacs%E3%81%AE%E3%82%A4%E3%83%B3%E3%82%B9%E3%83%88%E3%83%BC%E3%83%AB/
for rep in $(seq 1 $N_REPLICA);
do
    mkdir rep${rep}
    cp $STRUCTURE rep${rep}/gmx.gro
    cp $INDEX rep${rep}/index.ndx
    cp $ITP rep${rep}/
    cp $TEMPLATE_MDP rep${rep}/prd.mdp
    cp $SUBMIT_SCRIPT rep${rep}/
    touch rep${rep}/plumed.dat
    LAMBDA=$(awk -v i=$((rep-1)) -v n=$N_REPLICA -v tmin=$TMIN -v tmax=$TMAX 'BEGIN {
        k = log(tmax/tmin) / (n - 1)
        temp_i = tmin * exp(i * k)
        print tmin / temp_i
    }')
    gmx grompp -f prd.mdp -c gmx.gro -p gmx.top -pp rep${rep}/gmx_pre.top
    mdtbx partial_tempering -s "protein and resid 0 to 10" -p rep${rep}/gmx_pre.top -o rep${rep}/gmx_pre2.top
    plumed partial_tempering $LAMBDA < rep${rep}/gmx_pre2.top > rep${rep}/gmx.top
    gmx grompp \
        -f prd.mdp \
        -c gmx.gro \
        -n index.ndx \
        -p gmx.top \
        -maxwarn ${MAXWARN} \
        -o rest.tpr
done

echo done
