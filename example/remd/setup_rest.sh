#!/bin/bash
set -e

TEMPLATE_MDP="$TOOLS/mdtbx/example/mdp/solution/prd.mdp"
SELECTION_TEMPERING="protein and (resid 117 to 137)"
N_REPLICA=16
TMIN=310
TMAX=500
TOPOLOGY="gmx.top"
INDEX="index.ndx"
ITP="*.itp"
MAXWARN=10

source /home/apps/Modules/init/bash
modle purge
module load $TOOLS/plumed-2.10.0/build/lib/plumed/modulefile
module load gcc/13.3.0 cuda/12.9 cmake/3.31.6 openmpi/5.0.7
export PATH=$TOOLS/gromacs/2022.5-mpi-plumed/gromacs-2022.5/bin:$PATH

touch plumed.dat

# ref: https://www.ag.kagawa-u.ac.jp/charlesy/2020/02/17/plumed-patched-gromacs%E3%81%AE%E3%82%A4%E3%83%B3%E3%82%B9%E3%83%88%E3%83%BC%E3%83%AB/
for rep in $(seq 1 $N_REPLICA);
do
    mkdir rep${rep}
    cp $STRUCTURE rep${rep}/gmx.gro
    cp $INDEX rep${rep}/index.ndx
    cp $ITP rep${rep}/
    cp $TEMPLATE_MDP rep${rep}/prd.mdp
    touch rep${rep}/plumed.dat
    LAMBDA=$(awk -v i=$((rep-1)) -v n=$N_REPLICA -v tmin=$TMIN -v tmax=$TMAX 'BEGIN {
        k = log(tmax/tmin) / (n - 1)
        temp_i = tmin * exp(i * k)
        print tmin / temp_i
    }')
    gmx_mpi grompp -f rep${rep}/prd.mdp -c rep${rep}/gmx.gro -p $TOPOLOGY -n rep${rep}/index.ndx -pp rep${rep}/gmx_pre.top
    mdtbx partial_tempering -s "$SELECTION_TEMPERING" -p rep${rep}/gmx_pre.top -o rep${rep}/gmx_pre2.top
    plumed partial_tempering $LAMBDA < rep${rep}/gmx_pre2.top > rep${rep}/gmx.top
    gmx_mpi grompp \
        -f rep${rep}/prd.mdp \
        -c rep${rep}/gmx.gro \
        -n rep${rep}/index.ndx \
        -p rep${rep}/gmx.top \
        -maxwarn ${MAXWARN} \
        -o rep${rep}/rest.tpr
    rm -f rep${rep}/gmx_pre2.top rep${rep}/gmx_pre.top
    rm -f topol.tpr
done

echo done
