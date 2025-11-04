#!/bin/bash
set -e

# this scripts was used in flow-fx

function make_file_and_submit() {
    local part=$1
    local num_nodes=$2
    local mpi=$3
    local omp=$4
    local node_shape=$5

mkdir -p gridtest-$part-$num_nodes-$mpi-$omp-$node_shape
cd gridtest-$part-$num_nodes-$mpi-$omp-$node_shape
cp ../test.tpr .

    echo "#!/bin/bash
#PJM -L rscgrp=$part
#PJM -L elapse=00:60:00
#PJM -L node=$num_nodes:$node_shape
#PJM --mpi proc=$mpi
#PJM -j
#PJM -S
set -e

source /usr/share/Modules/init/bash
module purge
module load gcc/10.4.0
module load fjmpi-gcc/10.4.0
module load gromacs-gcc/2022.4

export OMP_NUM_THREADS=$omp
MDRUN_OPTION=\"-dlb no -ntomp \$OMP_NUM_THREADS\"

mpiexec gmx_mpi mdrun -v -deffnm test \${MDRUN_OPTION}

echo done
" > test_mdrun.sh

pjsub test_mdrun.sh
cd ..

}

CORE_PER_NODE=48

for part in fx-small;
do
    for num_nodes in 1 2 4 2x2 8 2x2x2 4x2 16 2x2x4 4x4 2x8 24 2x2x6 4x6 2x12;
    do
        for omp in 1 2 4 8 16 48;
        do
            if [[ "$num_nodes" == *"x"*"x"* ]]; then
                a=$(echo $num_nodes | awk -F "x" '{print $1 * $2 * $3}')
            elif [[ "$num_nodes" == *"x"* ]]; then
                a=$(echo $num_nodes | awk -F "x" '{print $1 * $2}')
            else
                a=$num_nodes
            fi

            mpi=$(echo $a $CORE_PER_NODE $omp | awk '{print $1 * $2 / $3 }')

            if [[ "$num_nodes" == *"x"* ]]; then
                for node_shape in torus mesh;
                do
                    make_file_and_submit $part $num_nodes $mpi $omp $node_shape
                done
            else
                node_shape=noncont
                make_file_and_submit $part $num_nodes $mpi $omp $node_shape
            fi
        done
    done
done


echo submit done
