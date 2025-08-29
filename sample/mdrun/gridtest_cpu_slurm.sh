#!/bin/bash
set -e

# this scripts was used in ISSP-B

function make_file_and_submit() {
    local part=$1
    local num_nodes=$2
    local mpi=$3
    local omp=$4

mkdir -p gridtest-$part-$num_nodes-$mpi-$omp
cd gridtest-$part-$num_nodes-$mpi-$omp
cp ../test.tpr .

    echo "#!/bin/bash
#SBATCH -N ${num_node}
#SBATCH -n ${mpi}
#SBATCH -c ${omp}
#SBATCH -p ${part}
#SBATCH -o gridtest-${part}-${num_nodes}-${mpi}-${omp}.log
#SBATCH -e gridtest-${part}-${num_nodes}-${mpi}-${omp}.err
#SBATCH -J gridtest-${part}-${num_nodes}-${mpi}-${omp}
#SBATCH --time=00:10:00
set -e

module purge
module load gromacs

export OMP_NUM_THREADS=$omp
MDRUN_OPTION=\"-dlb no -ntomp \$OMP_NUM_THREADS\"

srun -np $mpi gmx_mpi mdrun -v -deffnm test \${MDRUN_OPTION}

echo done
" > test_mdrun.sh

pjsub test_mdrun.sh
cd ..

}

CORE_PER_NODE=128

for part in L1cpu L4cpu L16cpu L36cpu L72cpu L144cpu;
do
    for num_node in {1..10};
    do
        case $part in
            F1cpu | B1cpu | L1cpu )
                if [ $num_node -ge 2 ]; then
                    continue
                fi
                ;;
            F4cpu | B4cpu | L4cpu )
                if [ $num_node -eq 1 ] || [ $num_node -ge 5 ]; then
                    continue
                fi
                ;;
            F16cpu | B16cpu | L16cpu )
                if [ $num_node -le 4 ] || [ $num_node -ge 17 ]; then
                    continue
                fi
                ;;
            F36cpu | B36cpu | L36cpu )
                if [ $num_node -le 16 ] || [ $num_node -ge 37 ]; then
                    continue
                fi
                ;;
            F72cpu | B72cpu | L72cpu )
                if [ $num_node -ne 72 ]; then
                    continue
                fi
                ;;
            F144cpu | B144cpu | L144cpu )
                if [ $num_node -ne 144 ]; then
                    continue
                fi
                ;;
            i8cpu )
                if [ $num_node -ge 9 ]; then
                    continue
                fi
                ;;
        esac

        for omp in 1 2 4 8 16 32 64;
        do
            if [[ "$num_nodes" == *"x"*"x"* ]]; then
                a=$(echo $num_nodes | awk -F "x" '{print $1 * $2 * $3}')
            elif [[ "$num_nodes" == *"x"* ]]; then
                a=$(echo $num_nodes | awk -F "x" '{print $1 * $2}')
            else
                a=$num_nodes
            fi

            mpi=$(echo $a $CORE_PER_NODE $omp | awk '{print $1 * $2 / $3 }')

            make_file_and_submit $part $num_nodes $mpi $omp
        done
    done
done


echo submit done
