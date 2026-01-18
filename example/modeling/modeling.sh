#!/bin/bash
#SBATCH --ntasks-per-node 16
#SBATCH -J modeling
#SBATCH -o modeling.out
#SBATCH -e modeling.err
#SBATCH -p q1
#SBATCH --gres=gpu:1
set -e

source /home/apps/Modules/init/bash
module purge
module load localcolabfold

seq=$(tail -n +2 target.fa)

echo "seq: $seq"

mdtbx modeling_cf -i reference.pdb -s $seq

echo done

