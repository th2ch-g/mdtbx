#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -e resp.log
#SBATCH -o resp.log
#SBATCH --mem-per-cpu=80G
set -e

source /home/apps/Modules/init/bash
module purge
module load gaussian16.C02

mdtbx gen_resp \
    -s NKP.mol2 \
    -r NKP \
    -m 1 \
    -c 0

echo done
