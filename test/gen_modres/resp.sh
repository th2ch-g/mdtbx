#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -e resp.log
#SBATCH -o resp.log
#SBATCH --mem-per-cpu=80G
set -e

# Tip1: NH(mainchain), CA, CB,..(sidechain), CO(mainchain)
# Tip2: run structure optimization
# Tip3: modify bond format in CGX.prep like (S => M)

mkdir CGX
cd CGX
mdtbx gen_modres_resp \
    -s CGX.mol2 \
    -r CGX \
    -m 1 \
    -c 0 \
    --sepbond1 N1 C2 \
    --sepbond2 C7 N2
cd ..

cat <<EOF | tleap -f -
source leaprc.protein.ff14SB
source leaprc.gaff2
loadAmberPrep   CGX/CGX.prep
loadamberparams CGX/CGX_gaff2.frcmod
loadamberparams CGX/CGX_parm10.frcmod
SYS = sequence { ALA CGX }
charge SYS
savepdb SYS leap.pdb
saveamberparm SYS leap.parm7 leap.rst7
quit
EOF


echo done
