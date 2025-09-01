#!/bin/bash
set -e

mdtbx addace -s template.pdb -o ace

mdtbx addnme -s ace.pdb -o ace_nme

mdtbx build_solution \
    -f ./ace_nme.pdb \
    -o ./ \
    --ion_conc 0.15 \
    --cation Na+ \
    --anion Cl- \
    --ligparam ./GXL.frcmod:./GXL.lib \
    --boxsize 100 100 100

mdtbx amb2gro -p leap.parm7 -x leap.rst7 --type parmed

mdtbx add_ndx -g gmx.gro

mdtbx centering_gro -f gmx.gro -p gmx.top -c Protein

mdtbx gen_posres -p gmx.top -s "(protein and backbone) or resname GXL" -o posres

mdtbx rmfile

mkdir gmx
mv gmx.gro gmx/
mv gmx.top gmx/
mv *itp gmx/
mv *.ndx gmx/
cp mdps/*.mdp gmx/
cp mdrun_slurm.sh gmx/

rm -f leap.parm7 leap.rst7

echo done
