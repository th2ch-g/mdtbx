#!/bin/bash
set -e

input_structure="template.pdb"
out_dir="${PWD}/gmx"

mdtbx addace -s ${input_structure} -o ace

mdtbx addnme -s ace.pdb -o ace_nme

mdtbx build_solution \
    -i ./ace_nme.pdb \
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

mkdir ${out_dir}
mv gmx.gro ${out_dir}/
mv gmx.top ${out_dir}/
mv *itp ${out_dir}/
mv *.ndx ${out_dir}/
cp mdps/*.mdp ${out_dir}/
cp mdrun_slurm.sh ${out_dir}/

rm -f leap.parm7 leap.rst7 leap.pdb gmx.pdb
rm -f ace_nme.pdb ace.pdb

echo done
