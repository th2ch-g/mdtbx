#!/bin/bash
set -e

input_structure="input.pdb"
out_dir="${PWD}/gmx"

mdtbx addh -s ${input_structure} -o input_h

mdtbx addace -s input_h.pdb -o ace

mdtbx addnme -s ace.pdb -o ace_nme

mdtbx find_bond -s ace_nme.pdb -o bonds.txt -op cym.pdb

# SS-bond有無で入力PDBとpostcmdを切り替え
if [ -s bonds.txt ]; then
    input_pdb="cym.pdb"
    mdtbx build_solution \
        -i ${input_pdb} \
        -o ./ \
        --ion_conc 0.15 \
        --cation Na+ \
        --anion Cl- \
        --boxsize 100 100 100 \
        --addpostcmd "$(cat bonds.txt)"
else
    input_pdb="ace_nme.pdb"
    mdtbx build_solution \
        -i ${input_pdb} \
        -o ./ \
        --ion_conc 0.15 \
        --cation Na+ \
        --anion Cl- \
        --boxsize 100 100 100
fi

mdtbx amb2gro -p leap.parm7 -x leap.rst7 --type parmed

mdtbx add_ndx -g gmx.gro

mdtbx centering_gro -f gmx.gro -p gmx.top -c Protein

mdtbx gen_posres -p gmx.top -s "protein and backbone" -o posres

mdtbx rmfile

mkdir ${out_dir}
mv gmx.gro ${out_dir}/
mv gmx.top ${out_dir}/
mv *itp ${out_dir}/
mv *.ndx ${out_dir}/
cp mdps/*.mdp ${out_dir}/
cp mdrun_slurm.sh ${out_dir}/

rm -f leap.parm7 leap.rst7 leap.pdb gmx.pdb
rm -f input_h.pdb ace.pdb ace_nme.pdb cym.pdb bonds.txt

echo "done"
