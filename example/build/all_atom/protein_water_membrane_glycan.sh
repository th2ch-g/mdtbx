#!/bin/bash
set -e

# グリカンがタンパク質と同一PDBに含まれている前提
input_structure="input.pdb"
out_dir="${PWD}/gmx"

mdtbx addh -s ${input_structure} -o input_h

mdtbx addace -s input_h.pdb -o ace

mdtbx addnme -s ace.pdb -o ace_nme

mdtbx find_bond -s ace_nme.pdb -o bonds.txt -op cym.pdb

# SS-bond有無で入力PDBを切り替え
if [ -s bonds.txt ]; then
    input_pdb="cym.pdb"
else
    input_pdb="ace_nme.pdb"
fi

# --keepligs でグリカン座標を保持、--gaff2 は使わずGLYCAM06-jをtleapでロード
mdtbx cmd packmol-memgen \
    --pdb ${input_pdb} \
    --lipids POPC:CHL1 \
    --ratio 4:1 \
    --salt \
    --salt_c Na+ \
    --salt_a Cl- \
    --saltcon 0.15 \
    --keepligs \
    --notprotonate \
    --dims 120 120 200 \
    --ffwat tip3p \
    --ffprot ff14SB \
    --fflip lipid21 \
    --leapline "source leaprc.GLYCAM_06j-1"

mdtbx amb2gro -p bilayer_${input_pdb%.pdb}_lipid.top -x bilayer_${input_pdb%.pdb}_lipid.crd --type parmed

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

rm -f input_h.pdb ace.pdb ace_nme.pdb cym.pdb bonds.txt \
    bilayer_*.pdb bilayer_*.pdb_FORCED bilayer_*.crd bilayer_*.top \
    *in_EMBED.pdb *in_memembed.log \
    leap_.log packmol-memgen.log packmol.inp packmol.log gmx.pdb

echo "done"
