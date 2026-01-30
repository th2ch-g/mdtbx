#!/bin/bash
set -e

input_structure="prot_lig.pdb"
out_dir="${PWD}/gmx"
lig_frcmod="./NKP.frcmod"
lig_lib="./NKP.lib"

mdtbx addace -s ${input_structure} -o ace

mdtbx addnme -s ace.pdb -o ace_nme

mdtbx cmd packmol-memgen \
    --pdb ace_nme.pdb \
    --lipids POPC:CHL1 \
    --ratio 4:1 \
    --salt \
    --salt_c Na+ \
    --salt_a Cl- \
    --saltcon 0.15 \
    --keepligs \
    --notprotonate \
    --dims 120 120 200 \
    --parametrize \
    --gaff2 \
    --ffwat tip3p \
    --ffprot ff14SB \
    --fflip lipid21 \
    --ligand_param ${lig_frcmod}:${lig_lib}

mdtbx amb2gro -p bilayer_ace_nme_lipid.top -x bilayer_ace_nme_lipid.crd --type parmed

mdtbx add_ndx -g gmx.gro

mdtbx centering_gro -f gmx.gro -p gmx.top -c Protein

mdtbx gen_posres -p gmx.top -s "(protein and backbone) or (resname NKP) or (resname PA or resname PC or resname OL or resname CHL)" -o posres

mdtbx rmfile

mkdir ${out_dir}
mv gmx.gro ${out_dir}
mv gmx.top ${out_dir}
mv *itp ${out_dir}/
mv *.ndx ${out_dir}/
cp mdps/*.mdp ${out_dir}/
cp mdrun_slurm.sh ${out_dir}/

rm -f ace.pdb \
    ace_nme.pdb \
    ace_nmein_EMBED.pdb \
    ace_nmein_memembed.log \
    bilayer_ace_nme.pdb \
    bilayer_ace_nme.pdb_FORCED \
    bilayer_ace_nme_lipid.crd \
    bilayer_ace_nme_lipid.pdb \
    bilayer_ace_nme_lipid.top \
    leap_.log \
    packmol-memgen.log \
    packmol.inp \
    packmol.log \
    gmx.pdb

echo done
