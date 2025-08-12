#!/bin/bash
set -e

mdtbx="pixi run mdtbx"
packmol_memgen="pixi run packmol-memgen"
input_pdb=./input.pdb

# prepare ligand parameter if needed
# ligad_mol2=./ligand.mol2
# $mdtbx gen_am1bcc -s $ligad_mol2 -r LIG -m 1 -c 0

# prepare system
$mdtbx addace -s input.pdb -o input_ace.pdb
$mdtbx addnme -s input_ace.pdb -o input_ace_nme.pdb
$packmol_memgen \
    --pdb input_ace_nme.pdb \
    --lipids POPC:CHL1 \
    --ratio 4:1 \
    --salt \
    --salt_c Na+ \
    --salt_a Cl- \
    --saltcon 0.15 \
    --keepligs \
    --dims 120 120 150 \
    --parametrize \
    --gaff2 \
    --ffprot ff14SB \
    --ffwat tip3p \
    --fflip lipid21

# convert from amber format to gromacs
$mdtbx amb2gro -p leap.prmtop -x leap.inpcrd --type parmed

# prepare index.ndx
$mdtbx add_ndx -g gmx.gro

# protein centering
$mdtbx centering_gro -f gmx.gro -p gmx.top -c Protein -o gmx.gro

# genrate posres
$mdtbx gen_posres -g gmx.gro -p gmx.top -s protein -o posres

echo done
