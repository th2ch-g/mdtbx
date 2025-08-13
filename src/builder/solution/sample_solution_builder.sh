#!/bin/bash
set -e

mdtbx="pixi run mdtbx"
input_pdb=./input.pdb

# prepare ligand parameter if needed
# ligad_mol2=./ligand.mol2
# $mdtbx gen_am1bcc -s $ligad_mol2 -r LIG -m 1 -c 0

# prepare system
$mdtbx addace -s input.pdb -o input_ace.pdb
$mdtbx addnme -s input_ace.pdb -o input_ace_nme.pdb
$mdtbx build_solution \
    -f input_ace_nme.pdb \
    -o ./ \
    --ion_conc 0.15 \
    --cation Na+ \
    --anion Cl- \
    --boxsize 100 100 100

# convert from amber format to gromacs
$mdtbx amb2gro -p leap.prmtop -x leap.inpcrd --type parmed

# prepare index.ndx
$mdtbx add_ndx -g gmx.gro

# protein centering
$mdtbx centering_gro -f gmx.gro -p gmx.top -c Protein -o gmx.gro

# genrate posres
$mdtbx gen_posres -g gmx.gro -p gmx.top -s protein -o posres

echo done
