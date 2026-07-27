#!/bin/bash
set -e

input_structure="input.pdb"

mdtbx addh -s ${input_structure} -o input_h

mdtbx addace -s input_h.pdb -o ace

mdtbx addnme -s ace.pdb -o ace_nme

mdtbx find_bond -s ace_nme.pdb -o bonds.txt -op cym.pdb

# Switch the input PDB and the postcmd depending on whether SS-bonds exist
if [ -s bonds.txt ]; then
    input_pdb="cym.pdb"
    mdtbx build_vacuum \
        -i ${input_pdb} \
        -o ./ \
        --addpostcmd "$(cat bonds.txt)"
else
    input_pdb="ace_nme.pdb"
    mdtbx build_vacuum \
        -i ${input_pdb} \
        -o ./
fi

rm -f input_h.pdb ace.pdb ace_nme.pdb cym.pdb bonds.txt

echo "done"
