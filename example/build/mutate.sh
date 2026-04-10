#!/bin/bash
# mutate.sh
# Mutate a selected residue with the PyMOL mutagenesis wizard.
#
# Usage:
#   bash mutate.sh

set -e

INPUT_STRUCTURE="input.pdb"
SELECTION="chain A and resi 42"
MUTANT="VAL"
OUTPUT_STRUCTURE="mutated_A42V.pdb"

mdtbx mutate \
    -s "${INPUT_STRUCTURE}" \
    --selection "${SELECTION}" \
    --mutant "${MUTANT}" \
    -o "${OUTPUT_STRUCTURE}"

echo "mutate done -> ${OUTPUT_STRUCTURE}"
