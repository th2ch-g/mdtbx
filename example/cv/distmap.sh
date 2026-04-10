#!/bin/bash
# distmap.sh
# Calculate an average residue distance matrix from representative atoms.
#
# Output: .npy (shape: [n_residues, n_residues], unit: angstrom)
#
# Usage:
#   bash distmap.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"
SELECTION="protein and name CA"

mkdir -p cvs

mdtbx distmap \
    -p "${TOPOLOGY}" \
    -t "${TRAJECTORY}" \
    -s "${SELECTION}" \
    -o cvs/distmap.npy

echo "distmap done -> cvs/distmap.npy"
