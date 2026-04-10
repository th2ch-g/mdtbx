#!/bin/bash
# contactmap.sh
# Calculate a residue contact matrix from representative atoms.
#
# Output: .npy (shape: [n_residues, n_residues])
# Value  : contact frequency in [0, 1]
#
# Usage:
#   bash contactmap.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"
SELECTION="protein and name CA"
CUTOFF_ANGSTROM=6.0

mkdir -p cvs

mdtbx contactmap \
    -p "${TOPOLOGY}" \
    -t "${TRAJECTORY}" \
    -s "${SELECTION}" \
    --cutoff "${CUTOFF_ANGSTROM}" \
    -o cvs/contactmap.npy

echo "contactmap done -> cvs/contactmap.npy"
