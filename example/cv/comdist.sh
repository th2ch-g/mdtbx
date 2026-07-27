#!/bin/bash
# comdist.sh
# Compute the centre-of-mass distance (COM distance) between two atom groups
#
# Output: .npy (shape: [n_frames], unit: nm)
#
# Example:
#   bash comdist.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# Basic case: computed with MDtraj
# Centre-of-mass distance between the protein and the ligand
# -----------------------------------------------------------------------
mdtbx comdist \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s1 "protein" \
    -s2 "resname LIG" \
    -o cvs/comdist.npy

echo "comdist done -> cvs/comdist.npy"

# -----------------------------------------------------------------------
# Gromacs interface (--gmx): faster for large systems
# Give ndx group names to -s1/-s2
# -----------------------------------------------------------------------
# mdtbx comdist \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -s1 "Protein" \
#     -s2 "LIG" \
#     --gmx \
#     -idx gmx.ndx \
#     -o cvs/comdist_gmx.npy

echo "All done."
