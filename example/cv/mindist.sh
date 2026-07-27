#!/bin/bash
# mindist.sh
# Compute the minimum distance over all pairs between two atom groups
#
# Unlike comdist (a centre-of-mass distance), this reflects whether the two
# groups are in contact; used for contact distances at a binding site or for
# protein-protein interfaces
#
# Output: .npy (shape: [n_frames], unit: nm)
#
# Example:
#   bash mindist.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# Basic case: minimum distance between active-site residues and the ligand
# -----------------------------------------------------------------------
mdtbx mindist \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s1 "resid 50 to 70" \
    -s2 "resname LIG" \
    -o cvs/mindist.npy

echo "mindist done -> cvs/mindist.npy"

# -----------------------------------------------------------------------
# Minimum distance across a protein-protein interface
# -----------------------------------------------------------------------
# mdtbx mindist \
#     -p ${TOPOLOGY} \
#     -t ${TRAJECTORY} \
#     -s1 "chainid 0" \
#     -s2 "chainid 1" \
#     -o cvs/mindist_pp.npy

echo "All done."
