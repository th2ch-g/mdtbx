#!/bin/bash
# comvec.sh
# Compute the centre-of-mass vector (COM vector) between two atom groups
#
# comdist returns a scalar (a distance); comvec returns a 3D vector, which
# suits directional motion such as channel permeation or membrane insertion
#
# Output: .npy (shape: [n_frames, 3], unit: nm)
#
# Example:
#   bash comvec.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# Basic case: computed with MDtraj
# -----------------------------------------------------------------------
mdtbx comvec \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s1 "protein" \
    -s2 "resname LIG" \
    -o cvs/comvec.npy

echo "comvec done -> cvs/comvec.npy"

# -----------------------------------------------------------------------
# Gromacs interface (--gmx)
# -----------------------------------------------------------------------
# mdtbx comvec \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -s1 "Protein" \
#     -s2 "LIG" \
#     --gmx \
#     -idx gmx.ndx \
#     -o cvs/comvec_gmx.npy

echo "All done."
