#!/bin/bash
# rmsf.sh
# Compute per-residue or per-atom fluctuations (RMSF)
#
# The counterpart of the B-factor; used to identify flexible and rigid regions
# --resolution residue: aggregate per residue (default)
# --resolution atom   : report per atom
#
# Output: .npy (shape: [n_residues] or [n_atoms], unit: nm)
#
# Example:
#   bash rmsf.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# Per-residue RMSF (the B-factor counterpart)
# -----------------------------------------------------------------------
mdtbx rmsf \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    --selection "protein" \
    --resolution residue \
    -o cvs/rmsf_residue.npy

echo "rmsf residue done -> cvs/rmsf_residue.npy"

# -----------------------------------------------------------------------
# Per-atom RMSF (backbone only)
# -----------------------------------------------------------------------
mdtbx rmsf \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    --selection "protein and backbone" \
    --resolution atom \
    -o cvs/rmsf_atom.npy

echo "rmsf atom done -> cvs/rmsf_atom.npy"

# -----------------------------------------------------------------------
# Using Gromacs gmx rmsf instead (--gmx)
# -----------------------------------------------------------------------
# mdtbx rmsf \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     --selection "Backbone" \
#     --resolution residue \
#     --gmx \
#     -o cvs/rmsf_gmx.npy

echo "All done."
