#!/bin/bash
# rmsd.sh
# Compute the RMSD against a reference structure
#
# The fitting selection (-sft/-sfr) and the RMSD selection (-sct/-scr) can be
# given separately, for example fit on the backbone and report the RMSD of an
# active-site loop
#
# Output: .npy (shape: [n_frames], unit: nm)
#
# Example:
#   bash rmsd.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"
REFERENCE="ref.gro"

mkdir -p cvs

# -----------------------------------------------------------------------
# RMSD of the whole backbone
# -----------------------------------------------------------------------
mdtbx rmsd \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -r ${REFERENCE} \
    -sft "protein and backbone" \
    -sfr "protein and backbone" \
    -sct "protein and backbone" \
    -scr "protein and backbone" \
    -o cvs/rmsd_backbone.npy

echo "rmsd backbone done -> cvs/rmsd_backbone.npy"

# -----------------------------------------------------------------------
# RMSD of an active-site loop (fit on the backbone, measure the loop)
# -sft/-sfr: selection used for fitting (a global fit)
# -sct/-scr: selection the RMSD is computed on (detects local changes)
# -----------------------------------------------------------------------
mdtbx rmsd \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -r ${REFERENCE} \
    -sft "protein and backbone" \
    -sfr "protein and backbone" \
    -sct "resid 100 to 120 and backbone" \
    -scr "resid 100 to 120 and backbone" \
    -o cvs/rmsd_loop.npy

echo "rmsd loop done -> cvs/rmsd_loop.npy"

echo "All done."
