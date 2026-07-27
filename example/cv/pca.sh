#!/bin/bash
# pca.sh
# Extract collective motion modes by principal component analysis (PCA)
#
# Used to visualise free-energy landscapes, or as a multi-dimensional CV for
# PaCS-MD. The fitting selection (-sft/-sfr) and the PCA selection
# (-sct/-scr) can be given separately
#
# Output:
#   - .npy: scores (shape: [n_frames, n_components])
#   - .npz: PCA metadata for PyMOL visualization
#   - .pdb: average structure after fitting
#
# Example:
#   bash pca.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"
REFERENCE="ref.gro"

mkdir -p cvs

# -----------------------------------------------------------------------
# Basic case: PCA with MDtraj + scikit-learn
# Fit on the backbone, then extract the principal components of the backbone
# -----------------------------------------------------------------------
mdtbx pca \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -r ${REFERENCE} \
    -sft "protein and backbone" \
    -sfr "protein and backbone" \
    -sct "protein and backbone" \
    -scr "protein and backbone" \
    -n 10 \
    -o cvs/pca.npy \
    -oz cvs/pca_backbone.npz \
    -oa cvs/pca_average.pdb

echo "pca done -> cvs/pca.npy  (shape: [n_frames, 10])"
echo "metadata done -> cvs/pca_backbone.npz"
echo "average structure done -> cvs/pca_average.pdb"

# -----------------------------------------------------------------------
# Using Gromacs gmx covar/anaeig instead (--gmx)
# Use this for large systems, or when consistency with gmx is required
# -----------------------------------------------------------------------
# mdtbx pca \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -sft "Backbone" \
#     -sct "Backbone" \
#     --gmx \
#     -idx gmx.ndx \
#     -n 10 \
#     -o cvs/pca_gmx.npy \
#     -oz cvs/pca_gmx.npz \
#     -oa cvs/pca_gmx_average.pdb

echo "All done."
