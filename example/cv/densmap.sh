#!/bin/bash
# densmap.sh
# Compute a 2D density map (histogram) for a set of atoms
#
# Used to visualise ligand distribution in membrane systems, or density
# projected onto a particular plane
# --axis: the plane to project onto (xy / xz / yz)
# --bins: number of bins per axis
#
# Output: .npy (object array [counts, edges0, edges1])
#   counts  : shape [bins, bins]  - frequency per cell
#   edges0/1: shape [bins+1]      - bin edges
#
# Example:
#   bash densmap.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# Projection onto the xy plane (in-plane distribution of the membrane)
# -----------------------------------------------------------------------
mdtbx densmap \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resname LIG" \
    --axis xy \
    --bins 100 \
    -o cvs/densmap_xy.npy

echo "densmap xy done -> cvs/densmap_xy.npy"

# -----------------------------------------------------------------------
# Projection onto the xz plane (a cross-section along the membrane normal)
# -----------------------------------------------------------------------
mdtbx densmap \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resname LIG" \
    --axis xz \
    --bins 100 \
    -o cvs/densmap_xz.npy

echo "densmap xz done -> cvs/densmap_xz.npy"

# -----------------------------------------------------------------------
# Using Gromacs gmx densmap instead (--gmx)
# -----------------------------------------------------------------------
# mdtbx densmap \
#     -p gmx.tpr \
#     -t ${TRAJECTORY} \
#     -s "LIG" \
#     --gmx \
#     -idx gmx.ndx \
#     -o cvs/densmap_gmx.npy

echo "All done."
