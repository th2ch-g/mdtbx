#!/bin/bash
# xyz.sh
# Extract the XYZ coordinate time series of a set of atoms
#
# Use when the spatial information a scalar CV discards is needed, or as
# preprocessing for a custom CV
#
# Output: .npy (shape: [n_frames, n_atoms, 3], unit: nm)
#
# Example:
#   bash xyz.sh

set -e

TOPOLOGY="gmx.gro"
TRAJECTORY="prd.xtc"

mkdir -p cvs

# -----------------------------------------------------------------------
# All ligand atom coordinates
# -----------------------------------------------------------------------
mdtbx xyz \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resname LIG" \
    -o cvs/lig_xyz.npy

echo "xyz done -> cvs/lig_xyz.npy"

# -----------------------------------------------------------------------
# Active-site CA coordinates (selected residues)
# -----------------------------------------------------------------------
mdtbx xyz \
    -p ${TOPOLOGY} \
    -t ${TRAJECTORY} \
    -s "resid 50 to 70 and name CA" \
    -o cvs/active_site_ca_xyz.npy

echo "xyz done -> cvs/active_site_ca_xyz.npy"

echo "All done."
