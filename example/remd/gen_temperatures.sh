#!/bin/bash
# gen_temperatures.sh
# Generate the temperature ladder for a REMD simulation
#
# Algorithm: the method of de Pablo et al.
# (virtualchemistry.org/remd-temperature-generator)
# Temperatures reaching the target exchange probability are computed
# iteratively from the degrees of freedom, the constraints and the solvent
# model of the system
#
# Prerequisite: read the atom and water counts off gmx.gro / gmx.top
#   water molecules: grep -c "SOL" gmx.gro (3 atoms per molecule, so divide by 3)
#   protein atoms:   grep -c "Protein" gmx.ndx or similar
#
# Example:
#   bash gen_temperatures.sh

set -e

# -----------------------------------------------------------------------
# Standard settings for a protein in aqueous solution (ff14SB + TIP3P/SPC)
# --pc 1 : constrain only the hydrogen bonds in the protein (typical
#          LINCS/SHAKE setting)
# --wc 3 : water is a rigid model (TIP3P/SPC/SPC-E)
# --hff 0 : an all-hydrogen force field (ff14SB)
# -----------------------------------------------------------------------
NW=10000    # water molecules
NP=3000     # protein atoms
TLOW=300.0  # lowest temperature [K]
THIGH=400.0 # highest temperature [K]
PDES=0.25   # target exchange probability (0.2-0.3 is typical)

echo "=== Standard protein-water REMD ==="
mdtbx gen_temperatures \
    --pdes ${PDES} \
    --tlow ${TLOW} \
    --thigh ${THIGH} \
    --nw ${NW} \
    --np ${NP} \
    --pc 1 \
    --wc 3 \
    --hff 0

# -----------------------------------------------------------------------
# Small molecules and peptides: few atoms, so the exchange probability tends
# to be high. That yields few replicas, so narrow the temperature range or
# lower pdes to compensate
# -----------------------------------------------------------------------
echo ""
echo "=== Small peptide REMD ==="
mdtbx gen_temperatures \
    --pdes 0.20 \
    --tlow 280.0 \
    --thigh 380.0 \
    --nw 3000 \
    --np 300 \
    --pc 1 \
    --wc 3 \
    --hff 0

# -----------------------------------------------------------------------
# Fully flexible (unconstrained) settings
# --pc 0 : no constraints on the protein
# --wc 0 : fully flexible water as well (SPC/Fw and similar)
# The extra degrees of freedom tend to increase the replica count
# -----------------------------------------------------------------------
echo ""
echo "=== Fully flexible (no constraints) ==="
mdtbx gen_temperatures \
    --pdes 0.25 \
    --tlow 300.0 \
    --thigh 400.0 \
    --nw 10000 \
    --np 3000 \
    --pc 0 \
    --wc 0 \
    --hff 0

echo "All done."
