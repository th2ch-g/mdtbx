# Membrane MDP templates

These files provide the minimization, six-stage equilibration, and production
defaults copied by `example/workflows/membrane_setup.sh`.

Packmol-Memgen systems can require gentler restraints and more equilibration
than systems prepared by other membrane builders. In particular, inspect the
packing, box dimensions, pressure-coupling behavior, and all restraint include
files before production. Do not assume that settings from a CHARMM-GUI system
transfer unchanged.

Some GROMACS/topology combinations have failed during `grompp` when multiple
position-restraint blocks were applied to an ACPYPE-converted topology. Prefer
the ParmEd conversion used by the canonical membrane workflow, keep separate
`moleculetype` sections where appropriate, and validate every equilibration
stage before removing restraints.

The equilibration templates define separate `POSRES`/`POSRES_FC` and
`POSRES_LIPID`/`POSRES_LIPID_FC` macros. The canonical membrane setup creates
both include sets. If the lipid selection produces no include file, correct
`LIPID_POSRES_SELECTION` before running `grompp`.
