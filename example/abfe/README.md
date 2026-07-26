# ABFE workflow

`mdtbx` implements a Boresch-restraint double-decoupling cycle with five BAR
legs:

- ligand charge removal in solvent,
- ligand Van der Waals removal in solvent,
- ligand charge removal in the restrained complex,
- ligand Van der Waals removal in the restrained complex,
- Boresch restraint release in the coupled complex.

The analytical standard-state restraint term is added during analysis.

## Inputs

Prepare and equilibrate two systems before setup:

- the receptor-ligand complex in solvent,
- the ligand alone in solvent.

Both topologies must use the same ligand molecule type. The production MDP
controls the integration and output protocol. mdtbx replaces only the
alchemical and Boresch settings and uses LJ-PME to avoid a separate
cutoff-dispersion long-range correction. The default
`lj-pme-comb-rule = Lorentz-Berthelot` is suitable for AMBER/GAFF; use
`--lj-pme-comb-rule Geometric` for a matching force field.

## Setup

Explicit one-based Boresch anchors are recommended:

```bash
pixi run mdtbx setup_abfe \
  -f run.mdp \
  --complex-topology complex.top \
  --complex-structure complex_equilibrated.gro \
  --solvent-topology ligand_solvent.top \
  --solvent-structure ligand_solvent_equilibrated.gro \
  --moltype LIG \
  --anchors 125 128 131 4201 4204 4207 \
  --charge-windows 12 \
  --vdw-windows 16 \
  --restraint-windows 12 \
  -o abfe
```

Alternatively, omit `--anchors` and provide an MDTraj ligand selection:

```bash
pixi run mdtbx setup_abfe \
  -f run.mdp \
  --complex-topology complex.top \
  --complex-structure complex_equilibrated.gro \
  --solvent-topology ligand_solvent.top \
  --solvent-structure ligand_solvent_equilibrated.gro \
  --moltype LIG \
  --ligand-selection 'resname LIG' \
  -o abfe
```

The automatic selector searches the receptor selection `protein` by default.
Inspect the selected anchors in `abfe/abfe_manifest.json` before production.

## Run and analyze

```bash
pixi run mdtbx run_abfe --path abfe --ntomp 1
pixi run mdtbx analyze_abfe --path abfe -b 2000
```

`analyze_abfe` combines the signed five-leg BAR estimates, the analytical
Boresch standard-state term, and the ligand symmetry number:

```bash
pixi run mdtbx analyze_abfe \
  --path abfe \
  -b 2000 \
  --symmetry-number 2
```

The result is written to `abfe/abfe_result.json` and
`abfe/abfe_result.txt`, including the binding free energy, uncertainty, and
derived dissociation constant.

## Charged ligands and extra corrections

For a charged ligand, calculate the electrostatic finite-size correction
separately, for example with a Rocklin/APBS workflow. Add each already-signed
contribution to the binding free energy with its uncertainty:

```bash
pixi run mdtbx analyze_abfe \
  --path abfe \
  -b 2000 \
  --correction finite_size_charge -1.25 0.20
```

`--correction NAME VALUE_KJ_MOL ERROR_KJ_MOL` can be repeated. mdtbx does not
run APBS automatically.
