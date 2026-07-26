# FEP workflow

`mdtbx` integrates a standard GROMACS alchemical workflow as three small
commands:

1. `setup_fep` creates lambda directories, MDP files, TPR files, and
   `fep_manifest.json`.
2. `run_fep` runs all or a selected range of lambda windows.
3. `analyze_fep` combines the generated Delta H files with `gmx bar`.

The input structure must already be equilibrated. The base MDP controls the
simulation protocol; `setup_fep` only replaces the alchemical settings.

## Molecule decoupling

The default `decouple` mode first removes electrostatics and then Van der Waals
interactions. It creates 23 windows from 12 electrostatic and 12 Van der Waals
states, sharing the intermediate state.

```bash
pixi run mdtbx setup_fep \
  -f example/mdp/solution/prd.mdp \
  -p gmx.top \
  -c equilibrated.gro \
  --mode decouple \
  --moltype LIG \
  --coul-windows 12 \
  --vdw-windows 12 \
  -o fep

pixi run mdtbx run_fep --path fep
pixi run mdtbx analyze_fep --path fep -b 2000
```

Use `--mode charge` or `--mode vdw` to prepare only one component. Custom
schedules are accepted with `--coul-lambdas` and `--vdw-lambdas`.

## Dual-state transformation

For a topology that already contains A-state and B-state atom parameters, use
`transform` mode:

```bash
pixi run mdtbx setup_fep \
  -f example/mdp/solution/prd.mdp \
  -p hybrid.top \
  -c equilibrated.gro \
  --mode transform \
  --fep-lambdas 0.0 0.05 0.1 0.2 0.35 0.5 0.65 0.8 0.9 0.95 1.0 \
  -o fep
```

`transform` mode does not modify or generate the dual-state topology.

## Hamiltonian replica exchange

Run all windows in one GROMACS invocation and enable exchange attempts:

```bash
pixi run mdtbx run_fep --path fep --multidir --replex 1000
```

This uses standard GROMACS Hamiltonian replica exchange. It does not implement
REST scaling.

## PLUMED FEP/REST

FEP/REST combines the three-stage alchemical path

1. remove charges from disappearing atoms,
2. transform Van der Waals parameters,
3. add charges to appearing atoms,

with REST2 scaling of a selected hot region. `setup_fep_rest` uses PLUMED
2.10.0 `partial_tempering`; `run_fep` uses the PLUMED-patched GROMACS
`-hrex` implementation.

Install both tools with `install_scripts/plumed.sh` and
`install_scripts/gmx_mpi_plumed.sh`. Load the generated PLUMED modulefile and
source the installed GROMACS `GMXRC` before running the commands, so that
`plumed`, `gmx_mpi`, and `mpirun` are available.

The topology must already contain the A-state and B-state atom parameters.
The automatic hot region contains complete non-water residues within the
specified distance of perturbed atoms:

```bash
pixi run mdtbx setup_fep_rest \
  -f run.mdp \
  -p hybrid.top \
  -c equilibrated.gro \
  --replicas 32 \
  --temperature 300 \
  --max-temperature 1200 \
  --hot-distance 0.4 \
  -o fep_rest
```

Use `--hot-selection` for an explicit mdtbx atom selection. The setup command
preprocesses the topology, calls `plumed partial_tempering` for every replica,
and creates the final TPR files with `gmx_mpi`.

FEP/REST requires one external MPI rank per replica by default:

```bash
pixi run mdtbx run_fep \
  --path fep_rest \
  --multidir \
  --replex 1000 \
  --ntomp 1
```

The default launcher is `mpirun -np {ntmpi}`, where `ntmpi` defaults to the
selected replica count. A scheduler launcher can be selected explicitly:

```bash
pixi run mdtbx run_fep \
  --path fep_rest \
  --multidir \
  --replex 1000 \
  --mpi-launcher 'srun --ntasks {ntmpi}' \
  --ntomp 1
```

Analyze the state trajectories by rerunning each trajectory with its own and
neighboring Hamiltonians, followed by neighbor BAR:

```bash
pixi run mdtbx analyze_fep_rest --path fep_rest -b 2000
```

PLUMED `partial_tempering` does not support CMAP or intermolecular-interaction
topologies. mdtbx also rejects explicit B-state parameters in `[ pairs ]`,
because PLUMED discards them, and rejects a hot selection that includes only
some copies of a repeated molecule type.

## Generated files

```text
fep/
  fep_manifest.json
  lambda_000/
    fep.mdp
    fep.tpr
  lambda_001/
    fep.mdp
    fep.tpr
  ...
```

After simulation, `analyze_fep` writes `bar.xvg`, `bar_integral.xvg`,
`bar_histogram.xvg`, `bar.log`, and `bar.json`. The JSON file reports the BAR
estimate and uncertainty in both kJ/mol and kcal/mol.

For absolute binding free energy, use the ABFE workflow documented in
`example/abfe/README.md`.
