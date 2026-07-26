# mdtbx

Toolbox for MD simulation

- Build system quickly
- Run conventional/enhanced MD simulation
- Analyze trajectory
- Calculate Free energy

<details> <summary> Supported Features </summary>

- Sampling

  - [x] cMD
  - [ ] Brownian Dynamics
  - [ ] Langevin Dynamics
  - [ ] Simulated Tempering
  - [ ] Simulated Annealing
  - [x] T-REMD
  - [ ] 2D-REMD
  - [x] REST
  - [x] REUS
  - [x] US
  - [x] FEP
  - [x] FEP/REST
  - [x] AWH
  - [ ] WT-Metadynamics
  - [ ] OPES
  - [ ] GaMD
  - [x] PaCS-MD
  - [ ] SMD
  - [ ] String method
  - [ ] Weighted Ensemble

- Free energy calculation

  - [ ] MMPBSA
  - [x] MBAR
  - [x] WHAM
  - [x] BAR
  - [x] ABFE
  - [ ] Zwanzig(FEP)
  - [ ] TI
  - [ ] Jarzynski Equality
  - [ ] ERmod

- Analysis

  - [x] trjconv/trjcat
  - [x] fit
  - [x] contactmap
  - [x] comdist
  - [x] comvec
  - [x] distmap
  - [x] mindist
  - [x] rmsd
  - [x] rmsf
  - [x] xyz
  - [x] PCA
  - [x] densmap
  - [ ] tICA
  - [x] RISM/3D-RISM
  - [ ] Elastic Network Model
  - [ ] Normal Mode analysis
  - [ ] Relaxation Mode analysis
  - [ ] Go model
  - [x] PCA Vector visualization

- Kinetic analysis

  - [x] MSM
  - [ ] TRAM

- Build system
  - [x] Vacuum
  - [x] Solution
  - [x] Membrane
  - [x] Protein Modeling
  - [x] Modified Residue
  - [x] Make index group
  - [x] Partial Charge
  - [ ] Martini
  - [ ] QM/MM
  - [ ] AutoDock Vina

</details>

## Assumptions

- System build Tools: AMBER, PyMOL, OpenBabel
- Simulation Tools: Gromacs, AMBER, Gaussian16
- ForceField: ff14SB, TIP3P, GAFF2, Lipid21, GLYCAM06-j

## Install

```bash
# for pixi users
export PIXI_FROZEN=false
pixi install
pixi run pymolrc
ln -s $PWD/.pixi/envs/default/bin/mdtbx $BIN
ln -s $PWD/.pixi/envs/default/bin/pymol $BIN

# for docker users
docker build -t mdtbx .

# for developers
pre-commit install
```

## Update

```bash
export PIXI_FROZEN=false; pixi run update
```

## Usage

```bash
# for pixi users
pixi run mdtbx ...
pixi run gmx ... # equal to pixi run mdtbx cmd gmx ...

# for docker users
docker run -it --rm mdtbx mdtbx ...
docker run -it --rm mdtbx gmx ... # equal to docker run -it --rm mdtbx mdtbx cmd gmx ...
```

### Enhanced workflows

```bash
# Build a solution with Packmol-randomized water coordinates.
pixi run mdtbx build_solution -i protein.pdb -o solution \
  --boxsize 100 100 100 --water-seed 42

# Retry deterministic placements until a non-clashing orientation is found.
pixi run mdtbx place -1 chain_a.pdb -2 chain_b.pdb -d 30 \
  --max-attempts 20 --seed 42 -o placed.pdb

# Match MOL2 atoms by atom index when atom names are duplicated.
pixi run mdtbx mv_crds_mol2 -r reference.mol2 -c coordinates.mol2 \
  --match-by index -o updated.mol2

# Choose the final trajectory path and retain the pre-PBC concatenation.
pixi run mdtbx trjcat -n 10 --prefix prd -o fitted.xtc \
  --keep-concatenated

# Generate and retain a reusable solvent susceptibility file.
pixi run mdtbx place_solvent -p leap.parm7 -x leap.rst7 \
  --solvent-model cSPCE --xvv-output cSPCE_300.xvv \
  -o placed_water.pdb

# Prepare, run, and analyze a standard GROMACS decoupling calculation.
pixi run mdtbx setup_fep -f example/mdp/solution/prd.mdp \
  -p gmx.top -c equilibrated.gro --moltype LIG -o fep
pixi run mdtbx run_fep --path fep
pixi run mdtbx analyze_fep --path fep -b 2000

# Prepare and run PLUMED FEP/REST with the installer-provided gmx_mpi.
pixi run mdtbx setup_fep_rest -f run.mdp -p hybrid.top \
  -c equilibrated.gro --replicas 32 -o fep_rest
pixi run mdtbx run_fep --path fep_rest --multidir --replex 1000
pixi run mdtbx analyze_fep_rest --path fep_rest -b 2000
```

See [`example/fep/README.md`](example/fep/README.md) for standard FEP and
FEP/REST, and [`example/abfe/README.md`](example/abfe/README.md) for ABFE.
