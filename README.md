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
  - [ ] T-REMD
  - [ ] 2D-REMD
  - [x] REST
  - [x] REUS
  - [x] US
  - [ ] FEP
  - [ ] FEP/REST
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
  - [ ] BAR
  - [ ] Zwanzig(FEP)
  - [ ] TI
  - [ ] Jarzynski Equality
  - [ ] ERmod

- Analysis
  - [x] trjconv/trjcat
  - [x] fit
  - [x] comdist
  - [x] comvec
  - [x] mindist
  - [x] rmsd
  - [x] rmsf
  - [x] xyz
  - [x] PCA
  - [x] densmap
  - [ ] tICA
  - [ ] RISM/3D-RISM
  - [ ] Elastic Network Model
  - [ ] Normal Mode analysis
  - [ ] Relaxation Mode analysis
  - [ ] Go model
  - [ ] PCA Vector visualization

- Kinetic analysis
  - [x] MSM
  - [ ] TRAM

- Build system
  - [ ] Vacuum
  - [x] Solution
  - [x] Membrane
  - [x] Protein Modeling
  - [x] Modified Residue
  - [x] Make index group
  - [x] Partial Chage
  - [ ] Martini
  - [ ] QM/MM

</details>

## Assumptions
- System build Tools: AMBER, PyMOL, OpenBabel
- Simulation Tools: Gromacs, AMBER, Gaussian16
- ForceField: ff14SB, TIP3P, GAFF2, Lipid21, GLYCAM06-j

## Install
~~~bash
# for pixi users
export PIXI_FROZEN=false
pixi install
pixi run pymolrc
ln -s $PWD/.pixi/envs/default/bin/mdtbx $BIN
ln -s $PWD/.pixi/envs/default/bin/pymol $BIN

# for docker users
docker build -t mdtbx .
~~~

## Update
```bash
pixi run update
```

## Usage
~~~bash
# for pixi users
pixi run mdtbx ...
pixi run gmx ... # equal to pixi run mdtbx cmd gmx ...

# for docker users
docker run -it --rm mdtbx mdtbx ...
docker run -it --rm mdtbx gmx ... # equal to docker run -it --rm mdtbx mdtbx cmd gmx ...
~~~
