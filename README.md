# mdtbx
Toolbox for MD simulation

- Build system quickly
- Run conventional/enhanced MD simulation
- Analyze trajectory
- Calculate Free energy

<details> <summary> Supported Features </summary>

- Sampling
  - [x] cMD
  - [x] PaCS-MD
  - [ ] REST
  - [ ] REUS
  - [ ] T-REMD
  - [ ] SMD
  - [ ] US
  - [ ] WT-Metadynamics
  - [ ] OPES
  - [ ] AWH
  - [ ] GaMD
  - [ ] String method
  - [ ] Weighted Ensemble
  - [ ] Adaptive Biasing Force

- Free energy calculation
  - [ ] MMPBSA
  - [ ] MBAR
  - [ ] FEP
  - [ ] TI
  - [ ] Jarzynski Equality

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

- Others
  - [ ] 3D-RISM
  - [ ] PCA Vector visualization

</details>

## Assumptions
- System build Tools: AMBER, PyMOL, OpenBabel
- Simulation Tools: Gromacs, AMBER, Gaussian16
- ForceField: ff14SB, TIP3P, GAFF2, Lipid21, GLYCAM06-j

## Install
~~~bash
# for pixi users
pixi install

# create pymolrc
pixi run pymolrc

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
