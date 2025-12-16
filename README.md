# mdtbx
Toolbox for MD simulation

- Build system quickly
- Run conventional/enhanced MD simulation
- Extract CVs from trajectory
- Calculate Free energy
- Run MSM/TRAM analysis

## Assumptions
- System build Tools: AMBER, PyMOL, OpenBabel
- Simulation Tools: Gromacs, AMBER, Gaussian16
- ForceField: ff14SB, TIP3P, GAFF2, Lipid21, GLYCAM06-j

## Install
~~~bash
# for pixi users
pixi install

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
