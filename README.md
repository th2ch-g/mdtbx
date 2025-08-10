# mdtbx
Toolbox for MD simulation

- Build system quickly and run MD simulation
- Extract CVs from trajectory
- Run MSM

## Assumptions
- System build Tools: Amber, PyMOL, OpenBabel
- Simulation Tools: Gromacs
- ForceField: ff14SB, TIP3P, Lipid21, GLYCAM06-j

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
pixi run gmx ...

# for docker users
docker run -it --rm mdtbx pixi run mdtbx ...
docker run -it --rm mdtbx pixi run gmx ...
~~~
