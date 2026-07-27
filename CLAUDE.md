# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code)
when working with code in this repository.

## Overview

`mdtbx` is a CLI toolbox for molecular dynamics simulations. It supports
system preparation, simulation execution, trajectory analysis, and free-energy
calculations.

External tools: AMBER, PyMOL, Open Babel, GROMACS, Gaussian 16
Force fields: ff14SB, TIP3P, GAFF2, Lipid21, GLYCAM06-j

## Design principles

Keep every subcommand small and focused on **one function**. Build complex
workflows by combining multiple subcommands.

```bash
# Good: combine small subcommands to build a pipeline
pixi run mdtbx addh --input protein.pdb --output protein_h.pdb
pixi run mdtbx gen_am1bcc --input ligand.mol2 --output ligand.mol2
pixi run mdtbx build_solution --receptor protein_h.pdb --ligand ligand.mol2
pixi run mdtbx fit --trajectory md.xtc --output fitted.xtc
pixi run mdtbx rmsd --trajectory fitted.xtc --output rmsd.npy
```

Before adding a subcommand, first consider whether the feature can be built by
combining existing commands.

## Development commands

```bash
# Set up the environment
pixi install

# Run the CLI
pixi run mdtbx <subcommand>
pixi run gmx ...           # Run the GROMACS binary from the pixi environment
                           # (preprocessing and analysis only -- see below)

# Test
pixi run test              # Run the full test suite
pixi run test-fast         # Stop after the first failure (-x)

# Format and lint
pixi run r                 # Run ruff format and ruff check
pixi run ruff-format       # Format only
pixi run ruff-lint         # Lint only
pre-commit run --all-files # Run all pre-commit hooks

# Update
pixi run update            # Run git pull and pixi install

# Configure PyMOL
pixi run pymolrc           # Generate ~/.pymolrc

# JupyterLab on a remote host
pixi run jupyter_remote
```

## Architecture

```text
src/
  __main__.py    # Entry point: main() -> cli()
  cli.py         # Automatic subcommand discovery and dispatch via pkgutil
  config.py      # Global constants such as water density, Gaussian settings, and MAXWARN
  logger.py      # Logger factory
  utils/         # General utilities and commands (mod_mdp, convert, rmfile, cmd, shell_hook,
                 # show_mdtraj, show_npy, partial_tempering, run_fep, run_abfe)
                 # Shared libraries without add_subcmd are not registered as commands:
               #   atom_selection_parser, parse_top, proc(run_cmd), tleap(run_tleap),
                 #   gmx(gmx_index_flag/to_gmx_index/gmx_tempfile), pymol_session,
                 #   common_args, fep, fep_rest, abfe
  build/         # System-building commands (addace, addh, addnme, add_ndx, mv_crds_mol2,
                 # calc_ion_conc, centering_gro, find_bond, gen_am1bcc, gen_resp,
                 # gen_modres_am1bcc, gen_modres_resp, gen_posres, gen_distres,
                 # modeling_cf, amb2gro, build_solution, build_vacuum, place_solvent,
                 # gen_temperatures, mutate, setup_fep, setup_fep_rest, setup_abfe)
  trajectory/    # Trajectory commands (fit, trjcat, pacs_trjcat, print_perf, opt_perf)
  analysis/      # Analysis commands (extract_str, extract_ave_str, analyze_fep,
                 # analyze_fep_rest, analyze_abfe)
  cv/            # Collective-variable commands (comdist, comvec, densmap, mindist, rmsd,
                 # rmsf, xyz, pca, contactmap, distmap)

tests/
  conftest.py      # Shared fixtures and PyMOL mocks
  fixtures/        # Test data such as sample.mdp, sample.top, and sample.pdb
  test_utils/      # Tests for src/utils/
  test_build/      # Tests for src/build/
  test_trajectory/ # Tests for src/trajectory/
  test_analysis/   # Tests for src/analysis/
  test_cv/         # Tests for src/cv/
  test_cli.py      # CLI registration tests for all subcommands

pymol-plugins/
  pymol_plugins/ # PyMOL plugins such as builders, visualizers, and selectors

example/         # Example notebooks and scripts organized by use case
install_scripts/ # Manual installation scripts for GROMACS, PLUMED, and related tools
```

### Subcommand pattern

Each command module under `src/build/`, `src/trajectory/`, `src/analysis/`,
`src/cv/`, or `src/utils/` implements these two functions:

```python
def add_subcmd(subparsers):
    # Define the argparse subcommand and its arguments

def run(args):
    # Implement the command
```

At startup, `cli.py` scans the category packages
(`build`/`trajectory`/`analysis`/`cv`/`utils`) and **automatically registers**
modules that expose `add_subcmd`. Adding a command only requires placing its
module in the correct directory; `cli.py` does not need to change. Libraries
and helpers without `add_subcmd` are skipped.

Import shared helpers and parsers through `..utils`:

```python
from ..utils.atom_selection_parser import AtomSelector
from ..utils.parse_top import GromacsTopologyParser
from ..utils.proc import run_cmd                 # Run subprocesses and log success
from ..utils.tleap import run_tleap              # Write tleap.in, run tleap, and clean up
from ..utils.gmx import gmx_index_flag, to_gmx_index, gmx_tempfile
from ..utils.pymol_session import pymol_session  # Reinitialize PyMOL and load a structure
from ..utils.common_args import add_topology_arg, add_trajectory_arg, add_output_arg
```

### Configuration (`src/config.py`)

- `MAXWARN`: maximum number of warnings accepted by `grompp`
- `GAUSSIAN_CMD`, `STRUCTURE_OPTIMIZATION`, and
  `SINGLE_POINT_CALCULATION`: Gaussian settings
- Density and molecular-volume constants for TIP3P, TIP4P, TIP5P, and OPC water
  models
- Side effects that add `.pixi/envs/default/bin` to `PATH` and load
  `pymol_plugins` are isolated in `src/__init__.py`; `config.py` contains
  constants only.

## GROMACS: two installations, two roles

| installation | role |
|---|---|
| `pixi run gmx` (conda-forge) | preprocessing and analysis: `grompp`, `trjconv`, `editconf`, ... |
| `$TOOLS/gromacs/*` (built by `install_scripts/gmx*.sh`) | production `mdrun`, including GPU and PLUMED runs |

The conda-forge build reports `GPU support: OpenCL` and does not detect NVIDIA
GPUs, so `pixi run gmx mdrun -nb gpu` fails with "no GPU is detected". This is
expected: production MD is not meant to run through the pixi GROMACS.

`src/__init__.py` prepends the pixi `bin` to `PATH` unconditionally, so a bare
`gmx` inside mdtbx resolves to the pixi build even when the caller exported
their own GROMACS first. Point the mdrun-running subcommands at the intended
binary explicitly:

```bash
mdtbx run_fep          --gmx $TOOLS/gromacs/2025.1/gromacs-2025.1/bin/gmx ...
mdtbx run_abfe         --gmx $TOOLS/gromacs/2025.1/gromacs-2025.1/bin/gmx ...
mdtbx analyze_fep_rest --gmx $TOOLS/gromacs/2022.5-mpi-plumed/gromacs-2022.5/bin/gmx_mpi ...
mdtbx opt_perf --mdrun-command "$TOOLS/gromacs/2025.1/gromacs-2025.1/bin/gmx mdrun -deffnm prd"
```

## Environment management

- Package management: `pixi` with both conda and pip dependencies
- Python version: fixed at 3.10
- Reproducibility: maintained through `pixi.lock`
- Container support: build with `Dockerfile`
