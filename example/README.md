# Examples

Use this directory as the starting point for calculation scripts. Copy a
script into the calculation directory, edit the configuration block at the
top, run its `--check` mode, and inspect the generated inputs before starting
the next stage.

## Canonical conventional-MD workflow

The scripts in [`workflows/`](workflows/) are the maintained, copy-ready
entry points:

| Script | Purpose |
| --- | --- |
| [`solution_setup.sh`](workflows/solution_setup.sh) | Build and convert a solvated system with `mdtbx` |
| [`membrane_setup.sh`](workflows/membrane_setup.sh) | Build and convert a membrane system with `mdtbx` and Packmol-Memgen |
| [`run_slurm.sh`](workflows/run_slurm.sh) | Run minimization, equilibration, and segmented production with installer GROMACS |
| [`analyze.sh`](workflows/analyze.sh) | Concatenate trajectories and calculate standard observables with `mdtbx` |

The scripts default to `MDTBX=(mdtbx)`, so `mdtbx` must be installed or on
`PATH`. When using a source checkout without activating its pixi environment,
replace that line with:

```bash
MDTBX=(pixi run --manifest-path /path/to/mdtbx/pyproject.toml mdtbx)
```

For a solution calculation:

```bash
mkdir my-calculation
cd my-calculation
cp /path/to/mdtbx/example/workflows/solution_setup.sh .

# Edit INPUT_PDB, MDP_SOURCE_DIR, box, ions, and optional ligand parameters.
bash solution_setup.sh --check
bash solution_setup.sh

cp /path/to/mdtbx/example/workflows/run_slurm.sh run/
cd run

# Edit GROMACS_BIN and the #SBATCH directives for the target cluster.
bash run_slurm.sh --check
# Review all generated inputs before explicitly submitting:
sbatch run_slurm.sh
```

`solution_setup.sh` expects an Amber-compatible prepared PDB. Confirm
protonation and tautomer states, terminal caps, disulfide bonds, residue and
atom names, and parameters for every nonstandard residue before running it.
Set both ligand parameter paths when the input contains a GAFF2 ligand. Use
`TLEAP_PRE_COMMAND` or `TLEAP_POST_COMMAND` for explicit disulfide,
coordination, or other covalent bonds that tleap cannot infer.

For a membrane calculation, replace `solution_setup.sh` with
`membrane_setup.sh` and edit the lipid, force-field, box, and restraint
settings. The input PDB must already be protonated, named, and oriented for
the intended membrane. The default `PREORIENTED=true` preserves that
orientation by passing `--preoriented` to Packmol-Memgen; set it to `false`
only when Packmol-Memgen should orient the protein with MEMEMBED. Set
`TLEAP_LINE` when Packmol-Memgen needs an explicit tleap bond or other
additional tleap command.

After production finishes:

```bash
cd my-calculation
cp /path/to/mdtbx/example/workflows/analyze.sh .

# Edit segment count and atom selections.
bash analyze.sh --check
bash analyze.sh
```

The setup and analysis scripts call `mdtbx`. The Slurm script requires an
explicit path to installer-provided GROMACS and rejects the project pixi
GROMACS. The pixi GROMACS is intended for preprocessing and analysis; it is
not the production GPU/PLUMED binary. None of the canonical scripts submits a
job, deletes intermediate files, or overwrites an existing setup directory.

Relative paths in the setup and analysis configuration blocks are resolved
from the copied script's directory. `run_slurm.sh` likewise defaults
`SYSTEM_DIR` to its own directory, so `sbatch /path/to/run/run_slurm.sh` and
submission from inside `run/` address the same inputs.

The selection settings use three different syntaxes:

- `CENTER_SELECTION` and `KEEP_SELECTION` are GROMACS index-group names and
  must exist in `index.ndx`.
- `ANALYSIS_SELECTION` and `CONTACT_SELECTION` use MDTraj selection syntax.
- `POSRES_SELECTION` and `LIPID_POSRES_SELECTION` use the `mdtbx`
  atom-selection syntax accepted by `mdtbx gen_posres`. Confirm that the
  lipid selection matches the residue names generated for the selected lipid
  mixture.

The analysis workflow passes `--preserve-time` to `mdtbx trjcat` because the
production segments continue their time from the preceding checkpoint.
The supplied solution and membrane `prd.mdp` templates run at 310 K for 10 ns
per segment. Review these protocol values rather than treating them as
universal defaults.

`rmsd.npy` and `rmsf.npy` are in nm. RMSD is measured against the structure
embedded in the first production TPR. `contactmap.npy` contains the fraction
of frames in contact for each selected representative-atom pair, using the
command's default 6 Å cutoff. Change the script arguments when a different
reference or contact definition is required.

## Restart and overwrite behavior

- Setup scripts require a new `OUTPUT_DIR`. If setup fails, all intermediate
  files remain for diagnosis and the script will refuse to reuse that
  directory. Move the incomplete directory aside or select a new output path
  before retrying.
- `run_slurm.sh` treats `stage.gro` as the completion marker, reuses an
  existing `stage.tpr`, and resumes from `stage.cpt` with `-append`. This makes
  wall-time resubmission safe for unchanged inputs. If an MDP, topology, or
  coordinate file changes, move the old stage outputs aside so a new TPR is
  generated and review the new run separately.
- `analyze.sh` writes to `OUTPUT_DIR` and may replace analysis results with the
  same names. Select a new output directory to retain an earlier analysis.

## Example status

- **Ready**: maintained as a copy-and-edit entry point with validation.
- **Specialized**: useful reference for a narrower method; inspect it and its
  assumptions before adapting it.
- **Site-specific**: contains environment, scheduler, software, or path
  assumptions that must be replaced.
- **Unsupported**: incomplete design notes or placeholders; do not use for a
  calculation.

| Area | Status | Notes |
| --- | --- | --- |
| [`workflows/`](workflows/) | Ready | Canonical solution, membrane, Slurm, and analysis scripts |
| [`abfe/`](abfe/) | Specialized | Absolute binding free energy workflow |
| [`build/`](build/) | Specialized | Older component examples; `coarse_grained/` and `oniom/` are Unsupported |
| [`cv/`](cv/) | Specialized | Individual collective-variable calculations |
| [`fep/`](fep/) | Specialized | FEP and FEP/REST workflows |
| [`gen_distres/`](gen_distres/) | Specialized | Distance-restraint generation |
| [`gen_modres/`](gen_modres/) | Specialized | Modified-residue parameterization |
| [`gen_resp/`](gen_resp/) | Site-specific | Gaussian command and compute environment assumptions |
| [`mbar/`](mbar/) | Specialized | MBAR analysis |
| [`mdp/`](mdp/) | Specialized | MDP templates consumed by the canonical workflows |
| [`mdrun/`](mdrun/) | Site-specific | Legacy scheduler and performance scripts; prefer `workflows/run_slurm.sh` |
| [`modeling/`](modeling/) | Site-specific | ColabFold/runtime environment assumptions |
| [`msm/`](msm/) | Specialized | tICA, clustering, and MSM workflow |
| [`pacs/`](pacs/) | Specialized | PaCS-MD analysis helpers |
| [`place_solvent/`](place_solvent/) | Specialized | 3D-RISM solvent placement |
| [`plumed/`](plumed/) | Specialized | PLUMED input fragments |
| [`remd/`](remd/) | Site-specific | Replica layout and scheduler assumptions |
| [`wham/`](wham/) | Specialized | WHAM analysis |

## Safety boundary

`--check` validates local files and prints the intended work without building
or simulating. It cannot validate scientific choices such as protonation,
force-field compatibility, equilibration length, atom selections, or sampling
convergence. Review those choices and the generated topology, coordinates,
index groups, MDP files, and scheduler resources before submitting a job.
