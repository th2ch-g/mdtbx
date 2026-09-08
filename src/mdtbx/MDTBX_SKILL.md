# Skill: mdtbx Workflow Operator

Use the `mdtbx` CLI to prepare molecular-dynamics systems, run simulation
workflows, and analyze trajectories and free energies.

**Follow the maintained workflows under `example/` as closely as possible.**
Read `example/README.md` first and prefer scripts marked **Ready**. Adapt an
existing workflow to the task; explain necessary deviations from its protocol.

## When to activate

- Preparing proteins, ligands, solvated systems, or membrane systems for MD.
- Preparing simulation inputs and scheduler scripts, or continuing a run.
- Processing trajectories and calculating structural observables.
- Setting up or analyzing FEP, FEP/REST, ABFE, or MSM calculations.

## Core principles

- Each subcommand performs one function. Compose existing commands and example
  scripts into workflows.
- Use `mdtbx --help` to discover commands and `mdtbx <command> --help` for exact
  flags, input formats, selections, and units. Options differ between commands.
- Work in a dedicated calculation directory. Keep inputs, configuration,
  random seeds, tool versions, and logs with the calculation.
- A successful command or `--check` does not establish scientific validity.
  Review protonation, force fields, equilibration, and sampling separately.

## Environment and example discovery

- With an installed CLI or activated environment, run `mdtbx <command>`.
  From a source checkout, use `pixi install` and `pixi run mdtbx <command>`;
  keep the committed `pixi.lock` for reproducible environments.
- Paths such as `example/workflows/` refer to the source checkout. Set
  `MDTBX_ROOT` to that checkout when working in a separate calculation
  directory. If an installation does not include `example/`, obtain the source
  checkout matching the installed revision and read its catalog. The
  [repository catalog](https://github.com/th2ch-g/mdtbx/blob/main/example/README.md)
  and [documentation](https://th2ch-g.github.io/mdtbx/) provide entry points;
  check their version before adapting a protocol.
- Copied workflow scripts default to `MDTBX=(mdtbx)`. Without an activated
  environment, edit that configuration line to use the checkout:

  ```bash
  MDTBX=(pixi run --manifest-path "$MDTBX_ROOT/pyproject.toml" mdtbx)
  ```

- **GROMACS has two roles.** The pixi GROMACS is for preprocessing and
  analysis. Production GPU or PLUMED `mdrun` requires an explicit
  installer-provided binary. The CLI exposes the project pixi executables on
  `PATH`, so a bare `gmx` within an mdtbx command can select the pixi build.
  Set `GROMACS_BIN` in the Slurm workflow; use `--gmx "$GROMACS_BIN"` for
  commands such as `run_fep`, `run_abfe`, and `analyze_fep_rest` that launch
  production or rerun calculations. Select a compatible MPI/PLUMED build
  when the method requires one.

## Workflow

### Choose and validate a template

- Read the example catalog. **Ready** scripts are maintained entry points;
  **Specialized** examples require method-specific review; **Site-specific**
  examples require adapting their environment assumptions. Do not use
  **Unsupported** placeholders for calculations.
- For conventional MD, use these scripts under `example/workflows/`:

  | Script              | Purpose                                                 |
  | ------------------- | ------------------------------------------------------- |
  | `solution_setup.sh` | Prepare and convert a solvated system                   |
  | `membrane_setup.sh` | Prepare and convert a membrane system                   |
  | `run_slurm.sh`      | Minimize, equilibrate, and run segmented production     |
  | `analyze.sh`        | Process trajectories and calculate standard observables |

- Copy the selected script into the calculation directory and edit its
  configuration block. For `run_slurm.sh`, also edit the scheduler directives.
  Run `bash <script> --check` before executing it. This checks local inputs and
  prints intended work; it does not build or simulate the system.

### Prepare a system

- The setup scripts expect a prepared, Amber-compatible PDB. Resolve
  protonation and tautomer states, caps, disulfides, residue/atom names, and
  parameters for nonstandard residues first. Use commands such as `addh`,
  `gen_am1bcc`, and `gen_resp` as appropriate; check their help and examples.
- For a solution, copy and configure the canonical script:

  ```bash
  cp "$MDTBX_ROOT/example/workflows/solution_setup.sh" .
  # Edit inputs, MDP_SOURCE_DIR, box, ions, and optional ligand parameters.
  bash solution_setup.sh --check
  bash solution_setup.sh
  ```

- Set both ligand parameter files when needed. Use the template's explicit
  tleap configuration for covalent bonds that cannot be inferred.
- For membranes, use `membrane_setup.sh` and review lipid composition,
  parameters, box, and restraints. Its `PREORIENTED=true` default preserves
  the input orientation; set it to `false` when Packmol-Memgen should orient
  the protein with MEMEMBED. Verify the intended orientation and packing.
- Inspect generated coordinates, topology, index groups, restraints, and MDP
  files before continuing. Template temperatures and durations are protocol
  choices to review for the requested calculation.

### Run conventional MD

- Copy `run_slurm.sh` into the generated run directory. Configure the
  installer GROMACS binary, scheduler resources, and production segment count.
- Run `bash run_slurm.sh --check` and inspect the generated inputs. Actual
  execution requires an existing Slurm allocation; the script never submits
  itself.
- Submit scheduler jobs only when the user has explicitly authorized
  submission. Apply any authorization and resource limits already given for
  the task.
- Follow minimization, equilibration, and production logs. Report completed,
  running, pending, and failed stages separately; a queued job has no results.

### Analyze trajectories

- Copy `analyze.sh` next to the run directory. Configure segment count, output
  directory, index groups, and selections; run its `--check` before analysis.
- The canonical workflow uses `trjcat`, `rmsd`, `rmsf`, `contactmap`, and
  `print_perf`. It preserves continued segment times with `--preserve-time`.
  Confirm time continuity, periodic-boundary handling, and matching atom order
  between the processed trajectory and its topology.
- Selection languages are different: `CENTER_SELECTION` and `KEEP_SELECTION`
  are GROMACS index-group names; `ANALYSIS_SELECTION` and `CONTACT_SELECTION`
  use MDTraj syntax; restraint selections passed to `gen_posres` use mdtbx
  atom-selection syntax. Verify the selected atoms or groups.
- RMSD and RMSF outputs are in nm. Review the RMSD reference and contact-map
  definition in `example/README.md`; consult command help for other units and
  formats. Execution success alone does not establish convergence.

## Specialized workflows

Read the corresponding examples and assumptions before composing commands:

| Task                         | Entry point              | Main commands                                                               |
| ---------------------------- | ------------------------ | --------------------------------------------------------------------------- |
| FEP and FEP/REST             | `example/fep/README.md`  | `setup_fep`, `setup_fep_rest`, `run_fep`, `analyze_fep`, `analyze_fep_rest` |
| Absolute binding free energy | `example/abfe/README.md` | `setup_abfe`, `run_abfe`, `analyze_abfe`                                    |
| Kinetic analysis             | `example/msm/`           | `tica`, `cluster`, `msm`                                                    |

Free-energy workflows require compatible equilibrated inputs. Review generated
manifests, lambda schedules, restraints, and any correction terms. For MSM
analysis, preserve independent-trajectory boundaries through the pipeline.
Consult the catalog for REMD, PLUMED, collective variables, and other methods.

## Recovery

- Locate the first failing stage and inspect its command, inputs, and log
  before retrying. Preserve intermediates for diagnosis.
- Setup scripts require a new output directory and retain failed
  intermediates. Choose a fresh output path or move the incomplete directory
  aside before retrying.
- `run_slurm.sh` treats a stage's `.gro` as its completion marker, reuses an
  existing `.tpr`, and resumes from `.cpt` with `-append`. This is suitable for
  unchanged inputs. After changing MDP, topology, or coordinates, move the old
  stage outputs aside so the inputs are regenerated and reviewed.
- `analyze.sh` may replace existing results with the same names. Select a new
  output directory to retain an earlier analysis.
- For a production GPU failure, confirm the selected GROMACS build first.
  The pixi binary is not the production GPU/PLUMED installation.
- Run `mdtbx skill` to reread this guide and `mdtbx <command> --help` to check
  the installed interface before changing a failing invocation.
