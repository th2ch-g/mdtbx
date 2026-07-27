Troubleshooting
===============

Command not found
-----------------

Run the CLI through the project environment:

.. code-block:: console

   $ pixi run mdtbx --help

If installation is incomplete, run ``pixi install`` from the repository root.
Do not use an unrelated Python environment.

PyMOL import errors
-------------------

Run ``pixi run pymolrc`` once after installation. Confirm that the generated
``~/.pymolrc`` refers to the current repository. Regenerate it after moving the
repository.

External executable not found
-----------------------------

Some workflows intentionally depend on separately installed software. Confirm
the expected executable with ``command -v`` and read the relevant workflow
page. FEP/REST requires the PLUMED-patched ``gmx_mpi`` rather than the default
``gmx`` binary.

GROMACS: no GPU is detected
---------------------------

.. code-block:: text

   Fatal error:
   Cannot run short-ranged nonbonded interactions on a GPU because no GPU is
   detected.

The bundled GROMACS is the conda-forge build, which reports
``GPU support: OpenCL`` and does not detect NVIDIA GPUs. It is intended for
preprocessing and analysis (``grompp``, ``trjconv``, ``editconf``), not for
production ``mdrun``. Run production MD with a GROMACS built by
``install_scripts/gmx*.sh``.

Note that mdtbx prepends its own environment to ``PATH``, so a bare ``gmx``
resolves to the bundled build even when another GROMACS was exported first.
Pass the intended binary explicitly with ``--gmx`` (``run_fep``, ``run_abfe``,
``analyze_fep_rest``) or spell it out in ``--mdrun-command`` (``opt_perf``).

Unexpected atom selections
--------------------------

Use ``show_mdtraj`` to inspect topology fields and atom indices. Verify whether
the command expects zero-based MDTraj selections, one-based GROMACS indices,
or another explicitly documented format.

Generated files overwrite prior work
------------------------------------

Run each workflow in a dedicated output directory. Inspect ``--output`` and
``--path`` before execution. Preserve manifests and inputs in versioned or
archived calculation directories.

Get more diagnostic output
--------------------------

Re-run the command with the same inputs and capture both standard output and
standard error. Record the command, ``mdtbx --version``, platform, and external
tool versions when reporting a problem.
