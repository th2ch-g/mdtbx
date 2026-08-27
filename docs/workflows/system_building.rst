System preparation
==================

``mdtbx`` provides focused commands for preparing proteins, ligands, solvent,
and GROMACS systems. Inspect each intermediate file before continuing to the
next stage.

Protein and ligand preparation
------------------------------

Typical protein preparation uses ``addh`` followed by optional terminal-cap or
mutation commands:

.. code-block:: console

   $ pixi run mdtbx addh -s protein.pdb -o protein_h
   $ pixi run mdtbx addace -s protein_h.pdb -o protein_ace
   $ pixi run mdtbx addnme -s protein_ace.pdb -o protein_capped

Use ``mutate`` for a single PyMOL mutagenesis operation and
``fill_chainname`` when a PDB file has blank chain identifiers.

For a small molecule, ``gen_am1bcc`` creates GAFF2-compatible AM1-BCC
parameters. ``gen_resp`` performs a Gaussian-based RESP workflow.
``gen_modres_am1bcc`` and ``gen_modres_resp`` target modified residues; the
AM1-BCC modified-residue workflow is currently work in progress.

Build a solvated system
-----------------------

``build_solution`` prepares a solvated all-atom system. Input structures must
already have the intended protonation, residue names, chain identifiers, and
ligand parameters:

.. code-block:: console

   $ pixi run mdtbx build_solution \
       --input protein_h.pdb \
       --ligparam ligand.frcmod:ligand.lib \
       --boxsize 100 100 100 \
       --water-seed 42 \
       --outdir solution

Set ``--water-seed`` when reproducible Packmol water placement is required.
Use ``build_vacuum`` when no solvent, ions, or periodic box should be added.

Copy-ready conventional MD
--------------------------

The repository's ``example/README.md`` catalog identifies maintained and
legacy examples. For a conventional calculation, copy a script from
``example/workflows/`` into the calculation directory and edit its
configuration block. ``solution_setup.sh`` and ``membrane_setup.sh`` prepare
GROMACS inputs with ``mdtbx``; both support a side-effect-free ``--check``
mode.

The companion ``run_slurm.sh`` uses an explicit installer-provided GROMACS
binary for production ``mdrun``. It only runs inside an existing Slurm
allocation and never submits itself. Review the generated topology,
coordinates, index groups, MDP files, and scheduler settings before an
explicit submission.

Place structures
----------------

``place`` applies a random rotation and places two single-chain structures at
a requested distance:

.. code-block:: console

   $ pixi run mdtbx place \
       -1 chain_a.pdb \
       -2 chain_b.pdb \
       -d 30 \
       --max-attempts 20 \
       --seed 42 \
       -o placed.pdb

The seed makes the trial sequence reproducible. Increase ``--max-attempts``
when placement fails because of clashes.

``place_solvent`` uses 3D-RISM solvent distributions. It requires compatible
Amber topology and coordinate files and can retain the generated
susceptibility file:

.. code-block:: console

   $ pixi run mdtbx place_solvent \
       -p leap.parm7 \
       -x leap.rst7 \
       --solvent-model cSPCE \
       --xvv-output cSPCE_300.xvv \
       -o placed_water.pdb

Supporting operations
---------------------

Use the remaining build commands for targeted transformations:

* ``add_ndx``, ``gen_posres``, and ``gen_distres`` generate GROMACS inputs.
* ``amb2gro`` converts Amber files to GROMACS files.
* ``mv_crds_mol2`` replaces MOL2 coordinates while retaining reference data.
* ``calc_ion_conc`` calculates ion counts or concentration.
* ``find_bond`` inspects candidate bonds such as disulfide bonds.
* ``gen_temperatures`` generates temperature ladders for replica exchange.

The canonical shell workflow is under ``example/workflows/``. Older component
examples remain under ``example/build/`` and input MDP templates under
``example/mdp/``; check their status in ``example/README.md`` before use.
