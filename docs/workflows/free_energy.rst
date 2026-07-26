Free-energy workflows
=====================

Standard FEP
------------

The standard GROMACS alchemical workflow has three stages:

1. ``setup_fep`` creates lambda directories, MDP files, TPR files, and a
   manifest.
2. ``run_fep`` runs all or selected lambda windows.
3. ``analyze_fep`` combines the generated derivative files with ``gmx bar``.

The input structure must already be equilibrated. The base MDP defines the
simulation protocol.

.. code-block:: console

   $ pixi run mdtbx setup_fep \
       -f example/mdp/solution/prd.mdp \
       -p gmx.top \
       -c equilibrated.gro \
       --moltype LIG \
       -o fep
   $ pixi run mdtbx run_fep --path fep
   $ pixi run mdtbx analyze_fep --path fep -b 2000

The default decoupling schedule removes electrostatics and then van der Waals
interactions. Use transformation mode only with a topology that already
contains valid A-state and B-state parameters.

FEP/REST
--------

``setup_fep_rest`` prepares a PLUMED FEP/REST Hamiltonian
replica-exchange calculation. It requires PLUMED ``partial_tempering``, a
PLUMED-patched ``gmx_mpi``, MPI, and a compatible dual-state topology.

.. code-block:: console

   $ pixi run mdtbx setup_fep_rest \
       -f run.mdp \
       -p hybrid.top \
       -c equilibrated.gro \
       --replicas 32 \
       --temperature 300 \
       --max-temperature 1200 \
       --hot-distance 0.4 \
       -o fep_rest
   $ pixi run mdtbx run_fep \
       --path fep_rest \
       --multidir \
       --replex 1000 \
       --ntomp 1
   $ pixi run mdtbx analyze_fep_rest --path fep_rest -b 2000

The topology restrictions and installation steps are documented in
``example/fep/README.md``. Review the generated manifest and hot-region
selection before production.

Absolute binding free energy
----------------------------

``setup_abfe`` creates a Boresch-restraint double-decoupling cycle from an
equilibrated receptor-ligand complex and an equilibrated solvated ligand.
Explicit one-based anchor atoms are recommended:

.. code-block:: console

   $ pixi run mdtbx setup_abfe \
       -f run.mdp \
       --complex-topology complex.top \
       --complex-structure complex_equilibrated.gro \
       --solvent-topology ligand_solvent.top \
       --solvent-structure ligand_solvent_equilibrated.gro \
       --moltype LIG \
       --anchors 125 128 131 4201 4204 4207 \
       -o abfe
   $ pixi run mdtbx run_abfe --path abfe --ntomp 1
   $ pixi run mdtbx analyze_abfe --path abfe -b 2000

Inspect ``abfe_manifest.json`` before running production. Charged-ligand
finite-size corrections are not calculated automatically; supply externally
calculated corrections to ``analyze_abfe``. See ``example/abfe/README.md`` for
the complete cycle, result files, symmetry numbers, and corrections.

Reproducibility
---------------

Keep the following with each calculation:

* the exact MDP and topology inputs;
* the generated JSON manifest;
* the ``mdtbx`` version and ``pixi.lock`` revision;
* scheduler, MPI, and thread settings;
* analysis begin time and any correction terms.
