mdtbx
=====

``mdtbx`` is a command-line toolbox for molecular dynamics workflows. It
combines small, focused commands for system preparation, simulation support,
trajectory processing, collective-variable analysis, and free-energy
calculations.

The project works with established molecular-simulation tools rather than
replacing them. Common integrations include GROMACS, AMBER Tools, PyMOL,
Open Babel, Gaussian, and PLUMED.

Start here
----------

* :doc:`getting_started` explains installation and the command-line model.
* :doc:`workflows/system_building` covers common system-preparation steps.
* :doc:`workflows/trajectory_analysis` shows trajectory and CV pipelines.
* :doc:`workflows/free_energy` covers FEP, FEP/REST, and ABFE.
* :doc:`agent_workflows` describes approved autonomous Agent workflows.
* :doc:`reference/index` lists every available command and option.

Design principle
----------------

Each subcommand performs one focused operation. Compose commands into a
pipeline instead of expecting one command to perform an entire study:

.. code-block:: console

   $ pixi run mdtbx addh --structure protein.pdb --output_prefix protein_h
   $ pixi run mdtbx build_solution --input protein_h.pdb --outdir solution
   $ pixi run mdtbx fit --file md.xtc --topology system.tpr --output fitted.xtc
   $ pixi run mdtbx contactmap --topology system.tpr \
       --trajectory fitted.xtc --output contactmap.npy

.. toctree::
   :maxdepth: 2
   :hidden:

   getting_started
   workflows/system_building
   workflows/trajectory_analysis
   workflows/free_energy
   agent_workflows
   pymol_ai
   troubleshooting
   reference/index
   development
   documentation
