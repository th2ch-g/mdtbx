Getting started
===============

Requirements
------------

``mdtbx`` supports Linux x86-64, macOS on Apple silicon, and macOS on Intel.
Python is fixed to version 3.10 in the project environment. Pixi installs most
Python packages and simulation tools from the lock file.

Some commands require software that is not supplied by the default
environment:

* Gaussian 16 for Gaussian-based RESP workflows.
* PLUMED and a PLUMED-patched GROMACS MPI build for FEP/REST.
* ColabFold for ``modeling_cf``.
* Scheduler commands such as ``srun`` when running on a cluster.

The relevant command reference and workflow page states when an external tool
is required.

Install with Pixi
-----------------

Clone the repository, install the locked environment, and configure the PyMOL
plugin path:

.. code-block:: console

   $ git clone https://github.com/th2ch-g/mdtbx.git
   $ cd mdtbx
   $ pixi install
   $ pixi run pymolrc

Run commands through Pixi:

.. code-block:: console

   $ pixi run mdtbx --help
   $ pixi run mdtbx --version
   $ pixi run mdtbx rmsd --help

Set ``PIXI_FROZEN=false`` only when intentionally updating the environment.
Normal reproducible installation should use the committed ``pixi.lock``.

Use Docker
----------

Build the local image:

.. code-block:: console

   $ docker build -t mdtbx .

Run ``mdtbx`` or an installed tool in the container:

.. code-block:: console

   $ docker run --rm -it -v "$PWD:/work" -w /work mdtbx mdtbx --help
   $ docker run --rm -it -v "$PWD:/work" -w /work mdtbx gmx --version

Mount the working directory explicitly so that input and output files persist
after the container exits.

Understand the CLI
------------------

Commands are grouped by purpose in the :doc:`reference/index`. Use ``--help``
at both levels:

.. code-block:: console

   $ pixi run mdtbx --help
   $ pixi run mdtbx build_solution --help

Most structure and trajectory commands accept explicit input and output paths.
Run commands from a dedicated calculation directory so generated files do not
mix with source files or other studies.

Minimal analysis example
------------------------

Fit a trajectory and calculate a contact map:

.. code-block:: console

   $ pixi run mdtbx fit \
       --topology system.tpr \
       --file md.xtc \
       --output fitted.xtc
   $ pixi run mdtbx contactmap \
       --topology system.tpr \
       --trajectory fitted.xtc \
       --output contactmap.npy
   $ pixi run mdtbx show_npy contactmap.npy

Confirm the exact options with ``--help`` because required selections and
accepted topology formats vary by command.
