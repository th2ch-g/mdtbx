Trajectory and collective-variable analysis
===========================================

Trajectory preparation
----------------------

Fit or concatenate trajectories before analysis when periodic-boundary
artifacts or separate production segments are present.

.. code-block:: console

   $ pixi run mdtbx trjcat \
       --num_of_step 10 \
       --prefix prd \
       --output fitted.xtc \
       --keep-concatenated

``trjcat`` can retain the concatenated trajectory before the final
periodic-boundary and fitting operations. ``pacs_trjcat`` performs the related
operation for PaCS-MD directory layouts. ``fit`` handles fitting as a separate
step.

Collective variables
--------------------

Commands in the CV category write NumPy arrays that can be inspected with
``show_npy`` or loaded by Python analysis tools:

* ``rmsd`` and ``rmsf`` calculate structural deviations and fluctuations.
* ``comdist`` and ``comvec`` calculate center-of-mass relationships.
* ``mindist`` calculates the minimum distance between selections.
* ``contactmap`` and ``distmap`` calculate residue-level matrices.
* ``densmap`` calculates a two-dimensional density map.
* ``xyz`` extracts coordinates.
* ``pca`` performs principal-component analysis.

Example:

.. code-block:: console

   $ pixi run mdtbx contactmap \
       --topology system.tpr \
       --trajectory fitted.xtc \
       --output contactmap.npy
   $ pixi run mdtbx show_npy contactmap.npy

Selection syntax and defaults differ between commands. Check the command
reference before assuming atom, residue, or mass-weighting behavior.

Extract structures
------------------

Use ``extract_str`` for a specific time and ``extract_ave_str`` for an average
structure. ``show_mdtraj`` prints MDTraj's topology table and helps verify atom
indices and selection fields before analysis.

Performance
-----------

``print_perf`` reads GROMACS log files and prints simulation performance.
``opt_perf`` uses Optuna to compare ``mdrun`` configurations. Performance
results depend on the GROMACS build, hardware, MPI layout, and scheduler
allocation, so repeat optimization on the target machine.

Runnable command examples are available under ``example/cv/`` and
``example/pacs/``.
