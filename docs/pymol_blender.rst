Blender / Molecular Nodes style in PyMOL
========================================

The ``blender_style`` command is provided by the standalone `blender_style_in_pymol package <https://github.com/th2ch-g/blender_style_in_pymol>`_. The ``pymol_plugins`` distribution installs it and registers both ``blender_style`` and ``blender_load_oxdna`` automatically. Both pixi locks record the resolved default-branch commit. Blender and Molecular Nodes are not runtime dependencies. The main CLI does not import the renderer or start PyMOL.

.. code-block:: text

   blender_style cartoon, color=secondary-structure
   blender_style glossy, representation=ball-and-stick
   blender_style ambient-occlusion, representation=surface
   blender_style density-iso-surface, selection=density
   blender_load_oxdna system.top, trajectory.dat
   blender_style oxdna-ribbon, selection=oxdna
   blender_style ray, filename=figure.png, width=1600, height=1200
   blender_style reset, name=all

The density example requires a loaded map object named ``density``. The oxDNA example requires local topology and trajectory files. Use ``blender_style list`` for all 16 representations and 8 materials, or ``blender_style gui`` for controls.

Managed views preserve source coordinates and colors, support native atom picking, trajectories, session save/load, refresh, and reset. Native molecular copies avoid Python mesh construction during playback; oxDNA base ellipsoids use CGO. Auto quality and a preparation-memory estimate limit the cost of larger structures. The estimate is not a hard process-memory limit.

See the `command guide <https://github.com/th2ch-g/blender_style_in_pymol/blob/main/docs/pymol_blender.md>`_, `GPU/ray gallery <https://github.com/th2ch-g/blender_style_in_pymol/tree/gallery>`_, `reference comparison <https://github.com/th2ch-g/blender_style_in_pymol/blob/main/docs/fidelity.md>`_, and `performance measurements <https://github.com/th2ch-g/blender_style_in_pymol/blob/main/docs/performance.md>`_. Cartoons, surfaces, materials, and ambient occlusion approximate Molecular Nodes. Outlines require ray export; Cycles path tracing and subsurface scattering are not reproduced.
