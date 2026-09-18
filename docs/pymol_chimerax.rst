ChimeraX style in PyMOL
=======================

The ``chimerax_style`` command is provided by the standalone
`chimerax_style_in_pymol package <https://github.com/th2ch-g/chimerax_style_in_pymol>`_.
The ``pymol_plugins`` distribution installs and registers it automatically.
Both pixi locks record the resolved default-branch commit. No sibling checkout
or ChimeraX installation is required. The main CLI does not load the renderer.

.. code-block:: text

   chimerax_style
   chimerax_style cartoon, color=secondary-structure
   chimerax_style soft, representation=surface
   chimerax_style volume-mesh, selection=density
   chimerax_style ray, filename=figure.png, width=1600, height=1200
   chimerax_style reset, name=all

The map example requires a map object named ``density`` loaded in PyMOL.
Use ``chimerax_style list`` to list the available styles.

The 84 named styles cover atoms, cartoons, surfaces, nucleotides, glycans,
thermal ellipsoids, density maps, annotations, shapes, presets, and lighting.
Managed views preserve source atoms and support trajectories, session save/load,
refresh, and reset. Interactive PyMOL and standard ray export are supported.

See the `English guide <https://github.com/th2ch-g/chimerax_style_in_pymol/blob/main/docs/pymol_chimerax.md>`_,
`Japanese guide <https://github.com/th2ch-g/chimerax_style_in_pymol/blob/main/docs/ja/pymol_chimerax.md>`_,
and `GPU/ray gallery <https://github.com/th2ch-g/chimerax_style_in_pymol/blob/main/docs/gallery.md>`_
for all styles, local data formats, and documented differences from ChimeraX
rendering. The gallery contains interactive GPU and ray images for every named
style.

Molecular surfaces use ``params.grid_spacing`` (default 0.5 Angstrom) at high quality; medium and low use 1.5 and 2 times that spacing. This follows the reference sampling scale, while surface triangulation remains native PyMOL. It reduces rotation and ray cost without changing source colors or material lighting. Set the spacing to 0.125 for the previous dense sampling. Run ``pixi run update`` and restart PyMOL after updating.
