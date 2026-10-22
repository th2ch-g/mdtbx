Blender / Molecular Nodes style in PyMOL
========================================

The ``blender_style`` command is provided by the standalone `blender_style_in_pymol package <https://github.com/th2ch-g/blender_style_in_pymol>`_. The ``pymol_plugins`` distribution installs it and registers both ``blender_style`` and ``blender_load_oxdna`` automatically. Both pixi locks record the resolved default-branch commit. Blender 4.0+ is an optional runtime dependency for Cycles refraction and subsurface scattering. Ordinary materials work without Blender, and Molecular Nodes is not required. The main CLI does not import the renderer or start PyMOL.

.. code-block:: text

   blender_style cartoon, color=secondary-structure
   blender_style glossy, representation=ball-and-stick
   blender_style ambient-occlusion, representation=surface
   blender_style density-iso-surface, selection=density
   blender_load_oxdna system.top, trajectory.dat
   blender_style oxdna-ribbon, selection=oxdna
   ray 1600, 1200
   png figure.png
   blender_style reset, name=all

The density example requires a loaded map object named ``density``. The oxDNA example requires local topology and trajectory files. Use ``blender_style list`` for all 16 representations and 10 materials, or ``blender_style gui`` for controls.

Managed views preserve source coordinates and colors and support trajectories, session save/load, refresh and reset. Automatic rendering uses smooth meshes for supported static selections within the memory estimate, and native PyMOL for trajectories and larger selections. Use ``backend=mesh`` to request meshes explicitly or ``backend=native`` for native atom picking. Mesh surfaces do not carry atom picking IDs. Mesh sessions require the package to be installed when reopened.

The mesh backend provides relaxed chain-separated surfaces, oriented protein cartoons, rectangular nucleic bases, geometric occlusion and linear-light GGX materials. Geometry is prepared once and stored in indexed GPU buffers with a bounded cache. Material changes reuse meshes. Explicit mesh trajectories store all selected states; automatic trajectories use native state updates. The preparation estimate is not a hard process-memory limit.

Standard ``ray``, Python ``cmd.ray()`` and ``png ..., ray=1`` automatically refresh current-state mesh shading for the camera. The hooks preserve the completed ray image for subsequent ``png`` and restore temporary settings. The ``blender_style ray`` helper additionally checks image-buffer estimates. GPU and ray outputs support approximate Fresnel outlines and layered transparency for physical meshes.

The ``glass`` and ``subsurface`` materials automatically render standard ``ray`` and ``png ..., ray=1`` commands through Blender Cycles. The viewport uses lightweight approximations. Blender is discovered on ``PATH`` or through macOS Spotlight; ``BLENDER_EXECUTABLE`` selects an executable explicitly. Cycles requires a complete managed mesh scene and reports visible unsupported native objects. Use ``samples=128`` for more samples, ``renderer=cycles`` for other materials, or ``renderer=pymol`` for the preview approximation. The ``params`` object accepts ``ior``, ``roughness``, ``transmission``, ``subsurface`` and ``subsurface_radius`` (in scene angstroms). See the `Cycles guide <https://github.com/th2ch-g/blender_style_in_pymol/blob/main/docs/cycles.md>`_.

.. code-block:: text

   blender_style glass, representation=surface
   ray 1200, 900
   png glass.png
   blender_style subsurface, representation=surface, samples=128
   png subsurface.png, 1200, 900, ray=1

See the `command guide <https://github.com/th2ch-g/blender_style_in_pymol/blob/main/docs/pymol_blender.md>`_, `GPU/ray gallery <https://github.com/th2ch-g/blender_style_in_pymol/tree/gallery>`_, `reference comparison <https://github.com/th2ch-g/blender_style_in_pymol/blob/main/docs/fidelity.md>`_, and `performance measurements <https://github.com/th2ch-g/blender_style_in_pymol/blob/main/docs/performance.md>`_. Molecular Nodes geometry and materials remain approximations. Mixed native/mesh ray scenes retain scene lighting and can differ from the dedicated GPU shader. Cycles computes actual refraction and random-walk subsurface scattering on the exported meshes; it does not replace the plugin geometry with Molecular Nodes node graphs.
