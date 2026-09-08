CueMol-inspired PyMOL styles
============================

The ``cuemol_style`` command is implemented by the standalone
`cuemol_style_in_pymol package <https://github.com/th2ch-g/cuemol_style_in_pymol>`_.
``pymol_plugins`` depends on that Git repository without an explicit revision
and registers the command when imported. ``pixi install`` uses the resolved
version recorded in ``pixi.lock``; no sibling source checkout is required.

Quick start
-----------

Run ``pixi install`` and ``pixi run pymolrc``, then restart PyMOL. After
loading a structure, use the PyMOL command line:

.. code-block:: text

   cuemol_style richardson
   cuemol_style richardson, representation=cpk
   cuemol_style toon1
   cuemol_style matte, representation=surface, transparency=0.4
   ray 2400, 1800
   png figure_ray.png
   cuemol_style png, filename=figure.png, width=2400, height=1800
   cuemol_style ray, filename=outlined_ray.png, width=2400, height=1800
   cuemol_style refresh
   cuemol_style reset
   cuemol_style list

Material and outline styles default to ``ribbon``, so helices remain
spiral-shaped. Geometry presets such as ``cartoon``, ``cpk``, and ``surface``
select their own geometry. Use ``representation=auto`` with a material or
outline style to inherit source layers.

Default colors follow CueMol GUI's initial painting: khaki helices,
SteelBlue sheets, FloralWhite coils, and yellow nucleic geometry. Carbon
inherits the molecular painting. ``color=keep`` uses existing atom colors.
Helix outside faces keep their base color and inside faces use a lighter
color. Applying or resetting styles preserves the existing background.

Standard ``ray`` and ``png, ray=1`` include the custom geometry. The dedicated
``cuemol_style ray`` export adds camera-dependent outlines and material
samples. The ``richardson`` GPU profile uses colored-pencil hatching on warm
paper. Use ``cuemol_style png`` to export those strokes. Ray and transparent
CGO approximate them with average vertex tones. PyMOL's source and standard
commands are unchanged; CueMol is not required.

See the standalone package's
`full guide <https://github.com/th2ch-g/cuemol_style_in_pymol/blob/main/docs/pymol_cuemol.md>`_
for installation requirements, all 26 presets, selection and trajectory
handling, session restoration, memory limits, rendering differences, and
validation commands. Its implementation, shaders, and dedicated tests are
maintained in that repository. mdtbx retains compatibility imports for
older PyMOL sessions that reference ``pymol_plugins.cuemol_style.gpu``.
