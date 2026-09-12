Mol* style in PyMOL
===================

The ``molstar_style`` command is provided by the standalone
`molstar_style_in_pymol package <https://github.com/th2ch-g/molstar_style_in_pymol>`_.
The ``pymol_plugins`` distribution installs and registers it automatically.
Both pixi locks record the resolved default-branch commit. The main CLI does not
import the renderer or start PyMOL.

.. code-block:: text

   molstar_style
   molstar_style cartoon, color=secondary-structure
   molstar_style glossy, representation=ball-and-stick
   molstar_style direct-volume, data=density.mrc
   molstar_style ray, filename=figure.png, width=1600, height=1200
   molstar_style reset, name=all

The independent Python implementation supports molecular, volumetric, particle,
measurement, glycan, annotation, orbital, and pairwise-metric visualization.
It reads local data only and has no Node.js, Mol*, CueMol, or service dependency.
Interactive rendering requires PyMOL 3.1 Qt and OpenGL 2.1 / GLSL 1.20;
headless PyMOL supports ray export.

See the `English guide <https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/guide.md>`_,
`Japanese guide <https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/ja/guide.md>`_,
and `coverage table <https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/coverage.md>`_
for input schemas, managed-view restoration, and documented numerical/rendering
approximations. This is not a full Mol* browser application or MVS protocol port.

Ellipsoid views require anisotropic displacement tensors; missing tensors are not replaced by van der Waals spheres. The corrected implementation also checks polymer radii, surface winding, colored cross sections, particle targets, and native ray output.

The `local GPU/ray gallery <https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/gallery.md>`_ and `geometry audit <https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/audit.md>`_ document the examples, regression checks, and remaining approximations. Gallery images are generated locally without CI.

Ribbon curves, cross sections, GGX materials, camera-relative lighting, chain colors and ambient occlusion now follow the pinned Mol* source. The `actual Mol* comparison <https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/fidelity.md>`_ checks matched-camera 8GNG, material and DNA images. Dedicated density ray integrates pixels directly in managed opaque scenes; mixed native scenes retain the documented approximation. All rendering checks run locally without CI.
