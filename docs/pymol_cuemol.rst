CueMol-inspired PyMOL styles
========================================

``cuemol_style`` adds molecular geometry and live GPU materials to the
``pymol_plugins`` package. It requires PyMOL 3.1 with Qt, NumPy, SciPy, and
PyOpenGL, and a compatibility OpenGL 2.1 / GLSL 1.20 context. ``pixi install``
installs the Python dependencies. CueMol itself is not required. PyMOL's
source code and standard commands are unchanged.

Quick start
-----------

Configure the plugin with ``pixi run pymolrc`` and restart PyMOL, or import
``pymol_plugins`` in an already configured Python environment. Run these
commands in the PyMOL command line after loading a structure:

.. code-block:: text

   cuemol_style richardson
   cuemol_style toon1
   cuemol_style toon2, selection=chain A
   cuemol_style matte, representation=surface, transparency=0.4
   cuemol_style list
   help cuemol_style

The default named view is ``cuemol``. Applying another style with the same
name replaces that view. ``richardson`` uses thin ribbons, sheet arrows,
lighter undersides, stepped lighting, and black silhouette/crease lines.
Rotation and zoom update the GPU rendering immediately.

Styles and controls
-------------------

.. list-table:: Style names
   :header-rows: 1
   :widths: 25 75

   * - Group
     - Names
   * - Protein geometry
     - ``richardson``, ``ribbon``, ``round_ribbon``, ``fancy_ribbon``,
       ``cartoon``, ``round_cartoon``, ``tube``
   * - Other molecular geometry
     - ``nucleic``, ``ballstick``, ``cpk``, ``surface``
   * - Lighting and shading
     - ``default``, ``shadow``, ``nolighting``, ``matte``, ``toon1``, ``toon2``
   * - Decorative materials
     - ``diff_metal``, ``spec_metal``, ``metallic_chrome``, ``metallic_copper``,
       ``stone35``, ``wood31``, ``wood14scl2``
   * - Outlines
     - ``outline``, ``silhouette``

Geometry presets choose a representation. Material presets use
``representation=auto``: existing sticks, spheres, surface, cartoon, and
ribbon layers are retained as custom geometry. Atoms shown only as lines or
nonbonded points become protein ribbons, nucleic backbones/base slabs, or
ball-and-stick geometry. Hidden atoms stay hidden in this automatic mode.
Maps, labels, and other unsupported layers remain native.

Override the representation with ``ribbon``, ``cartoon``, ``tube``,
``nucleic``, ``ballstick``, ``sticks``, ``cpk``, or ``surface``. ``cartoon``
uses helix cylinders; ``ribbon`` follows the backbone. Quality is ``low``,
``medium`` (default), or ``high``. Higher quality increases preparation time
and memory use.

``edge`` accepts ``auto``, ``none``, ``edges``, ``silhouette``, ``thin``,
``normal``, or ``thick``. ``edges`` includes sharp visible creases.
``edge_width`` is in angstroms, with a one-pixel minimum in GPU images;
``edge_color`` accepts a PyMOL color. ``transparency=keep`` preserves the
source representation's opacity; a number from 0 to 1 overrides it.

.. code-block:: text

   cuemol_style spec_metal, representation=cpk, edge=silhouette
   cuemol_style richardson, edge_width=0.12, edge_color=black
   cuemol_style wood31, representation=surface, quality=high

Default colors
--------------

``color=cuemol`` is the default. It combines CueMol's WoodyHSCPaint palette
for protein ribbons/cartoons/tubes with DarkCPKColoring for atomic and
surface representations on a white background. Nucleic backbone/base
geometry is white. This is a deliberate combination of CueMol styles.

.. list-table:: CueMol palette
   :header-rows: 1

   * - Feature
     - Color
   * - Helix / sheet / coil
     - ``#FF7F7F`` / ``#7FFF7F`` / ``#FFFFBF``
   * - Carbon / nitrogen / oxygen
     - ``#404040`` / blue / red
   * - Hydrogen / sulfur / phosphorus
     - cyan / green / yellow

Other modes are ``keep`` (existing atom colors), ``chain``, ``ss``
(WoodyHSC secondary-structure colors), ``rainbow``, and ``element`` (a
conventional element palette). Changing a PyMOL atom color while a style is
active takes effect after ``refresh`` and only when the color mode uses that
atom color. Source atom colors and secondary-structure assignments are
never changed by this command.

.. code-block:: text

   cuemol_style richardson, color=keep
   cuemol_style nucleic, color=chain

Selection and state management
------------------------------

Click custom geometry to select its original atoms in ``sele``. PyMOL's
mouse selection mode controls atom/residue/chain expansion; Shift adds to
the selection. Pink markers identify selected source atom positions.
Dragging retains native rotation/movement. Native opaque geometry blocks
clicks on custom geometry behind it. Picking requires the Qt GUI.

All loaded states are prepared before playback. Global state changes,
movie frame-to-state mappings, object state overrides, and ``all_states``
are supported. Coordinate, bond, color, secondary-structure, transparency,
or state-count edits require ``refresh``. To change the native representation
layers, reset the view, change the layers, and apply the style again. Source
visibility and object state settings are synchronized by a short Qt timer.

.. code-block:: text

   cuemol_style refresh
   cuemol_style reset
   cuemol_style reset, name=all

``reset`` restores the original native representation masks and background,
and removes generated objects and GUI hooks. Source coordinates and bonds
are untouched. Use distinct ``name`` values for disjoint selections;
overlapping managed selections are rejected. Deleting a managed group or
source object triggers cleanup. Saved PyMOL sessions store recipes and
rebuild the views on load when the plugin is available.

``cache_mb=2048`` limits prepared mesh storage and separately guards the
estimated input snapshots. It is not a cap on total process memory. The
GPU cache retains up to 256 MiB of vertex buffers and evicts old states;
evicted buffers are uploaded again from prepared CPU meshes. Large
trajectories, especially surfaces and transparent atomistic views, can be
expensive to prepare. ``cuemol_style list`` reports preparation time and
mesh storage for each active view.

Image export
------------

.. code-block:: text

   cuemol_style png, filename=figure.png, width=2400, height=1800
   cuemol_style ray, filename=figure_ray.png, width=2400, height=1800

Both operations save the complete visible scene and require a filename.
``png`` captures the GPU appearance and requires the Qt GUI. ``ray`` builds
temporary CGO geometry for the current state and camera; it also works in
headless PyMOL. Export restores visibility, playback, and temporary
settings even if rendering fails. Selection markers are omitted from
exports.

Opaque bodies use independent GPU shaders. Transparent bodies use native
CGO with lighting baked into vertex colors so that PyMOL can composite
them with native translucent objects. Their lighting is fixed in molecular
coordinates during rotation. Ray output samples materials into vertex
colors and uses PyMOL lighting and cylindrical outline geometry, so its
shading and line appearance differ from the GPU image. Wood, stone, and
metal are procedural approximations; they do not reproduce CueMol's
POV-Ray textures exactly. ``shadow`` is a flat shading material, not a
scene-shadow generator.

Standard PyMOL ``ray`` and ``png, ray=1`` cannot see opaque callback
geometry. Use the dedicated ``cuemol_style ray`` operation for a complete
image. Object transformation matrices, stereo/VR picking, editing atoms
by dragging, and headless interactive GPU rendering are outside the
supported interface; apply coordinate transforms to source atoms and
refresh when needed.

Validation
----------

The standalone harness exercises every style in real PyMOL, restoration
after failures, multiple states, session reload, Qt picking, rotation, and
native/custom transparency. GUI mode also writes GPU and ray images for
visual review. Run Python through the repository's pixi interpreter:

.. code-block:: console

   $ uv run --no-project --python .pixi/envs/default/bin/python python \
       tests/test_pymol_plugins/check_cuemol_style.py --output .cache/cuemol-headless
   $ uv run --no-project --python .pixi/envs/default/bin/python python \
       tests/test_pymol_plugins/check_cuemol_style.py --gui --benchmark \
       --output .cache/cuemol-gui

The benchmark uses 500 residues and 100 synthetic states at a 1280 by 720
viewport with medium-quality ribbons. It reports preparation time, CPU/GPU
mesh storage, rotation speed, explicit state-switch speed, actual movie
draw rate, and one-state surface preparation. The targets after preparation
are 30 FPS rotation and 15 FPS playback; results depend on the GPU and input
geometry. Generated images and reports under ``.cache`` are ignored by Git.

The geometry and names are inspired by the
`CueMol ribbon renderer <https://cuemol.github.io/cuemol2_docs/cuemol2/RibbonRenderer/>`_
and CueMol style definitions. This plugin implements its own mesh generation
and rendering and has no runtime dependency on the CueMol source tree.
