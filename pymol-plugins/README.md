# pymol-plugins

PyMOL commands for molecular preparation, analysis, and visualization.
Install the repository's pixi environment and run `pixi run pymolrc` to
configure plugin loading, then restart PyMOL.

`cuemol_style` comes from the standalone
[cuemol_style_in_pymol package](https://github.com/th2ch-g/cuemol_style_in_pymol),
installed from the repository's default branch without an explicit revision.
`pixi.lock` records the resolved version. No sibling checkout is required.
This integration registers its 26 geometry and material presets:

```text
cuemol_style richardson
cuemol_style toon1
cuemol_style matte, representation=surface, transparency=0.4
ray 2400, 1800
png figure_ray.png
cuemol_style png, filename=figure.png, width=2400, height=1800
cuemol_style ray, filename=figure_ray.png, width=2400, height=1800
cuemol_style reset
```

Material and outline styles default to `ribbon`, so helices remain spiral-shaped.
Geometry presets such as `cartoon`, `cpk`, and `surface` select their own geometry.
Use `representation=auto` with a material or outline style to inherit source layers.

Colors follow CueMol GUI defaults: khaki helices, SteelBlue sheets, FloralWhite
coils, and yellow nucleic geometry. Atomic representations use DefaultCPKColoring,
with carbon inheriting the molecular color. `color=keep` preserves existing colors;
`color=ss` selects the optional WoodyHSC palette.
The current background is preserved. Opaque geometry uses GPU shaders and
retained ray-only CGO; transparent geometry uses native CGO. Standard `ray`
and `png, ray=1` work directly. The dedicated ray export adds camera-dependent
outlines and material samples. Ray shading approximates the GPU appearance.
Richardson's colored-pencil strokes are available in opaque GPU rendering and
`cuemol_style png`; ray and transparent CGO use their average tone.
PyMOL's core and standard commands remain unchanged. CueMol is not required.

Requires PyMOL 3.1 Qt and compatibility OpenGL 2.1 / GLSL 1.20 for interactive
rendering; the dedicated ray export also works in headless PyMOL. See the
[full guide](https://github.com/th2ch-g/cuemol_style_in_pymol/blob/main/docs/pymol_cuemol.md)
for presets, selection, trajectories, session restoration, rendering differences,
and validation commands. The implementation, shaders, and dedicated tests live
in that repository. Compatibility imports remain here for older PyMOL sessions.


`molstar_style` comes from the standalone
[molstar_style_in_pymol package](https://github.com/th2ch-g/molstar_style_in_pymol).
It is installed from its default branch; both pixi locks record the resolved commit.
The command supports molecular representations, local volumes and particles,
measurements, SNFG glycans, annotation geometry, orbitals, PAE panels, and rendering
effects. Node.js, Mol*, CueMol, and network services are not runtime dependencies.

```text
molstar_style
molstar_style cartoon, color=secondary-structure
molstar_style glossy, representation=ball-and-stick
molstar_style direct-volume, data=density.mrc
molstar_style ray, filename=molstar_figure.png, width=1600, height=1200
molstar_style reset, name=all
```

See the [English guide](https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/guide.md)
and [Japanese guide](https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/ja/guide.md)
for local input schemas and the explicit differences from Mol* rendering.
The main `mdtbx` CLI does not import this renderer or start PyMOL.

Ellipsoid views require anisotropic displacement tensors. See the
[geometry audit](https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/audit.md)
for corrected defaults and remaining approximations, and the
[local GPU/ray gallery](https://github.com/th2ch-g/molstar_style_in_pymol/blob/main/docs/gallery.md)
for every named style. Gallery rendering runs locally without CI.
