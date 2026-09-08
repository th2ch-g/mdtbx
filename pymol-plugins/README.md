# pymol-plugins

PyMOL commands for molecular preparation, analysis, and visualization.
Install the repository's pixi environment and run `pixi run pymolrc` to
configure plugin loading, then restart PyMOL.

`cuemol_style` provides 26 CueMol-inspired geometry and material presets:

```text
cuemol_style richardson
cuemol_style toon1
cuemol_style matte, representation=surface, transparency=0.4
cuemol_style png, filename=figure.png, width=2400, height=1800
cuemol_style ray, filename=figure_ray.png, width=2400, height=1800
cuemol_style reset
```

Colors default to CueMol's WoodyHSC/DarkCPK palettes; `color=keep` preserves
existing colors. Opaque geometry uses GPU shaders, transparent geometry uses
native CGO, and the dedicated ray export converts the current view to CGO.
PyMOL's core and standard commands remain unchanged. CueMol is not required.

Requires PyMOL 3.1 Qt and compatibility OpenGL 2.1 / GLSL 1.20 for interactive
rendering; the dedicated ray export also works in headless PyMOL. See the
[full guide](../docs/pymol_cuemol.rst) for presets, selection, trajectories,
session restoration, rendering differences, and validation commands.
