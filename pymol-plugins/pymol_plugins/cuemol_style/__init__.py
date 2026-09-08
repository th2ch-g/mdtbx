"""CueMol-inspired molecular geometry and materials for PyMOL 3.1 Qt."""


def cuemol_style(
    style="richardson",
    selection="all",
    representation="auto",
    color="cuemol",
    quality="medium",
    name="cuemol",
    edge="auto",
    edge_width="",
    edge_color="black",
    transparency="keep",
    cache_mb=2048,
    filename="",
    width=0,
    height=0,
    quiet=0,
    _self=None,
):
    """
    DESCRIPTION

        Apply CueMol-inspired geometry and live GPU materials without changing
        PyMOL's source code or its standard commands. Requires PyMOL 3.1 Qt and
        compatibility OpenGL 2.1 / GLSL 1.20. CueMol itself is not required.

    USAGE

        cuemol_style style [, selection [, representation [, color [, quality [, name]]]]]
        cuemol_style list
        cuemol_style refresh [, name=cuemol]
        cuemol_style reset [, name=cuemol]
        cuemol_style png, filename=figure.png [, width=2400, height=1800]
        cuemol_style ray, filename=figure_ray.png [, width=2400, height=1800]

    ARGUMENTS

        style: richardson, ribbon, round_ribbon, fancy_ribbon, cartoon,
            round_cartoon, tube, nucleic, ballstick, cpk, surface, outline,
            silhouette, default, shadow, nolighting, matte, toon1, toon2,
            diff_metal, spec_metal, metallic_chrome, metallic_copper, stone35,
            wood31, wood14scl2. Use list to inspect profiles and active views.
        selection: molecular atoms to style {default: all}
        representation: auto, ribbon, cartoon, tube, nucleic, ballstick,
            sticks, cpk, surface {default: auto; profile may select one}
        color: cuemol, keep, chain, ss, rainbow, element {default: cuemol}
            cuemol uses WoodyHSCPaint for protein ribbons, solid white for
            nucleic ribbons, and DarkCPKColoring for atomic representations.
        quality: low, medium, high {default: medium}
        name: managed group; reapplying replaces it {default: cuemol}
        edge: auto, none, edges, silhouette, thin, normal, thick
        edge_width: positive line width in angstroms (minimum one GPU pixel)
        edge_color: a PyMOL color name {default: black}
        transparency: keep, or 0 (opaque) to 1 (invisible) {default: keep}
        cache_mb: maximum prepared mesh cache in MiB {default: 2048}
        filename: required output PNG path for png/ray
        width, height: output pixels; 0 preserves current dimensions

    NOTES

        All loaded states are prepared before playback. Use refresh after
        coordinate, topology, color, secondary-structure, or state edits.
        Opaque bodies and edges use GLSL. Transparent bodies use native CGO
        with baked lighting. The dedicated ray operation exports temporary
        CGO with approximate materials; standard ray is unchanged and cannot
        see opaque callback geometry. Reset restores native representations
        and background. Use name=all to reset or refresh all managed views.
        Maps, labels, and unsupported representations remain native.
    """
    if _self is None:
        from pymol import cmd as _self
    if not _self.is_gui_thread() and callable(
        getattr(_self, "_call_in_gui_thread", None)
    ):
        parameters = locals().copy()
        return _self._call_in_gui_thread(lambda: cuemol_style(**parameters))
    from .controller import manager_for
    from .presets import PROFILES

    manager = manager_for(_self)
    manager.maintenance()
    style = str(style).strip().lower()
    try:
        if style == "list":
            print(" cuemol_style profiles: " + ", ".join(PROFILES))
            for entry in manager.entries.values():
                states = max(len(ds) for ds in entry.drawings.values())
                print(
                    f" {entry.name}: {entry.options['style']}, {states} states, prepared {entry.seconds:.2f} s, mesh cache {entry.nbytes / 1024**2:.1f} MiB"
                )
            return tuple(PROFILES)
        if style == "reset":
            manager.reset(str(name))
            return
        if style == "refresh":
            manager.refresh(str(name))
            return
        if style in ("png", "ray"):
            from .export import image

            return image(manager, str(filename), width, height, style == "ray")
        entry = manager.apply(
            style,
            str(selection),
            str(representation),
            str(color),
            str(quality),
            str(name),
            str(edge),
            None if edge_width in ("", None) else float(edge_width),
            edge_color,
            None if transparency in ("keep", "", None) else float(transparency),
            cache_mb=float(cache_mb),
        )
        if not int(quiet):
            states = max(len(ds) for ds in entry.drawings.values())
            print(
                f" cuemol_style: {name} ({style}), {states} states prepared in {entry.seconds:.2f} s; {entry.nbytes / 1024**2:.1f} MiB meshes."
            )
        return entry
    except (ValueError, RuntimeError) as exc:
        from pymol import CmdException

        raise CmdException(str(exc)) from exc


def __init_plugin__(app=None):
    from pymol import cmd
    from pymol.shortcut import Shortcut
    from .presets import COLORS, PROFILES, REPRESENTATIONS
    from .controller import manager_for

    cmd.extend("cuemol_style", cuemol_style)
    cmd.auto_arg[0]["cuemol_style"] = [
        Shortcut([*PROFILES, "list", "refresh", "reset", "png", "ray"]),
        "style or operation",
        "",
    ]
    cmd.auto_arg[1]["cuemol_style"] = [cmd.selection_sc, "selection", ""]
    cmd.auto_arg[2]["cuemol_style"] = [Shortcut(REPRESENTATIONS), "representation", ""]
    cmd.auto_arg[3]["cuemol_style"] = [Shortcut(COLORS), "color mode", ""]
    manager_for(cmd)
