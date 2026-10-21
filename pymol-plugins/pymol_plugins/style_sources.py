"""Exclude managed molecular copies from renderer source readers."""

from functools import wraps


def _source_selection(cmd, selection):
    owner = cmd._pymol
    session = owner.session
    generated = set()
    for attribute, key in (
        ("chimerax_style_state", "views"),
        ("blender_style_in_pymol", "entries"),
    ):
        state = getattr(session, attribute, {})
        for entry in state.get(key, {}).values():
            generated.update(entry.get("objects", ()))
    for attribute in ("_cuemol_style_manager", "_molstar_style_manager"):
        manager = getattr(owner, attribute, None)
        if manager is not None:
            for entry in manager.entries.values():
                generated.update(entry.generated)
    existing = set(cmd.get_names("objects"))
    molecules = sorted(
        name for name in generated & existing if cmd.get_type(name) == "object:molecule"
    )
    if not molecules:
        return selection
    excluded = " or ".join("%" + name for name in molecules)
    return f"({selection}) and not ({excluded})"


def _filter_reader(reader, *, method=False):
    if getattr(reader, "_mdtbx_source_filter", False):
        return reader

    @wraps(reader)
    def filtered(owner, selection, *args, **kwargs):
        cmd = owner.cmd if method else owner
        return reader(owner, _source_selection(cmd, selection), *args, **kwargs)

    filtered._mdtbx_source_filter = True
    return filtered


def install():
    """Adapt the pinned source readers, including GUI and refresh entry points."""
    from blender_style_in_pymol.controller import Manager
    from chimerax_style_in_pymol import controller as chimerax
    from cuemol_style_in_pymol import source as cuemol
    from molstar_style_in_pymol import source as molstar

    Manager._source = _filter_reader(Manager._source, method=True)
    chimerax.source_selection = _filter_reader(chimerax.source_selection)
    cuemol.read = _filter_reader(cuemol.read)
    molstar.read = _filter_reader(molstar.read)
