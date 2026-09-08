"""Run real PyMOL checks in a fresh process, outside the suite's PyMOL mocks.

Use --gui for OpenGL, Qt mouse events, and image generation. Output belongs in
an ignored directory. No network access or user startup files are required.
"""

import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import sys
from time import monotonic, perf_counter, sleep
from unittest.mock import patch

import numpy as np


def snapshot(cmd, selection="all"):
    rows = []
    cmd.iterate(
        selection, "rows.append((model,index,reps,color,ss))", space={"rows": rows}
    )
    return rows, cmd.get_coords(selection).copy(), cmd.get_setting_tuple("bg_rgb")


def same_snapshot(cmd, before, selection="all"):
    after = snapshot(cmd, selection)
    assert before[0] == after[0]
    np.testing.assert_array_equal(before[1], after[1])
    assert before[2] == after[2]


def peptide(cmd, structure=None):
    if structure:
        cmd.load(str(structure), "peptide")
    else:
        cmd.fab("ACDEFGHIKLMN", name="peptide", ss=1)
        cmd.alter("peptide and resi 8-10", 'ss="S"')
    cmd.remove("hydro")
    cmd.show_as("cartoon", "peptide")
    cmd.color("salmon", "peptide")
    cmd.bg_color("grey30")
    cmd.orient("peptide")
    cmd.zoom("peptide", 3)
    cmd.set("orthoscopic", 1)


def headless_checks(cmd, output, style, report):
    from pymol import CmdException
    from pymol_plugins.cuemol_style import export, geometry
    from pymol_plugins.cuemol_style.controller import manager_for
    from pymol_plugins.cuemol_style.presets import PROFILES

    native_commands = (cmd.ray, cmd.png, cmd.draw, cmd.do)
    peptide(cmd)
    baseline = snapshot(cmd)
    profiles = []
    for profile in PROFILES:
        entry = style(profile, quality="low", quiet=1, _self=cmd)
        assert all(
            np.isfinite(p.mesh.vertices).all()
            for ds in entry.drawings.values()
            for d in ds
            for p in d.pieces
        )
        style(
            "ray",
            filename=str(output / f"headless_{profile}.png"),
            width=240,
            height=180,
            _self=cmd,
        )
        style("reset", _self=cmd)
        same_snapshot(cmd, baseline)
        profiles.append(profile)
    report["headless_profiles"] = profiles
    style("toon1", quiet=1, _self=cmd)
    manager = manager_for(cmd)
    previous = manager.entries["cuemol"]
    with patch.object(
        geometry, "build", side_effect=RuntimeError("injected mesh failure")
    ):
        try:
            style("toon2", quiet=1, _self=cmd)
        except CmdException:
            pass
        else:
            raise AssertionError("A failed build must raise")
    assert manager.entries["cuemol"] is previous
    enabled = set(cmd.get_names("objects", enabled_only=1))
    with patch.object(
        export, "settings", side_effect=RuntimeError("injected export failure")
    ):
        try:
            style("ray", filename=str(output / "must_not_exist.png"), _self=cmd)
        except CmdException:
            pass
        else:
            raise AssertionError("A failed export must raise")
    assert set(cmd.get_names("objects", enabled_only=1)) == enabled
    assert not (output / "must_not_exist.png").exists()
    assert not any(n.startswith("_cuemol_ray_") for n in cmd.get_names())
    cmd.pseudoatom("cuemol_alpha_1", pos=[100, 100, 100])
    try:
        style("toon1", selection="peptide", transparency=0.4, quiet=1, _self=cmd)
    except CmdException:
        pass
    else:
        raise AssertionError("Generated names must not overwrite source objects")
    assert cmd.count_atoms("cuemol_alpha_1") == 1
    assert manager.entries["cuemol"] is previous
    assert "cuemol" in cmd.get_names("objects")
    cmd.delete("cuemol_alpha_1")
    try:
        style("matte", name="overlap", quiet=1, _self=cmd)
    except CmdException:
        pass
    else:
        raise AssertionError("Overlapping managed selections must be rejected")
    session = output / "session.pse"
    cmd.save(str(session))
    style("reset", _self=cmd)
    cmd.load(str(session))
    assert list(manager.entries) == ["cuemol"]
    style("reset", _self=cmd)
    same_snapshot(cmd, baseline)
    style("richardson", quiet=1, _self=cmd)
    cmd.delete("cuemol_shape_1")
    manager.maintenance()
    same_snapshot(cmd, baseline)
    assert manager.entries == {}
    assert (cmd.ray, cmd.png, cmd.draw, cmd.do) == native_commands
    cmd.fragment("ala", "ligand")
    cmd.show_as("lines", "ligand")
    cmd.hide("everything", "ligand and name CB")
    mixed = snapshot(cmd)
    entry = style("toon1", quiet=1, _self=cmd)
    assert set(entry.drawings) == {"peptide", "ligand"}
    assert entry.drawings["peptide"][0].pieces
    atoms = [a for p in entry.drawings["ligand"][0].pieces for a in p.atoms]
    assert atoms and all(a.name != "CB" for a in atoms)
    style("reset", _self=cmd)
    same_snapshot(cmd, mixed)
    report["headless_lifecycle"] = (
        "apply/reset, failed replacement/export, overlap, saved session, generated-object deletion, unchanged standard commands"
    )
    cmd.delete("all")
    state_checks(cmd, style)
    report["headless_states"] = (
        "nonsequential states, movie mapping, per-object overrides, all states, source enable/disable, added states, transparent state count"
    )


def state_checks(cmd, style, pump=None):
    from pymol_plugins.cuemol_style.controller import manager_for

    manager = manager_for(cmd)
    cmd.fragment("ala", "trajectory")
    cmd.show_as("spheres")
    base = cmd.get_coords("trajectory")
    for state in (2, 3):
        cmd.create("trajectory", "trajectory", 1, state)
        cmd.load_coords(base + [state * 4, 0, 0], "trajectory", state)
    cmd.frame(2)
    style("toon1", representation="cpk", quiet=1, _self=cmd)
    for state in (1, 3, 2):
        cmd.frame(state)
        manager.maintenance()
        if pump:
            pump(0.12)
        assert cmd.get_object_state("cuemol_shape_1") == state
        active = list(manager.active_drawings())
        assert [d.state for d in active] == [state]
        if pump:
            assert active[0].draws > 0
    cmd.mset("1 3 2")
    cmd.frame(2)
    if pump:
        pump(0.12)
    assert cmd.get_state() == 3
    assert [d.state for d in manager.active_drawings()] == [3]
    cmd.set("state", 1, "trajectory")
    manager.maintenance()
    assert [d.state for d in manager.active_drawings()] == [1]
    cmd.unset("state", "trajectory")
    manager.maintenance()
    assert [d.state for d in manager.active_drawings()] == [3]
    cmd.set("all_states", 1, "trajectory")
    manager.maintenance()
    assert [d.state for d in manager.active_drawings()] == [1, 2, 3]
    cmd.unset("all_states", "trajectory")
    cmd.disable("trajectory")
    manager.maintenance()
    assert list(manager.active_drawings()) == []
    cmd.enable("trajectory")
    manager.maintenance()
    cmd.mset("")
    cmd.create("trajectory", "trajectory", 1, 4)
    style("refresh", _self=cmd)
    assert cmd.count_states("cuemol_shape_1") == 4
    style("toon1", representation="cpk", transparency=0.4, quiet=1, _self=cmd)
    assert cmd.count_states("cuemol_alpha_1") == 4
    cmd.frame(4)
    if pump:
        pump(0.12)
    assert cmd.get_object_state("cuemol_alpha_1") == 4
    style("reset", _self=cmd)
    cmd.delete("all")
    cmd.frame(1)


@contextmanager
def qt_session():
    import pymol

    pymol.invocation.options.plugins = 0
    pymol.invocation.options.deferred = []
    pymol.invocation.options.pymolrc = None
    pymol.invocation.options.no_gui = 0
    pymol.invocation.options.show_splash = 0
    from pmg_qt.pymol_qt_gui import PyMOLApplication, PyMOLQtGUI

    app = PyMOLApplication(["cuemol-style-checks"])
    window = PyMOLQtGUI()
    widget = window.pymolwidget
    cmd = widget.cmd

    def in_context(func):
        with widget:
            return func()

    # Match PyMOL's normal execapp context dispatch for this standalone harness.
    cmd._call_with_opengl_context = in_context

    def pump(seconds=0.1):
        end = monotonic() + seconds
        while monotonic() < end:
            app.processEvents()
            sleep(0.002)

    window.resize(1000, 850)
    window.show()
    pump(0.3)
    try:
        yield cmd, widget, pump
    finally:
        manager = vars(cmd._pymol).get("_cuemol_style_manager")
        if manager:
            manager.reset()
        window.hide()
        cmd.delete("all")


def gui_checks(cmd, widget, pump, output, style, report, structure):
    from OpenGL import GL as gl
    from pymol.Qt import QtCore, QtGui, QtWidgets
    from PyQt5.QtTest import QTest
    from pymol_plugins.cuemol_style.controller import manager_for
    from pymol_plugins.cuemol_style.export import view_matrix
    from pymol_plugins.cuemol_style.picking import hit_at
    from pymol_plugins.cuemol_style.presets import PROFILES

    manager = manager_for(cmd)
    with widget:
        report["opengl"] = gl.glGetString(gl.GL_VERSION).decode()
        report["renderer"] = gl.glGetString(gl.GL_RENDERER).decode()
    peptide(cmd, structure)
    cmd.set("internal_gui", 0)
    cmd.set("internal_feedback", 0)
    baseline = snapshot(cmd)
    profiles = []
    for profile in PROFILES:
        entry = style(profile, quality="medium", quiet=1, _self=cmd)
        pump()
        style(
            "png",
            filename=str(output / f"gpu_{profile}.png"),
            width=640,
            height=480,
            _self=cmd,
        )
        assert all(
            d.draws > 0 and not d.error for ds in entry.drawings.values() for d in ds
        )
        if profile == "richardson":
            style(
                "png",
                filename=str(output / "gpu_high_resolution.png"),
                width=2400,
                height=1800,
                _self=cmd,
            )
            style(
                "ray",
                filename=str(output / "ray_richardson.png"),
                width=640,
                height=480,
                _self=cmd,
            )
            cmd.turn("y", 65)
            pump()
            style(
                "png",
                filename=str(output / "gpu_richardson_rotated.png"),
                width=640,
                height=480,
                _self=cmd,
            )
            cmd.turn("y", -65)
        style("reset", _self=cmd)
        same_snapshot(cmd, baseline)
        profiles.append(profile)
    report["gpu_profiles"] = profiles
    cmd.delete("all")
    cmd.pseudoatom("ligand", name="C1", resi="1", resn="LIG", elem="C", pos=[-2, 0, 0])
    cmd.pseudoatom("ligand", name="O1", resi="1", resn="LIG", elem="O", pos=[-2, 1, 0])
    cmd.pseudoatom("ligand", name="N1", resi="2", resn="LIG", elem="N", pos=[2, 0, 0])
    cmd.show_as("spheres", "ligand")
    cmd.color("salmon", "ligand")
    cmd.reset()
    cmd.zoom("ligand", 3)
    entry = style("toon1", representation="cpk", quiet=1, _self=cmd)
    pump(0.3)
    drawing = next(manager.active_drawings())
    np.testing.assert_allclose(drawing.matrices[0], view_matrix(cmd), atol=1e-5)

    def point(world):
        d = next(manager.active_drawings())
        mv, projection, viewport = d.matrices
        clip = projection @ mv @ [*world, 1]
        ndc = clip[:3] / clip[3]
        x = viewport[0] + (ndc[0] + 1) * viewport[2] / 2
        y = viewport[1] + (ndc[1] + 1) * viewport[3] / 2
        return x, y

    def click(world, shift=False):
        x, y = point(world)
        ratio = widget.devicePixelRatioF()
        pos = QtCore.QPoint(round(x / ratio), round(widget.height() - y / ratio))
        QTest.mouseClick(
            widget,
            QtCore.Qt.LeftButton,
            QtCore.Qt.ShiftModifier if shift else QtCore.Qt.NoModifier,
            pos,
        )
        pump(0.2)

    cmd.set("mouse_selection_mode", 0)
    click([-2, 0, 0])
    assert cmd.count_atoms("sele") == 1
    assert cmd.count_atoms("sele and name C1") == 1
    assert len(next(manager.active_drawings()).selected) == 1
    widget.grabFramebuffer().save(str(output / "selection_live.png"))
    cmd.set("mouse_selection_mode", 1)
    click([-2, 0, 0])
    assert cmd.count_atoms("sele") == 2
    click([2, 0, 0], shift=True)
    assert cmd.count_atoms("sele") == 3
    before = np.asarray(cmd.get_view())
    x, y = point([-2, 0, 0])
    start = QtCore.QPoint(
        round(x / widget.devicePixelRatioF()),
        round(widget.height() - y / widget.devicePixelRatioF()),
    )
    QTest.mousePress(widget, QtCore.Qt.LeftButton, pos=start)
    QtWidgets.QApplication.sendEvent(
        widget,
        QtGui.QMouseEvent(
            QtCore.QEvent.MouseMove,
            QtCore.QPointF(start + QtCore.QPoint(90, 40)),
            QtCore.Qt.NoButton,
            QtCore.Qt.LeftButton,
            QtCore.Qt.NoModifier,
        ),
    )
    QTest.mouseRelease(widget, QtCore.Qt.LeftButton, pos=start + QtCore.QPoint(90, 40))
    pump()
    assert not np.allclose(cmd.get_view(), before)
    cmd.set_view(before.tolist())
    cmd.deselect()
    pump()
    # An opaque native sphere in front must block the custom geometry's picker.
    cmd.pseudoatom("occluder", pos=[-2, 0, 4], vdw=2.5)
    cmd.show_as("spheres", "occluder")
    cmd.color("blue", "occluder")
    pump()
    x, y = point([-2, 0, 0])
    depth = manager.events.depth(x, y)
    assert hit_at(list(manager.active_drawings()), x, y, depth) is None
    style(
        "png",
        filename=str(output / "native_occlusion.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    cmd.delete("occluder")
    style("reset", _self=cmd)
    style("toon1", representation="cpk", transparency=0.45, quiet=1, _self=cmd)
    cmd.pseudoatom("native_alpha", pos=[0, 0, 1], vdw=2.7)
    cmd.show_as("spheres", "native_alpha")
    cmd.color("cyan", "native_alpha")
    cmd.set("sphere_transparency", 0.5, "native_alpha")
    pump()
    style(
        "png",
        filename=str(output / "transparency_native.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    cmd.delete("native_alpha")
    style("reset", _self=cmd)
    # Two custom surfaces use PyMOL's native triangle transparency pass together.
    cmd.create("copy", "ligand")
    cmd.translate([1, 0, 2], "copy")
    cmd.zoom("ligand or copy", 3)
    style(
        "toon1",
        selection="ligand",
        representation="surface",
        transparency=0.4,
        name="front_view",
        quiet=1,
        _self=cmd,
    )
    style(
        "matte",
        selection="copy",
        representation="surface",
        transparency=0.55,
        name="back_view",
        color="keep",
        quiet=1,
        _self=cmd,
    )
    cmd.color("blue", "copy")
    style("refresh", name="back_view", _self=cmd)
    pump()
    assert set(manager.entries) == {"front_view", "back_view"}
    assert len(list(manager.active_drawings())) == 2
    style(
        "png",
        filename=str(output / "transparency_surfaces.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    style(
        "ray",
        filename=str(output / "transparency_surfaces_ray.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    style("reset", name="all", _self=cmd)
    cmd.delete("all")
    cmd.fnab("ATGCA", name="dna")
    cmd.show_as("cartoon")
    cmd.orient()
    cmd.zoom("dna", 4)
    style("nucleic", color="chain", quiet=1, _self=cmd)
    pump()
    style("png", filename=str(output / "nucleic.png"), width=640, height=480, _self=cmd)
    style("reset", _self=cmd)
    cmd.delete("all")
    report["gui_interaction"] = (
        "atom/residue/shift selection, native drag rotation, native depth occlusion, native/custom transparency, two transparent surfaces, nucleic slabs"
    )
    state_checks(cmd, style, pump)
    report["gui_states"] = "native callback and CGO state switching plus movie mapping"
    cmd.fragment("ala", "peptide")
    cmd.show_as("spheres")
    baseline = snapshot(cmd)
    entry = style("toon1", quiet=1, _self=cmd)
    pump()
    keys = (
        gl.GL_CURRENT_PROGRAM,
        gl.GL_ARRAY_BUFFER_BINDING,
        gl.GL_ELEMENT_ARRAY_BUFFER_BINDING,
        gl.GL_CLIENT_ACTIVE_TEXTURE,
        gl.GL_FRAMEBUFFER_BINDING,
        gl.GL_DEPTH_WRITEMASK,
        gl.GL_DEPTH_FUNC,
        gl.GL_BLEND,
        gl.GL_CULL_FACE,
        gl.GL_FOG,
    )
    with widget:
        before = [int(gl.glGetIntegerv(key)) for key in keys]
        manager.pool.draw(entry.drawings["peptide"][0])
        assert before == [int(gl.glGetIntegerv(key)) for key in keys]
        assert gl.glGetError() == gl.GL_NO_ERROR
    style("reset", _self=cmd)
    with patch.object(
        manager.pool, "program", side_effect=RuntimeError("injected shader failure")
    ):
        style("toon1", quiet=1, _self=cmd)
        pump(0.3)
    assert not manager.entries
    same_snapshot(cmd, baseline)
    cmd.delete("all")
    report["gui_gl_state"] = (
        "GL state preserved; native view restored after shader failure"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--structure", type=Path)
    parser.add_argument("--output", type=Path, default=Path(".cache/cuemol-checks"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pymol-plugins"))
    try:
        import pymol
    except ImportError:
        print("PyMOL is not installed")
        return 77
    report = {"pymol": pymol.cmd.get_version()[0], "python": sys.version.split()[0]}
    start = perf_counter()
    if args.gui:
        with qt_session() as (cmd, widget, pump):
            from pymol_plugins.cuemol_style import cuemol_style

            gui_checks(
                cmd, widget, pump, args.output, cuemol_style, report, args.structure
            )
            if args.benchmark:
                from cuemol_benchmark import benchmark

                benchmark(
                    cmd, widget, pump, args.output, cuemol_style, report, args.structure
                )
    else:
        from pymol_plugins.cuemol_style import cuemol_style

        headless_checks(pymol.cmd, args.output, cuemol_style, report)
    report["elapsed_seconds"] = perf_counter() - start
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
