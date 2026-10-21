"""Exercise all registered renderers in a process without the suite's mocks."""

import subprocess
import sys
import textwrap

import pytest


@pytest.mark.parametrize("mode", ["cartoon", "atoms"])
def test_style_combinations_preserve_original_sources(mode):
    script = textwrap.dedent(
        """
        from itertools import permutations
        import sys

        import numpy as np
        import pymol2
        from blender_style_in_pymol.controller import manager_for as blender_manager
        from chimerax_style_in_pymol.controller import storage
        from cuemol_style_in_pymol.controller import manager_for as cue_manager
        from molstar_style_in_pymol.controller import manager_for as mol_manager

        with pymol2.SingletonPyMOL() as instance:
            from pymol import cmd
            import pymol_plugins

            atomic = sys.argv[1] == "atoms"
            styles = {"chimerax": "sphere", "cuemol": "cpk", "molstar": "spacefill", "blender": "spheres"}
            assert cmd.keyword["blender_load_oxdna"][0] is pymol_plugins.blender_load_oxdna
            for order in permutations(("chimerax", "cuemol", "molstar", "blender")):
                cmd.reinitialize()
                cmd.fab("AAAAAAA", name="sample", ss=1)
                cmd.show_as("cartoon", "sample")
                coordinates = cmd.get_coords("sample").copy()
                before = []
                cmd.iterate("sample", "before.append((reps, color, ss))", space={"before": before})
                for package in order:
                    command = package + "_style"
                    assert cmd.keyword[command][0] is getattr(pymol_plugins, command)
                    if atomic:
                        options = {"quality": "low", "quiet": 1}
                        if package == "cuemol":
                            options["atomic_mode"] = "native"
                        elif package == "molstar":
                            options["params"] = {"atomicMode": "native"}
                        getattr(pymol_plugins, command)(styles[package], **options)
                    else:
                        cmd.do(command + " cartoon, quality=low, quiet=1")
                    if package == "blender":
                        entry = blender_manager(cmd).entries["blender"]
                        keys = cmd.index("%" + entry["source_selection"])
                        assert entry["objects"]
                    elif package == "chimerax":
                        keys = storage(cmd)["views"]["chimerax"]["keys"]
                    else:
                        manager = cue_manager(cmd) if package == "cuemol" else mol_manager(cmd)
                        keys = manager.entries[package].saved
                        assert list(manager.active_drawings())
                        if atomic:
                            assert manager.entries[package].native_objects
                    assert {key[0] for key in keys} == {"sample"}, (order, package)
                    np.testing.assert_array_equal(cmd.get_coords("sample"), coordinates)
                pymol_plugins.blender_style("refresh", name="all", quiet=1)
                keys = cmd.index("%" + blender_manager(cmd).entries["blender"]["source_selection"])
                assert {key[0] for key in keys} == {"sample"}, order
                for package in reversed(order):
                    cmd.do(package + "_style reset, name=all")
                assert cmd.get_names("objects") == ["sample"], order
                after = []
                cmd.iterate("sample", "after.append((reps, color, ss))", space={"after": after})
                assert after == before, order
                np.testing.assert_array_equal(cmd.get_coords("sample"), coordinates)
                print("PASS:", " -> ".join(order))
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", script, mode],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.count("PASS:") == 24


def test_compact_transparent_styles_share_and_restore_native_mode():
    script = textwrap.dedent(
        """
        import pymol2
        from cuemol_style_in_pymol.controller import manager_for as cue_manager
        from molstar_style_in_pymol.controller import manager_for as mol_manager

        with pymol2.SingletonPyMOL() as instance:
            from pymol import cmd
            import pymol_plugins

            for first in ("cuemol", "molstar"):
                cmd.reinitialize()
                cmd.fragment("ala", "sample")
                cmd.create("other", "sample")
                cmd.set("transparency_mode", 2)
                cmd.do("cuemol_style ballstick, selection=sample, transparency=0.4, atomic_mode=native, quiet=1")
                pymol_plugins.molstar_style("ball-and-stick", "other", transparency=0.4, params={"atomicMode": "native"}, quiet=1)
                assert cue_manager(cmd).entries["cuemol"].native_objects
                assert mol_manager(cmd).entries["molstar"].native_objects
                assert cmd.get_setting_int("transparency_mode") == 3
                cmd.do(first + "_style reset, name=all")
                assert cmd.get_setting_int("transparency_mode") == 3
                second = "molstar" if first == "cuemol" else "cuemol"
                cmd.do(second + "_style reset, name=all")
                assert cmd.get_setting_int("transparency_mode") == 2
                assert sorted(cmd.get_names("objects")) == ["other", "sample"]
                print("PASS:", first)
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.count("PASS:") == 2
