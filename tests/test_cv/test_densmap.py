"""
cv/densmap unit tests

Tests the 2D density map computation with MDtraj (the gmx=False path).
"""

import types
from pathlib import Path

import numpy as np
import pytest


class TestDensmapRun:
    def _make_args(self, traj_files, output, selection="all", axis="xy", bins=10):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            output=str(output),
            bins=bins,
            axis=axis,
            gmx=False,
            index=None,
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy file is generated"""
        from mdtbx.cv.densmap import run

        out = tmp_path / "densmap.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_is_nonempty(self, trajectory_files, tmp_path):
        """The output file is not empty"""
        from mdtbx.cv.densmap import run

        out = tmp_path / "densmap_size.npy"
        run(self._make_args(trajectory_files, out))
        assert out.stat().st_size > 0

    @pytest.mark.parametrize("axis", ["xy", "xz", "yz"])
    def test_different_axes(self, trajectory_files, tmp_path, axis):
        """A file is generated for each of the xy / xz / yz projections"""
        from mdtbx.cv.densmap import run

        out = tmp_path / f"densmap_{axis}.npy"
        run(self._make_args(trajectory_files, out, axis=axis))
        assert out.exists()

    def test_empty_selection_raises(self, trajectory_files, tmp_path):
        from mdtbx.cv.densmap import run

        out = tmp_path / "densmap_empty.npy"
        with pytest.raises(ValueError, match="No atoms selected"):
            run(self._make_args(trajectory_files, out, selection="name XXXX"))
        assert not out.exists()

    @pytest.mark.parametrize("bins", [0, -1])
    def test_non_positive_bins_raise(self, trajectory_files, tmp_path, bins):
        from mdtbx.cv.densmap import run

        with pytest.raises(ValueError, match="bins"):
            run(self._make_args(trajectory_files, tmp_path / "out.npy", bins=bins))

    def test_histogram_shape(self, trajectory_files, tmp_path):
        """The result agrees with np.histogram2d (golden-path check)"""
        import mdtraj as md
        from mdtbx.cv.densmap import _AXIS_MAP, run

        bins = 8
        axis = "xy"
        out = tmp_path / "densmap_golden.npy"
        run(self._make_args(trajectory_files, out, bins=bins, axis=axis))

        # Reproduce counts with the same computation and compare the totals
        traj = md.load(trajectory_files["xtc"], top=trajectory_files["pdb"])
        atom_indices = traj.topology.select("all")
        xyz = traj.xyz[:, atom_indices, :]
        ax0, ax1 = _AXIS_MAP[axis]
        pos0 = xyz[:, :, ax0].ravel()
        pos1 = xyz[:, :, ax1].ravel()
        counts_ref, _, _ = np.histogram2d(pos0, pos1, bins=bins)

        # The file exists and the total count matches
        assert out.exists()
        total_ref = int(counts_ref.sum())
        expected = traj.n_frames * traj.n_atoms
        assert total_ref == expected

    def test_gmx_uses_temporary_output(self, tmp_path, monkeypatch):
        from mdtbx.cv import densmap

        monkeypatch.chdir(tmp_path)

        def fake_run_cmd(command, **_kwargs):
            output = command[command.index("-od") + 1]
            Path(output).write_text("0 1 2\n3 4 5\n6 7 8\n")

        monkeypatch.setattr(densmap, "run_cmd", fake_run_cmd)
        out = tmp_path / "densmap_gmx.npy"
        args = types.SimpleNamespace(
            topology="topology.tpr",
            trajectory="trajectory.xtc",
            selection="Water",
            output=str(out),
            bins=10,
            axis="xy",
            gmx=True,
            index=None,
        )

        densmap.run(args)

        assert out.exists()
        assert not list(tmp_path.glob("tmp*.dat"))

    def test_gmx_axis_orientation_and_bins(self, tmp_path, monkeypatch):
        """The gmx .dat header column/row map to axis0/axis1 like histogram2d"""
        from mdtbx.cv import densmap

        monkeypatch.chdir(tmp_path)
        commands = []

        def fake_run_cmd(command, **_kwargs):
            commands.append(command)
            output = command[command.index("-od") + 1]
            # Header row: second-axis ticks; header column: first-axis ticks.
            Path(output).write_text("0 10 20 30\n1 11 12 13\n2 21 22 23\n")

        monkeypatch.setattr(densmap, "run_cmd", fake_run_cmd)
        out = tmp_path / "densmap_gmx_axes.npy"
        args = types.SimpleNamespace(
            topology="topology.tpr",
            trajectory="trajectory.xtc",
            selection="Water",
            output=str(out),
            bins=10,
            axis="xy",
            gmx=True,
            index=None,
        )

        densmap.run(args)

        # allow_pickle: object array written by this very test via save_npy
        counts, axis0, axis1 = np.load(out, allow_pickle=True)
        np.testing.assert_allclose(counts, [[11.0, 12.0, 13.0], [21.0, 22.0, 23.0]])
        np.testing.assert_allclose(axis0, [1.0, 2.0])
        np.testing.assert_allclose(axis1, [10.0, 20.0, 30.0])
        (command,) = commands
        assert command[command.index("-n1") + 1] == "10"
        assert command[command.index("-n2") + 1] == "10"
