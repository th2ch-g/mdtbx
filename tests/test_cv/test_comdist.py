"""
cv/comdist unit tests

Verifies the centre-of-mass distance (COM distance) on a synthetic trajectory.
"""

import types
from pathlib import Path

import numpy as np


class TestComdistRun:
    def _make_args(self, traj_files, output, sel1="resid 0", sel2="resid 1"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection1=sel1,
            selection2=sel2,
            output=str(output),
            gmx=False,
            index=None,
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy file is generated"""
        from mdtbx.cv.comdist import run

        out = tmp_path / "comdist.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        """The output length matches the frame count"""
        from mdtbx.cv.comdist import run

        out = tmp_path / "comdist.npy"
        run(self._make_args(trajectory_files, out))

        dist = np.load(str(out))
        n_frames = trajectory_files["traj"].n_frames
        assert dist.shape == (n_frames,)

    def test_output_nonnegative(self, trajectory_files, tmp_path):
        """Distances are always non-negative"""
        from mdtbx.cv.comdist import run

        out = tmp_path / "comdist.npy"
        run(self._make_args(trajectory_files, out))

        dist = np.load(str(out))
        assert np.all(dist >= 0)

    def test_same_selection_gives_zero(self, trajectory_files, tmp_path):
        """The COM distance between identical atom groups is 0"""
        from mdtbx.cv.comdist import run

        out = tmp_path / "comdist_zero.npy"
        args = self._make_args(trajectory_files, out, sel1="resid 0", sel2="resid 0")
        run(args)

        dist = np.load(str(out))
        assert np.allclose(dist, 0.0, atol=1e-6)

    def test_gmx_single_frame_output(self, tmp_path, monkeypatch):
        from mdtbx.cv import comdist

        monkeypatch.chdir(tmp_path)

        def fake_run_cmd(command, **_kwargs):
            output = command[command.index("-oxyz") + 1]
            Path(output).write_text("0.0 1.0 2.0 2.0\n")

        monkeypatch.setattr(comdist, "run_cmd", fake_run_cmd)
        out = tmp_path / "comdist_gmx.npy"
        args = types.SimpleNamespace(
            topology="topol.tpr",
            trajectory="traj.xtc",
            selection1="Group1",
            selection2="Group2",
            output=str(out),
            gmx=True,
            index=None,
        )

        comdist.run(args)

        assert np.load(out).tolist() == [3.0]
