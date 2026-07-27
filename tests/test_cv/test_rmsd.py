"""
cv/rmsd unit tests

Verifies the accuracy of the RMSD computation on synthetic trajectory data.
"""

import types

import numpy as np
import pytest


class TestRmsdRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            reference=traj_files["pdb"],
            selection_fit_trj=selection,
            selection_fit_ref=selection,
            selection_cal_trj=selection,
            selection_cal_ref=selection,
            output=str(output),
        )

    def test_rmsd_output_file_created(self, trajectory_files, tmp_path):
        """run() generates a .npy file"""
        from mdtbx.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)
        assert out.exists()

    def test_rmsd_shape(self, trajectory_files, tmp_path):
        """The RMSD array length matches the trajectory frame count"""
        from mdtbx.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        n_frames = trajectory_files["traj"].n_frames
        assert rmsd.shape == (n_frames,)

    def test_rmsd_nonnegative(self, trajectory_files, tmp_path):
        """RMSD is always non-negative"""
        from mdtbx.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        assert np.all(rmsd >= 0)

    def test_rmsd_first_frame_near_zero(self, trajectory_files, tmp_path):
        """
        The reference is the PDB of frame 0, so the RMSD of frame 0 is
        close to 0 after fitting.
        """
        from mdtbx.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        # XTC is float32, so a round trip through PDB (float64) loses precision
        assert rmsd[0] == pytest.approx(0.0, abs=1e-3)

    def test_rmsd_other_frames_nonzero(self, trajectory_files, tmp_path):
        """Other frames hold random coordinates, so their RMSD exceeds 0"""
        from mdtbx.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        # Check that some frame other than frame 0 exceeds 0
        assert np.any(rmsd[1:] > 0)

    def test_empty_fit_selection_raises(self, trajectory_files, tmp_path):
        from mdtbx.cv.rmsd import run

        args = self._make_args(
            trajectory_files, tmp_path / "rmsd.npy", selection="name XXXX"
        )
        with pytest.raises(ValueError, match="No atoms selected"):
            run(args)
