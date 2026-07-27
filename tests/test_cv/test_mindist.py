"""
cv/mindist unit tests

Verifies the minimum interatomic distance on a synthetic trajectory.
"""

import types

import numpy as np


class TestMindistRun:
    def _make_args(self, traj_files, output, sel1="resid 0", sel2="resid 1"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection1=sel1,
            selection2=sel2,
            output=str(output),
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy file is generated"""
        from src.cv.mindist import run

        out = tmp_path / "mindist.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        """The output length matches the frame count"""
        from src.cv.mindist import run

        out = tmp_path / "mindist_shape.npy"
        run(self._make_args(trajectory_files, out))
        dist = np.load(str(out))

        n_frames = trajectory_files["traj"].n_frames
        assert dist.shape == (n_frames,)

    def test_output_nonnegative(self, trajectory_files, tmp_path):
        """Distances are always non-negative"""
        from src.cv.mindist import run

        out = tmp_path / "mindist_nn.npy"
        run(self._make_args(trajectory_files, out))
        dist = np.load(str(out))
        assert np.all(dist >= 0)

    def test_mindist_consistent_across_runs(self, trajectory_files, tmp_path):
        """Running twice on the same input gives the same result (deterministic)"""
        from src.cv.mindist import run

        out1 = tmp_path / "mindist_run1.npy"
        out2 = tmp_path / "mindist_run2.npy"
        run(self._make_args(trajectory_files, out1))
        run(self._make_args(trajectory_files, out2))

        dist1 = np.load(str(out1))
        dist2 = np.load(str(out2))
        assert np.allclose(dist1, dist2)
