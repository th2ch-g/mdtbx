"""
cv/xyz のユニットテスト

MDtraj を使った XYZ 座標抽出をテストする（gmx=False パス）。
"""

import types

import numpy as np


class TestXyzRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            output=str(output),
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy ファイルが生成されること"""
        from src.cv.xyz import run

        out = tmp_path / "xyz.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape_all_atoms(self, trajectory_files, tmp_path):
        """全原子選択時の形状が (n_frames, n_atoms, 3) であること"""
        from src.cv.xyz import run

        out = tmp_path / "xyz_all.npy"
        run(self._make_args(trajectory_files, out))
        xyz = np.load(str(out))

        n_frames = trajectory_files["traj"].n_frames
        n_atoms = trajectory_files["traj"].n_atoms
        assert xyz.shape == (n_frames, n_atoms, 3)

    def test_output_shape_subset(self, trajectory_files, tmp_path):
        """部分選択時の形状が (n_frames, n_selected, 3) であること"""
        from src.cv.xyz import run

        out = tmp_path / "xyz_subset.npy"
        run(self._make_args(trajectory_files, out, selection="resid 0"))
        xyz = np.load(str(out))

        n_frames = trajectory_files["traj"].n_frames
        n_selected = len(trajectory_files["traj"].topology.select("resid 0"))
        assert xyz.shape == (n_frames, n_selected, 3)

    def test_coordinates_are_finite(self, trajectory_files, tmp_path):
        """出力座標に NaN / Inf が含まれないこと"""
        from src.cv.xyz import run

        out = tmp_path / "xyz_finite.npy"
        run(self._make_args(trajectory_files, out))
        xyz = np.load(str(out))
        assert np.all(np.isfinite(xyz))
