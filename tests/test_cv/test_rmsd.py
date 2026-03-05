"""
cv/rmsd のユニットテスト

合成軌跡データを使って RMSD 計算の正確性を検証する。
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
        """run() が .npy ファイルを生成すること"""
        from src.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)
        assert out.exists()

    def test_rmsd_shape(self, trajectory_files, tmp_path):
        """RMSD 配列の長さが軌跡のフレーム数と一致すること"""
        from src.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        n_frames = trajectory_files["traj"].n_frames
        assert rmsd.shape == (n_frames,)

    def test_rmsd_nonnegative(self, trajectory_files, tmp_path):
        """RMSD は常に非負であること"""
        from src.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        assert np.all(rmsd >= 0)

    def test_rmsd_first_frame_near_zero(self, trajectory_files, tmp_path):
        """
        reference は frame 0 の PDB であるため、
        frame 0 の RMSD は fitting 後にほぼ 0 になること。
        """
        from src.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        # XTC は float32 のため PDB (float64) との往復で精度損失がある
        assert rmsd[0] == pytest.approx(0.0, abs=1e-3)

    def test_rmsd_other_frames_nonzero(self, trajectory_files, tmp_path):
        """ランダムな座標を持つ他のフレームの RMSD は 0 より大きいこと"""
        from src.cv.rmsd import run

        out = tmp_path / "rmsd.npy"
        args = self._make_args(trajectory_files, out)
        run(args)

        rmsd = np.load(str(out))
        # frame 0 以外のどこかが 0 より大きいことを確認
        assert np.any(rmsd[1:] > 0)
