"""
cv/comdist のユニットテスト

合成軌跡を使って重心間距離（COM distance）の計算を検証する。
"""

import types

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
        """.npy ファイルが生成されること"""
        from src.cv.comdist import run

        out = tmp_path / "comdist.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        """出力配列の長さがフレーム数と一致すること"""
        from src.cv.comdist import run

        out = tmp_path / "comdist.npy"
        run(self._make_args(trajectory_files, out))

        dist = np.load(str(out))
        n_frames = trajectory_files["traj"].n_frames
        assert dist.shape == (n_frames,)

    def test_output_nonnegative(self, trajectory_files, tmp_path):
        """距離は常に非負であること"""
        from src.cv.comdist import run

        out = tmp_path / "comdist.npy"
        run(self._make_args(trajectory_files, out))

        dist = np.load(str(out))
        assert np.all(dist >= 0)

    def test_same_selection_gives_zero(self, trajectory_files, tmp_path):
        """同じ原子群の重心間距離は 0 になること"""
        from src.cv.comdist import run

        out = tmp_path / "comdist_zero.npy"
        args = self._make_args(trajectory_files, out, sel1="resid 0", sel2="resid 0")
        run(args)

        dist = np.load(str(out))
        assert np.allclose(dist, 0.0, atol=1e-6)
