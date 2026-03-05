"""
cv/comvec のユニットテスト

重心ベクトル（COM vector）の計算を合成軌跡で検証する。
"""

import types

import numpy as np


class TestComvecRun:
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
        from src.cv.comvec import run

        out = tmp_path / "comvec.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        """出力配列の形状が (n_frames, 3) であること"""
        from src.cv.comvec import run

        out = tmp_path / "comvec_shape.npy"
        run(self._make_args(trajectory_files, out))
        vec = np.load(str(out))

        n_frames = trajectory_files["traj"].n_frames
        assert vec.shape == (n_frames, 3)

    def test_same_selection_gives_zero_vector(self, trajectory_files, tmp_path):
        """同じ原子群のベクトルは零ベクトルになること"""
        from src.cv.comvec import run

        out = tmp_path / "comvec_zero.npy"
        args = self._make_args(trajectory_files, out, sel1="resid 0", sel2="resid 0")
        run(args)
        vec = np.load(str(out))
        assert np.allclose(vec, 0.0, atol=1e-6)

    def test_antisymmetry(self, trajectory_files, tmp_path):
        """sel1/sel2 を入れ替えると符号が反転すること"""
        from src.cv.comvec import run

        out1 = tmp_path / "comvec_ab.npy"
        out2 = tmp_path / "comvec_ba.npy"
        run(self._make_args(trajectory_files, out1, sel1="resid 0", sel2="resid 1"))
        run(self._make_args(trajectory_files, out2, sel1="resid 1", sel2="resid 0"))

        vec_ab = np.load(str(out1))
        vec_ba = np.load(str(out2))
        assert np.allclose(vec_ab, -vec_ba, atol=1e-6)
