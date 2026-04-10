"""
cv/pca のユニットテスト

MDtraj + scikit-learn を使った PCA をテストする（gmx=False パス）。
"""

import types

import numpy as np


class TestPcaRun:
    def _make_args(self, traj_files, output, n_components=2):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            reference=traj_files["pdb"],
            selection_cal_trj="all",
            selection_cal_ref="all",
            selection_fit_trj="all",
            selection_fit_ref="all",
            output=str(output),
            gmx=False,
            n_components=n_components,
            index=None,
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy ファイルが生成されること"""
        from src.cv.pca import run

        out = tmp_path / "pca.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        """出力形状が (n_frames, n_components) であること"""
        from src.cv.pca import run

        n_components = 3
        out = tmp_path / "pca_shape.npy"
        run(self._make_args(trajectory_files, out, n_components=n_components))
        pc = np.load(str(out))

        n_frames = trajectory_files["traj"].n_frames
        assert pc.shape == (n_frames, n_components)

    def test_first_pc_has_max_variance(self, trajectory_files, tmp_path):
        """PC1 の分散が PC2 以上であること（PCA の定義）"""
        from src.cv.pca import run

        out = tmp_path / "pca_var.npy"
        run(self._make_args(trajectory_files, out, n_components=2))
        pc = np.load(str(out))

        var1 = np.var(pc[:, 0])
        var2 = np.var(pc[:, 1])
        assert var1 >= var2

    def test_subset_selection(self, trajectory_files, tmp_path):
        """部分選択でも正常に計算できること"""
        from src.cv.pca import run

        out = tmp_path / "pca_subset.npy"
        args = types.SimpleNamespace(
            topology=trajectory_files["pdb"],
            trajectory=trajectory_files["xtc"],
            reference=trajectory_files["pdb"],
            selection_cal_trj="resid 0",
            selection_cal_ref="resid 0",
            selection_fit_trj="resid 0",
            selection_fit_ref="resid 0",
            output=str(out),
            gmx=False,
            n_components=2,
            index=None,
        )
        run(args)
        pc = np.load(str(out))
        assert pc.shape[0] == trajectory_files["traj"].n_frames
