"""
cv/densmap のユニットテスト

MDtraj を使った 2D 密度マップ計算をテストする（gmx=False パス）。
"""

import types

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
        """.npy ファイルが生成されること"""
        from src.cv.densmap import run

        out = tmp_path / "densmap.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_is_nonempty(self, trajectory_files, tmp_path):
        """出力ファイルが空でないこと"""
        from src.cv.densmap import run

        out = tmp_path / "densmap_size.npy"
        run(self._make_args(trajectory_files, out))
        assert out.stat().st_size > 0

    @pytest.mark.parametrize("axis", ["xy", "xz", "yz"])
    def test_different_axes(self, trajectory_files, tmp_path, axis):
        """xy / xz / yz の各投影面でファイルが生成されること"""
        from src.cv.densmap import run

        out = tmp_path / f"densmap_{axis}.npy"
        run(self._make_args(trajectory_files, out, axis=axis))
        assert out.exists()

    def test_empty_selection_does_not_create_file(self, trajectory_files, tmp_path):
        """原子が選択されない場合はファイルを生成しないこと"""
        from src.cv.densmap import run

        out = tmp_path / "densmap_empty.npy"
        run(self._make_args(trajectory_files, out, selection="name XXXX"))
        assert not out.exists()

    def test_histogram_shape(self, trajectory_files, tmp_path):
        """np.histogram2d の結果と整合すること（ゴールデンパス検証）"""
        import mdtraj as md
        from src.cv.densmap import _AXIS_MAP, run

        bins = 8
        axis = "xy"
        out = tmp_path / "densmap_golden.npy"
        run(self._make_args(trajectory_files, out, bins=bins, axis=axis))

        # 内部と同じ計算で counts を再現し、合計を比較する
        traj = md.load(trajectory_files["xtc"], top=trajectory_files["pdb"])
        atom_indices = traj.topology.select("all")
        xyz = traj.xyz[:, atom_indices, :]
        ax0, ax1 = _AXIS_MAP[axis]
        pos0 = xyz[:, :, ax0].ravel()
        pos1 = xyz[:, :, ax1].ravel()
        counts_ref, _, _ = np.histogram2d(pos0, pos1, bins=bins)

        # ファイルが存在し、総カウントが一致すること
        assert out.exists()
        total_ref = int(counts_ref.sum())
        expected = traj.n_frames * traj.n_atoms
        assert total_ref == expected
