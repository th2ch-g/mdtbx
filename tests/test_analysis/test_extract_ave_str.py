"""
analysis/extract_ave_str のユニットテスト

MDtraj を使った平均構造抽出をテストする（gmx=False パス）。
"""

import types

import numpy as np


class TestExtractAveStrRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            output=str(output),
            gmx=False,
            index=None,
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """run() が PDB ファイルを生成すること"""
        from src.analysis.extract_ave_str import run

        out = tmp_path / "ave.pdb"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_has_single_frame(self, trajectory_files, tmp_path):
        """平均構造は 1 フレームのみであること"""
        import mdtraj as md

        from src.analysis.extract_ave_str import run

        out = tmp_path / "ave_single.pdb"
        run(self._make_args(trajectory_files, out))

        loaded = md.load_pdb(str(out))
        assert loaded.n_frames == 1

    def test_output_has_correct_atom_count_all(self, trajectory_files, tmp_path):
        """全原子選択時の原子数が元の軌跡と一致すること"""
        import mdtraj as md

        from src.analysis.extract_ave_str import run

        out = tmp_path / "ave_all.pdb"
        run(self._make_args(trajectory_files, out, selection="all"))

        loaded = md.load_pdb(str(out))
        assert loaded.n_atoms == trajectory_files["traj"].n_atoms

    def test_output_has_correct_atom_count_subset(self, trajectory_files, tmp_path):
        """部分選択時の原子数が選択した原子数と一致すること"""
        import mdtraj as md

        from src.analysis.extract_ave_str import run

        selection = "resid 0"
        out = tmp_path / "ave_subset.pdb"
        run(self._make_args(trajectory_files, out, selection=selection))

        loaded = md.load_pdb(str(out))
        n_selected = len(trajectory_files["traj"].topology.select(selection))
        assert loaded.n_atoms == n_selected

    def test_average_xyz_is_within_trajectory_range(self, trajectory_files, tmp_path):
        """平均座標が元の軌跡の各軸の min/max の範囲内であること"""
        import mdtraj as md

        from src.analysis.extract_ave_str import run

        out = tmp_path / "ave_range.pdb"
        run(self._make_args(trajectory_files, out))

        ave = md.load_pdb(str(out))
        traj = trajectory_files["traj"]

        assert np.all(ave.xyz >= traj.xyz.min(axis=0) - 1e-6)
        assert np.all(ave.xyz <= traj.xyz.max(axis=0) + 1e-6)
