"""
analysis/extract_str のユニットテスト

MDtraj を使った特定フレームの構造抽出をテストする（gmx=False パス）。
"""

import types


class TestExtractStrRun:
    def _make_args(self, traj_files, output, time=1, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            time=time,
            output=str(output),
            gmx=False,
            index=None,
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """run() が PDB ファイルを生成すること"""
        from src.analysis.extract_str import run

        out = tmp_path / "frame1.pdb"
        run(self._make_args(trajectory_files, out, time=1))
        assert out.exists()

    def test_output_is_valid_pdb(self, trajectory_files, tmp_path):
        """生成された PDB が MDtraj で読み込めること"""
        import mdtraj as md

        from src.analysis.extract_str import run

        out = tmp_path / "frame1.pdb"
        run(self._make_args(trajectory_files, out, time=1))

        loaded = md.load_pdb(str(out))
        assert loaded.n_frames == 1

    def test_output_atom_count(self, trajectory_files, tmp_path):
        """抽出した構造の原子数が元の軌跡と一致すること"""
        import mdtraj as md

        from src.analysis.extract_str import run

        out = tmp_path / "frame1.pdb"
        run(self._make_args(trajectory_files, out, time=1))

        loaded = md.load_pdb(str(out))
        assert loaded.n_atoms == trajectory_files["traj"].n_atoms

    def test_extract_different_frames(self, trajectory_files, tmp_path):
        """異なるフレームで異なる座標が得られること"""
        import mdtraj as md

        from src.analysis.extract_str import run

        out1 = tmp_path / "frame1.pdb"
        out2 = tmp_path / "frame2.pdb"
        run(self._make_args(trajectory_files, out1, time=1))
        run(self._make_args(trajectory_files, out2, time=2))

        t1 = md.load_pdb(str(out1))
        t2 = md.load_pdb(str(out2))
        # ランダム座標なので 2 フレームの座標は異なるはず
        import numpy as np

        assert not np.allclose(t1.xyz, t2.xyz)
