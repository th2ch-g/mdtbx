"""
trajectory/fit のユニットテスト

MDtraj を使った軌跡フィッティングをテストする（gmx=False パス）。
"""

import types


class TestFitRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            file=traj_files["xtc"],
            topology=traj_files["pdb"],
            output=str(output),
            selection=selection,
            gmx=False,
            pbc="mol",
            index="index.ndx",
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """run() がフィット済み軌跡ファイルを生成すること"""

        from src.trajectory.fit import run

        out = tmp_path / "fitted.xtc"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_has_same_frames(self, trajectory_files, tmp_path):
        """フィット後の軌跡のフレーム数が元と同じであること"""
        import mdtraj as md

        from src.trajectory.fit import run

        out = tmp_path / "fitted.xtc"
        run(self._make_args(trajectory_files, out))

        fitted = md.load(str(out), top=trajectory_files["pdb"])
        original_n_frames = trajectory_files["traj"].n_frames
        assert fitted.n_frames == original_n_frames

    def test_output_has_same_atoms(self, trajectory_files, tmp_path):
        """フィット後の軌跡の原子数が元と同じであること"""
        import mdtraj as md

        from src.trajectory.fit import run

        out = tmp_path / "fitted.xtc"
        run(self._make_args(trajectory_files, out))

        fitted = md.load(str(out), top=trajectory_files["pdb"])
        original_n_atoms = trajectory_files["traj"].n_atoms
        assert fitted.n_atoms == original_n_atoms
