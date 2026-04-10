"""
cv/rmsf のユニットテスト

MDtraj を使った RMSF 計算をテストする（gmx=False パス）。
"""

import types

import numpy as np


class TestRmsfRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            output=str(output),
            gmx=False,
            resolution="atom",
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy ファイルが生成されること"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape_all_atoms(self, trajectory_files, tmp_path):
        """全原子選択時の出力長が原子数と一致すること"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_all.npy"
        run(self._make_args(trajectory_files, out))
        rmsf = np.load(str(out))

        n_atoms = trajectory_files["traj"].n_atoms
        assert rmsf.shape == (n_atoms,)

    def test_output_shape_subset(self, trajectory_files, tmp_path):
        """部分選択時の出力長が選択原子数と一致すること"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_subset.npy"
        run(self._make_args(trajectory_files, out, selection="resid 0"))
        rmsf = np.load(str(out))

        n_selected = len(trajectory_files["traj"].topology.select("resid 0"))
        assert rmsf.shape == (n_selected,)

    def test_output_nonnegative(self, trajectory_files, tmp_path):
        """RMSF は常に非負であること"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_nn.npy"
        run(self._make_args(trajectory_files, out))
        rmsf = np.load(str(out))
        assert np.all(rmsf >= 0)

    def test_static_trajectory_gives_zero_rmsf(self, tmp_path_factory):
        """全フレームが同一座標の軌跡では RMSF が 0 になること"""
        import mdtraj as md

        from src.cv.rmsf import run

        tmp = tmp_path_factory.mktemp("static")

        # 静止軌跡を作成
        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("ALA", chain)
        top.add_atom("CA", md.element.carbon, res)

        xyz = np.zeros((5, 1, 3))
        traj = md.Trajectory(xyz, top)

        pdb_path = str(tmp / "static.pdb")
        xtc_path = str(tmp / "static.xtc")
        traj[0].save_pdb(pdb_path)
        traj.save_xtc(xtc_path)

        out = tmp / "rmsf_static.npy"
        run(
            types.SimpleNamespace(
                topology=pdb_path,
                trajectory=xtc_path,
                selection="all",
                output=str(out),
                gmx=False,
                resolution="atom",
            )
        )
        rmsf = np.load(str(out))
        assert np.allclose(rmsf, 0.0, atol=1e-6)
