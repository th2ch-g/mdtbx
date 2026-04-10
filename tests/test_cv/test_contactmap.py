"""
cv/contactmap のユニットテスト

MDtraj ベースの residue contact matrix を検証する。
"""

import types

import numpy as np


class TestContactmapRun:
    def _make_args(self, traj_files, output, trajectory=True):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"] if trajectory else None,
            selection="protein and name CA",
            cutoff=100.0,
            output=str(output),
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        from src.cv.contactmap import run

        out = tmp_path / "contactmap.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        from src.cv.contactmap import run

        out = tmp_path / "contactmap_shape.npy"
        run(self._make_args(trajectory_files, out))

        contactmap = np.load(out)
        assert contactmap.shape == (2, 2)

    def test_diagonal_is_zero(self, trajectory_files, tmp_path):
        from src.cv.contactmap import run

        out = tmp_path / "contactmap_diag.npy"
        run(self._make_args(trajectory_files, out))

        contactmap = np.load(out)
        assert np.allclose(np.diag(contactmap), 0.0)

    def test_single_structure_path_supported(self, trajectory_files, tmp_path):
        from src.cv.contactmap import run

        out = tmp_path / "contactmap_pdb.npy"
        run(self._make_args(trajectory_files, out, trajectory=False))

        contactmap = np.load(out)
        assert contactmap.shape == (2, 2)
