"""
cv/distmap のユニットテスト

MDtraj ベースの residue distance matrix を検証する。
"""

import types

import numpy as np
import pytest


class TestDistmapRun:
    def _make_args(self, traj_files, output, trajectory=True):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"] if trajectory else None,
            selection="protein and name CA",
            output=str(output),
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        from src.cv.distmap import run

        out = tmp_path / "distmap.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        from src.cv.distmap import run

        out = tmp_path / "distmap_shape.npy"
        run(self._make_args(trajectory_files, out))

        distmap = np.load(out)
        assert distmap.shape == (2, 2)

    def test_output_is_symmetric(self, trajectory_files, tmp_path):
        from src.cv.distmap import run

        out = tmp_path / "distmap_sym.npy"
        run(self._make_args(trajectory_files, out))

        distmap = np.load(out)
        assert np.allclose(distmap, distmap.T)

    def test_single_structure_path_supported(self, trajectory_files, tmp_path):
        from src.cv.distmap import run

        out = tmp_path / "distmap_pdb.npy"
        run(self._make_args(trajectory_files, out, trajectory=False))

        distmap = np.load(out)
        assert distmap.shape == (2, 2)


@pytest.mark.parametrize(
    "coordinates",
    [
        np.empty((0, 2, 3)),
        np.empty((1, 0, 3)),
        np.empty((1, 2, 2)),
        np.array([[[np.nan, 0.0, 0.0]]]),
    ],
)
def test_pairwise_distances_rejects_invalid_coordinates(coordinates):
    from src.cv.distmap import pairwise_distances

    with pytest.raises(ValueError, match="coordinates"):
        pairwise_distances(coordinates)
