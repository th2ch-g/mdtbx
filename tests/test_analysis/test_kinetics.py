from argparse import _SubParsersAction
from types import SimpleNamespace

import numpy as np
import pytest

from mdtbx.analysis import cluster, msm, tica
from mdtbx.cli import create_parser


def _features(tmp_path):
    rng = np.random.default_rng(7)
    first = np.cumsum(rng.normal(size=(100, 3)), axis=0)
    second = np.cumsum(rng.normal(size=(80, 3)), axis=0)
    first_path = tmp_path / "first.npy"
    second_path = tmp_path / "second.npy"
    np.save(first_path, first)
    np.save(second_path, second)
    return first_path, second_path


def test_kinetic_commands_are_registered():
    parser = create_parser()
    subparsers = next(
        action for action in parser._actions if isinstance(action, _SubParsersAction)
    )
    assert {"tica", "cluster", "msm"} <= set(subparsers.choices)


def test_tica_cluster_msm_pipeline_preserves_trajectory_boundaries(tmp_path):
    first, second = _features(tmp_path)
    tica_path = tmp_path / "tica.npz"
    result = tica.run(
        SimpleNamespace(
            input=[str(first), str(second)],
            lagtime=2,
            n_components=2,
            epsilon=1e-6,
            scaling="kinetic_map",
            output=str(tica_path),
        )
    )
    assert result["trajectories"] == 2
    with np.load(tica_path, allow_pickle=False) as archive:
        assert archive["kind"].item() == "mdtbx.tica"
        assert archive["scores"].shape == (180, 2)
        assert archive["lengths"].tolist() == [100, 80]
        assert archive["components"].shape == (3, 2)

    cluster_path = tmp_path / "clusters.npz"
    result = cluster.run(
        SimpleNamespace(
            input=[str(tica_path)],
            n_clusters=4,
            seed=11,
            max_iter=100,
            n_jobs=1,
            output=str(cluster_path),
        )
    )
    assert result["frames"] == 180
    with np.load(cluster_path, allow_pickle=False) as archive:
        assert archive["kind"].item() == "mdtbx.cluster"
        assert archive["dtrajs"].shape == (180,)
        assert archive["lengths"].tolist() == [100, 80]
        assert archive["centers"].shape == (4, 2)

    msm_path = tmp_path / "msm.npz"
    result = msm.run(
        SimpleNamespace(
            input=[str(cluster_path)],
            lagtime=2,
            count_mode="effective",
            nonreversible=False,
            output=str(msm_path),
        )
    )
    assert 1 <= result["states"] <= 4
    with np.load(msm_path, allow_pickle=False) as archive:
        transition = archive["transition_matrix"]
        assert archive["kind"].item() == "mdtbx.msm"
        assert np.allclose(transition.sum(axis=1), 1.0)
        assert archive["reversible"].item() is True
        assert archive["state_symbols"].ndim == 1


def test_tica_rejects_nonfinite_features(tmp_path):
    path = tmp_path / "bad.npy"
    np.save(path, np.array([[0.0, np.nan], [1.0, 2.0], [2.0, 3.0]]))
    with pytest.raises(ValueError, match="NaN"):
        tica.run(
            SimpleNamespace(
                input=[str(path)],
                lagtime=1,
                n_components=1,
                epsilon=1e-6,
                scaling="none",
                output=str(tmp_path / "out.npz"),
            )
        )


def test_msm_rejects_negative_state_labels(tmp_path):
    path = tmp_path / "bad.npy"
    np.save(path, np.array([0, 1, -1, 0], dtype=np.int64))
    with pytest.raises(ValueError, match="negative"):
        msm.run(
            SimpleNamespace(
                input=[str(path)],
                lagtime=1,
                count_mode="sliding",
                nonreversible=False,
                output=str(tmp_path / "out.npz"),
            )
        )
