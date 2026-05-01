import argparse
from unittest.mock import MagicMock

import numpy as np
import pytest

from src.build import place


def _parse_args(add_subcmd, argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    add_subcmd(subparsers)
    return parser.parse_args(argv)


def test_random_rotation_matrix_is_orthogonal():
    R = place.random_rotation_matrix(seed=42)
    assert R.shape == (3, 3)
    # R @ R.T == I
    np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-10)
    # det(R) == +1 (proper rotation, not reflection)
    np.testing.assert_allclose(np.linalg.det(R), 1.0, atol=1e-10)


def test_random_rotation_matrix_is_deterministic():
    R1 = place.random_rotation_matrix(seed=42)
    R2 = place.random_rotation_matrix(seed=42)
    np.testing.assert_array_equal(R1, R2)


def test_random_rotation_matrix_seed_changes_result():
    R1 = place.random_rotation_matrix(seed=42)
    R2 = place.random_rotation_matrix(seed=43)
    assert not np.allclose(R1, R2)


def test_random_unit_vector_is_unit_norm():
    v = place.random_unit_vector(seed=42)
    assert v.shape == (3,)
    np.testing.assert_allclose(np.linalg.norm(v), 1.0, atol=1e-12)


def test_random_unit_vector_is_deterministic():
    v1 = place.random_unit_vector(seed=42)
    v2 = place.random_unit_vector(seed=42)
    np.testing.assert_array_equal(v1, v2)


def test_random_unit_vector_seed_changes_result():
    v1 = place.random_unit_vector(seed=42)
    v2 = place.random_unit_vector(seed=43)
    assert not np.allclose(v1, v2)


def test_rotation_and_direction_are_independent_streams():
    # Same seed but different consumers must not share their underlying samples.
    R = place.random_rotation_matrix(seed=42)
    v = place.random_unit_vector(seed=42)
    # rotation acts on direction differently than identity would
    assert not np.allclose(R @ v, v)


def test_default_seed_is_42():
    args = _parse_args(
        place.add_subcmd,
        ["place", "-1", "a.pdb", "-2", "b.pdb", "-d", "30.0"],
    )
    assert args.seed == 42
    assert args.output == "placed.pdb"


def test_required_args_missing():
    with pytest.raises(SystemExit):
        _parse_args(place.add_subcmd, ["place", "-1", "a.pdb"])


def test_min_interchain_distance_basic():
    coords1 = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    coords2 = np.array([[3.0, 0.0, 0.0], [11.0, 0.0, 0.0]])
    # min pair: (10,0,0) and (11,0,0) -> 1.0
    assert place.min_interchain_distance(coords1, coords2) == pytest.approx(1.0)


def _install_fake_pymol(monkeypatch, coords1, coords2, *, distance=30.0):
    """Install a MagicMock pymol module with deterministic COM/coords."""
    fake_cmd = MagicMock()
    fake_cmd.centerofmass.side_effect = [
        [1.0, 2.0, 3.0],  # mol1 input
        [4.0, 5.0, 6.0],  # mol2 input
        [0.0, 0.0, 0.0],  # mol1 placed
        [distance, 0.0, 0.0],  # mol2 placed
    ]
    fake_cmd.get_coords.side_effect = [coords1, coords2]
    fake_pymol = MagicMock()
    fake_pymol.cmd = fake_cmd
    monkeypatch.setitem(__import__("sys").modules, "pymol", fake_pymol)
    return fake_cmd


def test_run_invokes_pymol_pipeline(monkeypatch, tmp_path):
    # well-separated atoms -> no clash
    coords1 = np.array([[0.0, 0.0, 0.0]])
    coords2 = np.array([[30.0, 0.0, 0.0]])
    fake_cmd = _install_fake_pymol(monkeypatch, coords1, coords2)

    args = _parse_args(
        place.add_subcmd,
        [
            "place",
            "-1",
            "a.pdb",
            "-2",
            "b.pdb",
            "-d",
            "30.0",
            "--seed",
            "7",
            "-o",
            str(tmp_path / "out.pdb"),
        ],
    )
    place.run(args)

    fake_cmd.load.assert_any_call("a.pdb", "mol1")
    fake_cmd.load.assert_any_call("b.pdb", "mol2")
    fake_cmd.transform_selection.assert_called_once()
    sel, ttt = fake_cmd.transform_selection.call_args[0]
    assert sel == "mol2"
    assert len(ttt) == 16
    assert ttt[12:16] == [0.0, 0.0, 0.0, 1.0]

    # Final translate call uses the seeded random unit vector * distance.
    expected_dir = place.random_unit_vector(seed=7)
    expected_translation = (expected_dir * 30.0).tolist()
    final_translate_args, final_translate_kwargs = fake_cmd.translate.call_args
    np.testing.assert_allclose(final_translate_args[0], expected_translation)
    assert final_translate_args[1] == "mol2"
    assert final_translate_kwargs == {"camera": 0}
    # Magnitude must equal d
    assert np.linalg.norm(final_translate_args[0]) == pytest.approx(30.0)

    fake_cmd.save.assert_called_once_with(str(tmp_path / "out.pdb"), "placed")


def test_run_aborts_on_clash_by_default(monkeypatch, tmp_path):
    # overlapping atoms -> serious clash (min distance = 0.5 A)
    coords1 = np.array([[0.0, 0.0, 0.0]])
    coords2 = np.array([[0.5, 0.0, 0.0]])
    fake_cmd = _install_fake_pymol(monkeypatch, coords1, coords2)

    args = _parse_args(
        place.add_subcmd,
        [
            "place",
            "-1",
            "a.pdb",
            "-2",
            "b.pdb",
            "-d",
            "30.0",
            "-o",
            str(tmp_path / "out.pdb"),
        ],
    )
    with pytest.raises(SystemExit) as exc_info:
        place.run(args)
    assert exc_info.value.code == 1
    # PDB must NOT be written
    fake_cmd.save.assert_not_called()
    fake_cmd.create.assert_not_called()


def test_run_warns_on_clash_with_ignore(monkeypatch, tmp_path, caplog):
    coords1 = np.array([[0.0, 0.0, 0.0]])
    coords2 = np.array([[0.5, 0.0, 0.0]])
    fake_cmd = _install_fake_pymol(monkeypatch, coords1, coords2)

    args = _parse_args(
        place.add_subcmd,
        [
            "place",
            "-1",
            "a.pdb",
            "-2",
            "b.pdb",
            "-d",
            "30.0",
            "--ignore",
            "-o",
            str(tmp_path / "out.pdb"),
        ],
    )
    place.run(args)
    # PDB must be written despite the clash
    fake_cmd.save.assert_called_once_with(str(tmp_path / "out.pdb"), "placed")


def test_clash_cutoff_is_configurable(monkeypatch, tmp_path):
    # min distance = 2.5 A; default cutoff (2.0) would pass, but 3.0 should abort
    coords1 = np.array([[0.0, 0.0, 0.0]])
    coords2 = np.array([[2.5, 0.0, 0.0]])
    _install_fake_pymol(monkeypatch, coords1, coords2)

    args = _parse_args(
        place.add_subcmd,
        [
            "place",
            "-1",
            "a.pdb",
            "-2",
            "b.pdb",
            "-d",
            "30.0",
            "--clash_cutoff",
            "3.0",
            "-o",
            str(tmp_path / "out.pdb"),
        ],
    )
    with pytest.raises(SystemExit):
        place.run(args)
