import types
from pathlib import Path

import numpy as np
import parmed as pmd
import pytest

from mdtbx.build import packmol_water


def test_build_packmol_input_uses_box_bounds_and_fixed_solute(tmp_path):
    text = packmol_water._build_packmol_input(
        fixed_path=tmp_path / "fixed.pdb",
        water_path=tmp_path / "water.pdb",
        output_path=tmp_path / "packed.pdb",
        water_count=790,
        box=[32.0, 33.0, 34.0, 90.0, 90.0, 90.0],
        seed=42,
        tolerance=2.0,
    )

    assert "tolerance 2.0" in text
    assert "seed 42" in text
    assert "inside box 1.0 1.0 1.0 31.0 32.0 33.0" in text
    assert "number 790" in text
    assert "fixed 0.0 0.0 0.0 0.0 0.0 0.0" in text


def test_build_packmol_input_supports_pure_water(tmp_path):
    text = packmol_water._build_packmol_input(
        fixed_path=None,
        water_path=tmp_path / "water.pdb",
        output_path=tmp_path / "packed.pdb",
        water_count=100,
        box=[20.0, 20.0, 20.0, 90.0, 90.0, 90.0],
        seed=1,
        tolerance=2.0,
    )

    assert "fixed" not in text
    assert text.count("structure ") == 1


def test_water_atom_groups_uses_residue_names():
    residues = [
        types.SimpleNamespace(
            name="ALA",
            atoms=[types.SimpleNamespace(idx=0)],
        ),
        types.SimpleNamespace(
            name="WAT",
            atoms=[
                types.SimpleNamespace(idx=1),
                types.SimpleNamespace(idx=2),
                types.SimpleNamespace(idx=3),
            ],
        ),
    ]
    structure = types.SimpleNamespace(residues=residues)

    assert packmol_water._water_atom_groups(structure, {"WAT"}) == [[1, 2, 3]]


class _FakeStructure:
    def __init__(self, coordinates, residues=None):
        self._coordinates = np.asarray(coordinates, dtype=float).copy()
        self.atoms = [types.SimpleNamespace(idx=i) for i in range(len(coordinates))]
        self.residues = residues or []
        self.box = [20.0, 20.0, 20.0, 90.0, 90.0, 90.0]
        self.coordinate_set_count = 0
        self.bonds = []

    @property
    def coordinates(self):
        return self._coordinates.copy()

    @coordinates.setter
    def coordinates(self, value):
        self._coordinates = np.asarray(value, dtype=float).copy()
        self.coordinate_set_count += 1

    def __getitem__(self, indices):
        selection = list(indices)
        if len(selection) == len(self.atoms):
            selection = [index for index, keep in enumerate(selection) if keep]
        return _FakeStructure(self._coordinates[selection])

    def save(self, path, **_kwargs):
        Path(path).write_text("fake structure\n")


def test_repack_water_writes_packed_coordinates_back_once(tmp_path, monkeypatch):
    residues = [
        types.SimpleNamespace(name="ALA", atoms=[]),
        types.SimpleNamespace(
            name="WAT",
            atoms=[types.SimpleNamespace(idx=i) for i in (2, 3, 4)],
        ),
        types.SimpleNamespace(
            name="WAT",
            atoms=[types.SimpleNamespace(idx=i) for i in (5, 6, 7)],
        ),
    ]
    original_coordinates = np.arange(24, dtype=float).reshape(8, 3)
    packed_coordinates = original_coordinates.copy()
    packed_coordinates[2:] += 100.0
    original = _FakeStructure(original_coordinates, residues)
    packed = _FakeStructure(packed_coordinates)

    def fake_load_file(path, **_kwargs):
        return packed if Path(path).name == "packed.pdb" else original

    def fake_run_cmd(_command, **kwargs):
        (Path(kwargs["cwd"]) / "packed.pdb").write_text("fake packed structure\n")
        return types.SimpleNamespace(returncode=0, stdout="Success!", stderr="")

    monkeypatch.setattr(pmd, "load_file", fake_load_file)
    monkeypatch.setattr(packmol_water, "run_cmd", fake_run_cmd)

    report = packmol_water.repack_water(
        tmp_path / "system.parm7",
        tmp_path / "system.rst7",
        tmp_path / "system.pdb",
    )

    assert original.coordinate_set_count == 1
    np.testing.assert_allclose(original._coordinates[:2], original_coordinates[:2])
    np.testing.assert_allclose(original._coordinates[2:], packed_coordinates[2:])
    assert report["transfer_max_abs_error_A"] == 0.0
    assert report["saved_max_abs_error_A"] == 0.0
    assert report["vectorized_transfer"] is True


def test_atom_selection_mask_keeps_atom_zero():
    assert packmol_water._atom_selection_mask(3, [0, 1, 2]) == [True, True, True]


def test_wrap_fixed_coordinates_wraps_single_atom_ions():
    structure = _FakeStructure([[-1.0, 21.0, 2.0], [5.0, 5.0, 5.0]])

    coordinates, moved = packmol_water._wrap_fixed_coordinates(
        structure, [0, 1], np.asarray(structure.box[:3])
    )

    np.testing.assert_allclose(coordinates[0], [19.0, 1.0, 2.0])
    np.testing.assert_allclose(coordinates[1], [5.0, 5.0, 5.0])
    assert moved == 1


def test_wrap_fixed_coordinates_rejects_split_bond():
    structure = _FakeStructure([[-1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    structure.bonds = [
        types.SimpleNamespace(atom1=structure.atoms[0], atom2=structure.atoms[1])
    ]

    with pytest.raises(ValueError, match="split a bonded component"):
        packmol_water._wrap_fixed_coordinates(
            structure, [0, 1], np.asarray(structure.box[:3])
        )


def test_align_coordinates_to_reference_removes_uniform_box_shift():
    reference = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    restart = reference + np.asarray([10.0, 11.0, 12.0])

    aligned, translation, max_error = packmol_water._align_coordinates_to_reference(
        restart, reference
    )

    np.testing.assert_allclose(aligned, reference)
    np.testing.assert_allclose(translation, [10.0, 11.0, 12.0])
    assert max_error == 0.0
