import types
from pathlib import Path

import numpy as np
import parmed as pmd

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

    @property
    def coordinates(self):
        return self._coordinates.copy()

    @coordinates.setter
    def coordinates(self, value):
        self._coordinates = np.asarray(value, dtype=float).copy()
        self.coordinate_set_count += 1

    def __getitem__(self, indices):
        return _FakeStructure(self._coordinates[list(indices)])

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

    packmol_water.repack_water(
        tmp_path / "system.parm7",
        tmp_path / "system.rst7",
        tmp_path / "system.pdb",
    )

    assert original.coordinate_set_count == 1
    np.testing.assert_allclose(original._coordinates[:2], original_coordinates[:2])
    np.testing.assert_allclose(original._coordinates[2:], packed_coordinates[2:])
