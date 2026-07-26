import types

from src.build import packmol_water


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
