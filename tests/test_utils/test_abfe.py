import json
import math

import pytest

from src.utils.abfe import (
    ABFE_MANIFEST,
    boresch_pull_settings,
    boresch_standard_state_correction,
    calculate_anchor_geometry,
    load_abfe_manifest,
    write_anchor_index,
)


def _anchor_pdb(path):
    coordinates = [
        (0.0, 10.0, 0.0),
        (0.0, 0.0, 0.0),
        (10.0, 0.0, 0.0),
        (10.0, 10.0, 0.0),
        (20.0, 10.0, 0.0),
        (20.0, 10.0, 10.0),
    ]
    lines = []
    for index, (x, y, z) in enumerate(coordinates, start=1):
        residue = "REC" if index <= 3 else "LIG"
        lines.append(
            f"ATOM  {index:5d}  C{index:<2} {residue:>3} A{index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
    path.write_text("".join(lines) + "END\n")
    return path


def test_anchor_geometry_and_standard_state_correction(tmp_path):
    structure = _anchor_pdb(tmp_path / "anchors.pdb")

    geometry = calculate_anchor_geometry(structure, [1, 2, 3, 4, 5, 6])
    correction = boresch_standard_state_correction(
        geometry,
        temperature=300.0,
        distance_spring=4184.0,
        angle_spring=41.84,
        dihedral_spring=41.84,
    )

    assert geometry["distance_nm"] == pytest.approx(1.0)
    assert geometry["angles_rad"] == pytest.approx([math.pi / 2, math.pi / 2])
    assert math.isfinite(correction)


def test_pull_settings_release_all_six_restraints():
    geometry = {
        "distance_nm": 0.4,
        "angles_rad": [math.pi / 2, math.pi / 2],
        "dihedrals_rad": [0.1, 0.2, 0.3],
    }

    settings = boresch_pull_settings(
        geometry,
        distance_spring=4184.0,
        angle_spring=41.84,
        dihedral_spring=41.84,
        release=True,
    )

    assert settings["pull-ncoords"] == 6
    assert settings["pull-coord1-geometry"] == "distance"
    assert settings["pull-coord6-geometry"] == "dihedral"
    assert all(settings[f"pull-coord{index}-kB"] == "0" for index in range(1, 7))


def test_write_anchor_index_preserves_existing_groups(tmp_path):
    source = tmp_path / "base.ndx"
    source.write_text("[ System ]\n1 2 3 4 5 6\n")
    output = tmp_path / "abfe.ndx"

    write_anchor_index(output, [1, 2, 3, 4, 5, 6], source)

    text = output.read_text()
    assert "[ System ]" in text
    assert "[ ABFE_6 ]\n6" in text


def test_load_abfe_manifest(tmp_path):
    legs = {
        "solvent_charge": "solvent_charge",
        "solvent_vdw": "solvent_vdw",
        "complex_charge": "complex_charge",
        "complex_vdw": "complex_vdw",
        "complex_restraint": "complex_restraint",
    }
    manifest = {"schema_version": 1, "legs": legs}
    (tmp_path / ABFE_MANIFEST).write_text(json.dumps(manifest))

    assert load_abfe_manifest(tmp_path) == (tmp_path, manifest)
