import json
import math

import pytest

from mdtbx.utils.abfe import (
    ABFE_MANIFEST,
    boresch_pull_settings,
    boresch_standard_state_correction,
    calculate_anchor_geometry,
    calculate_trajectory_anchor_geometry,
    load_abfe_manifest,
    select_trajectory_boresch_anchors,
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


def test_trajectory_anchor_geometry_uses_all_strided_frames(tmp_path):
    import mdtraj as md
    import numpy as np

    structure = _anchor_pdb(tmp_path / "anchors.pdb")
    frame = md.load(str(structure))
    xyz = np.repeat(frame.xyz, 3, axis=0)
    xyz[1, 3, 1] += 0.02
    xyz[2, 3, 1] -= 0.02
    trajectory = tmp_path / "anchors.xtc"
    md.Trajectory(xyz, frame.topology).save_xtc(str(trajectory))

    geometry, deviations, frame_count = calculate_trajectory_anchor_geometry(
        trajectory,
        structure,
        [1, 2, 3, 4, 5, 6],
        stride=2,
    )

    assert frame_count == 2
    assert geometry["distance_nm"] == pytest.approx(0.99, abs=1e-3)
    assert deviations["distance_nm"] >= 0.0


def test_trajectory_anchor_selection_returns_low_variance_candidate(tmp_path):
    import mdtraj as md
    import numpy as np

    atoms = [
        ("N", "REC", 1, (0.0, 1.0, 0.0)),
        ("CA", "REC", 1, (0.0, 0.0, 0.0)),
        ("C", "REC", 1, (1.0, 0.0, 0.0)),
        ("C1", "LIG", 2, (1.0, 1.0, 0.0)),
        ("C2", "LIG", 2, (2.0, 1.0, 0.0)),
        ("C3", "LIG", 2, (2.0, 1.0, 1.0)),
    ]
    structure = tmp_path / "selection.pdb"
    lines = []
    for index, (name, residue, residue_id, xyz) in enumerate(atoms, start=1):
        x, y, z = xyz
        lines.append(
            f"ATOM  {index:5d} {name:^4s} {residue:>3s} A{residue_id:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
    structure.write_text("".join(lines) + "END\n")
    frame = md.load(str(structure))
    trajectory = tmp_path / "selection.xtc"
    md.Trajectory(np.repeat(frame.xyz, 2, axis=0), frame.topology).save_xtc(
        str(trajectory)
    )

    anchors, _geometry, deviations, frame_count, score = (
        select_trajectory_boresch_anchors(
            trajectory,
            structure,
            receptor_selection="resname REC",
            ligand_selection="resname LIG",
            distance_spring=4184.0,
            angle_spring=41.84,
            dihedral_spring=41.84,
            search_distance=0.25,
        )
    )

    assert len(set(anchors)) == 6
    assert frame_count == 2
    assert score == pytest.approx(0.0, abs=1e-5)
    assert deviations["distance_nm"] == pytest.approx(0.0, abs=1e-6)
