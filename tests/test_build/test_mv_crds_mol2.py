import types

import pytest

from src.build.mv_crds_mol2 import run


def _mol2(atom_lines):
    return (
        "@<TRIPOS>MOLECULE\n"
        "LIG\n"
        "@<TRIPOS>ATOM\n" + "".join(atom_lines) + "@<TRIPOS>BOND\n"
    )


def test_run_replaces_coordinates_by_atom_name(tmp_path):
    reference = tmp_path / "reference.mol2"
    coordinates = tmp_path / "coordinates.mol2"
    output = tmp_path / "output.mol2"
    reference.write_text(
        _mol2(
            [
                "1 C1 0.000 0.000 0.000 C.3 1 LIG\n",
                "\n",
                "2 O1 1.000 0.000 0.000 O.2 1 LIG\n",
            ]
        )
    )
    coordinates.write_text(
        _mol2(
            [
                "1 O1 4.000 5.000 6.000 O.2 1 LIG\n",
                "2 C1 1.000 2.000 3.000 C.3 1 LIG\n",
            ]
        )
    )

    run(
        types.SimpleNamespace(
            reference=str(reference),
            coordinates=str(coordinates),
            output=str(output),
        )
    )

    text = output.read_text()
    assert "C1 1.000 2.000 3.000" in text
    assert "O1 4.000 5.000 6.000" in text


def test_run_rejects_duplicate_atom_names(tmp_path):
    reference = tmp_path / "reference.mol2"
    coordinates = tmp_path / "coordinates.mol2"
    reference.write_text(_mol2(["1 C1 0 0 0 C.3 1 LIG\n"]))
    coordinates.write_text(
        _mol2(
            [
                "1 C1 0 0 0 C.3 1 LIG\n",
                "2 C1 1 1 1 C.3 1 LIG\n",
            ]
        )
    )

    with pytest.raises(ValueError, match="Duplicate atom name"):
        run(
            types.SimpleNamespace(
                reference=str(reference),
                coordinates=str(coordinates),
                output=str(tmp_path / "output.mol2"),
            )
        )


def test_run_rejects_mismatched_atom_sets(tmp_path):
    reference = tmp_path / "reference.mol2"
    coordinates = tmp_path / "coordinates.mol2"
    reference.write_text(_mol2(["1 C1 0 0 0 C.3 1 LIG\n"]))
    coordinates.write_text(_mol2(["1 O1 0 0 0 O.2 1 LIG\n"]))

    with pytest.raises(ValueError, match="missing from coordinate"):
        run(
            types.SimpleNamespace(
                reference=str(reference),
                coordinates=str(coordinates),
                output=str(tmp_path / "output.mol2"),
            )
        )


def test_run_can_match_duplicate_names_by_atom_index(tmp_path):
    reference = tmp_path / "reference.mol2"
    coordinates = tmp_path / "coordinates.mol2"
    output = tmp_path / "output.mol2"
    reference.write_text(
        _mol2(
            [
                "1 C 0 0 0 C.3 1 LIG\n",
                "2 C 0 0 0 C.3 1 LIG\n",
            ]
        )
    )
    coordinates.write_text(
        _mol2(
            [
                "2 C 4 5 6 C.3 1 LIG\n",
                "1 C 1 2 3 C.3 1 LIG\n",
            ]
        )
    )

    run(
        types.SimpleNamespace(
            reference=str(reference),
            coordinates=str(coordinates),
            output=str(output),
            match_by="index",
        )
    )

    atom_lines = [
        line
        for line in output.read_text().splitlines()
        if line.startswith(("1 ", "2 "))
    ]
    assert "1 C 1 2 3" in atom_lines[0]
    assert "2 C 4 5 6" in atom_lines[1]
