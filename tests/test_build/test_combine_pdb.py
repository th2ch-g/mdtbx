import types


def test_combine_pdb_preserves_atom_fields_and_coordinates(tmp_path):
    from mdtbx.build import combine_pdb

    first = tmp_path / "protein.pdb"
    second = tmp_path / "ligand.pdb"
    output = tmp_path / "complex.pdb"
    first.write_text(
        "ATOM     90  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n"
    )
    second.write_text(
        "HETATM  900  C1  LIG L   1       4.000   5.000   6.000  1.00  0.00           C\nEND\n"
    )

    combine_pdb.run(
        types.SimpleNamespace(pdb1=str(first), pdb2=str(second), output=str(output))
    )

    lines = output.read_text().splitlines()
    assert lines[0].startswith("ATOM      1  CA  ALA")
    assert lines[2].startswith("HETATM    2  C1  LIG")
    assert "   4.000   5.000   6.000" in lines[2]
    assert lines[-1] == "END"


def test_renumber_wraps_serials_at_pdb_column_limit():
    from mdtbx.build import combine_pdb

    line = (
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C"
    )
    records, next_serial = combine_pdb._renumber([line, line], 99999)

    assert records[0][6:11] == "99999"
    # 100000 does not fit the 5-column serial field; it must wrap to 1
    # instead of shifting every downstream column by one character.
    assert records[1][6:11] == "    1"
    assert all(len(record) == len(line) for record in records)
    assert next_serial == 100001


def test_combine_pdb_renames_explicit_amber_histidines(tmp_path):
    from mdtbx.build import combine_pdb

    first = tmp_path / "protein.pdb"
    second = tmp_path / "ligand.pdb"
    output = tmp_path / "complex.pdb"
    first.write_text(
        "ATOM      1  AND1 HIS A   7       1.000   2.000   3.000  1.00  0.00           N\n"
        "ATOM      2  HD1 HIS A   7       1.000   2.000   3.000  1.00  0.00           H\n"
        "ATOM      3  NE2 HIS A   8       1.000   2.000   3.000  1.00  0.00           N\n"
        "ATOM      4  HE2 HIS A   8       1.000   2.000   3.000  1.00  0.00           H\n"
        "END\n"
    )
    second.write_text(
        "HETATM    5  C1  LIG L   1       4.000   5.000   6.000  1.00  0.00           C\nEND\n"
    )

    combine_pdb.run(
        types.SimpleNamespace(
            pdb1=str(first),
            pdb2=str(second),
            output=str(output),
            amber_histidines=True,
        )
    )

    lines = output.read_text().splitlines()
    assert lines[0][17:20] == "HID"
    assert lines[2][17:20] == "HIE"
