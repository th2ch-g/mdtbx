import pytest

from mdtbx.build.pdb_caps import normalize_methyl_hydrogen_names


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("ATOM      1 1HH3 ACE A   1\n", "ATOM      1  H1  ACE A   1\n"),
        ("ATOM      2 2HH3 NME A   2\n", "ATOM      2  H2  NME A   2\n"),
        ("ATOM      3 HH33 ACE A   1\n", "ATOM      3  H3  ACE A   1\n"),
    ],
)
def test_normalize_methyl_hydrogen_names(source, expected):
    assert normalize_methyl_hydrogen_names(source) == expected


def test_unrelated_atom_name_is_unchanged():
    line = "ATOM      1  CA  ALA A   1\n"

    assert normalize_methyl_hydrogen_names(line) == line
