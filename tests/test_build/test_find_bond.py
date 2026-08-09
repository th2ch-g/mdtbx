from mdtbx.build.find_bond import _tleap_residue_numbers


def test_tleap_residue_numbers_are_sequential_across_chains():
    records = [
        (1, "", "A", "10"),
        (2, "", "A", "10"),
        (3, "", "A", "11"),
        (4, "", "B", "10"),
        (5, "", "B", "10"),
    ]

    assert _tleap_residue_numbers(records) == {
        1: 1,
        2: 1,
        3: 2,
        4: 3,
        5: 3,
    }


def test_tleap_residue_numbers_follow_atom_order_not_pdb_number():
    records = [(4, "", "B", "1"), (1, "", "A", "100"), (2, "", "A", "100")]

    assert _tleap_residue_numbers(records) == {1: 1, 2: 1, 4: 2}
