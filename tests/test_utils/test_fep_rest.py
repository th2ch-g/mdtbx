import pytest

from mdtbx.utils.fep_rest import (
    build_fep_rest_schedule,
    hot_global_indices,
    parse_preprocessed_topology,
    prepare_plumed_dual_state,
    restore_plumed_dual_state,
    underline_global_atoms,
    unify_fep_rest_charges,
)


TOPOLOGY = """\
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
C    6 12.0 0.0 A 0.34 0.45
DUM  0 12.0 0.0 A 0.00 0.00

[ moleculetype ]
LIG 3

[ atoms ]
1 C   1 LIG C1 1  0.5 12.0 DUM 0.0 12.0
2 C   1 LIG C2 2 -0.5 12.0 C  -0.2 12.0

[ dihedrals ]
1 2 2 1 1 0.0 2.0 3 0.0 4.0 3

[ system ]
Test

[ molecules ]
LIG 1
"""


def test_schedule_contains_three_alchemical_phases_and_symmetric_rest():
    schedule = build_fep_rest_schedule(7, 300.0, 1200.0)

    assert schedule["coordinates"] == [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    assert schedule["vdw"] == [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
    assert schedule["charge_a"] == [0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0]
    assert schedule["charge_b"] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0]
    assert schedule["rest_temperatures"][0] == pytest.approx(300.0)
    assert schedule["rest_temperatures"][3] == pytest.approx(1200.0)
    assert schedule["rest_temperatures"][-1] == pytest.approx(300.0)


def test_schedule_rejects_too_few_replicas():
    with pytest.raises(ValueError, match="at least 4"):
        build_fep_rest_schedule(3, 300.0, 1200.0)


def test_even_replica_schedule_reaches_requested_rest_temperature():
    schedule = build_fep_rest_schedule(4, 300.0, 900.0)

    assert schedule["rest_temperatures"] == pytest.approx([300.0, 900.0, 900.0, 300.0])


def test_parse_and_underline_topology_atoms():
    parsed = parse_preprocessed_topology(TOPOLOGY)

    assert list(parsed.atoms) == ["LIG"]
    assert parsed.atoms["LIG"][0].type_b == "DUM"
    underlined = underline_global_atoms(TOPOLOGY, {0})
    assert "1 C_" in underlined
    assert hot_global_indices(underlined) == [0]


def test_underline_rejects_partial_selection_of_repeated_molecule():
    repeated = TOPOLOGY.replace("LIG 1\n", "LIG 2\n")

    with pytest.raises(ValueError, match="only some copies"):
        underline_global_atoms(repeated, {0})


def test_charge_unification_uses_separate_disappearing_schedule():
    underlined = underline_global_atoms(TOPOLOGY, {0})

    halfway = unify_fep_rest_charges(
        underlined,
        general_lambda=0.1,
        charge_a_lambda=0.5,
        charge_b_lambda=0.0,
    )
    record = parse_preprocessed_topology(halfway).atoms["LIG"][0]

    assert record.charge_a == pytest.approx(0.25)
    assert record.charge_b == pytest.approx(0.25)


def test_plumed_compatibility_restores_hot_b_state_and_dual_dihedral():
    underlined = underline_global_atoms(TOPOLOGY, {0, 1})
    compatible, records = prepare_plumed_dual_state(underlined)

    assert "MDTBX_DUAL_ATOM" in compatible
    dihedral_line = next(
        line for line in compatible.splitlines() if "MDTBX_DUAL_DIHEDRAL" in line
    )
    assert len(dihedral_line.partition(";")[0].split()) == 8

    restored = restore_plumed_dual_state(compatible, records, scale=0.25)
    atoms = parse_preprocessed_topology(restored).atoms["LIG"]
    assert atoms[0].type_b == "DUM_"
    assert atoms[0].charge_b == pytest.approx(0.0)
    restored_dihedral = next(
        line
        for line in restored.splitlines()
        if line.partition(";")[0].split()[:5] == ["1", "2", "2", "1", "1"]
    )
    parameters = restored_dihedral.partition(";")[0].split()[5:]
    assert parameters == ["0.0", "2.0", "3", "0.0", "1", "3"]


def test_plumed_compatibility_rejects_cmap():
    with pytest.raises(ValueError, match="CMAP"):
        prepare_plumed_dual_state(TOPOLOGY.replace("[ system ]", "[ cmap ]"))


def test_plumed_compatibility_rejects_intermolecular_interactions():
    with pytest.raises(ValueError, match="intermolecular"):
        prepare_plumed_dual_state(
            TOPOLOGY.replace(
                "[ system ]",
                "[ intermolecular_interactions ]",
            )
        )


def test_plumed_compatibility_rejects_dual_state_pairs():
    topology = TOPOLOGY.replace(
        "[ dihedrals ]",
        "[ pairs ]\n1 2 1 0.1 0.2 0.3 0.4\n\n[ dihedrals ]",
    )

    with pytest.raises(ValueError, match="B-state pair"):
        prepare_plumed_dual_state(topology)
