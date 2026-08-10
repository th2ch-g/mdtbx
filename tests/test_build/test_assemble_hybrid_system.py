from types import SimpleNamespace

import numpy as np
import pytest

from mdtbx.build import assemble_hybrid_system


def test_run_replaces_ligand_topology_and_coordinates(tmp_path):
    topology = tmp_path / "gmx.top"
    topology.write_text(
        "[ defaults ]\n1 2 yes 0.5 0.8333\n"
        "[ atomtypes ]\nc 6 12 0 A 0.3 0.2\n"
        "[ moleculetype ]\nPRO 3\n[ atoms ]\n1 c 1 PRO C 1 0 12\n"
        "[ moleculetype ]\nOLD 3\n[ atoms ]\n1 c 1 LIG C1 1 0 12\n"
        "2 c 1 LIG C2 2 0 12\n"
        "3 c 1 LIG C3 3 0 12\n"
        "[ system ]\nTest\n[ molecules ]\nPRO 1\nOLD 1\n"
    )
    structure = tmp_path / "gmx.gro"
    structure.write_text(
        "Test\n4\n"
        "    1PRO      C    1   0.000   0.000   0.000\n"
        "    2LIG     C1    2   0.600   0.800   0.800\n"
        "    2LIG     C2    3   0.500   0.700   0.800\n"
        "    2LIG     C3    4   0.600   0.700   0.900\n"
        "1.0 1.0 1.0\n"
    )
    hybrid_topology = tmp_path / "merged.itp"
    hybrid_topology.write_text(
        "[ moleculetype ]\nHYB 3\n[ atoms ]\n"
        "1 c 1 LIG C1 1 0 12\n2 c 1 LIG C2 2 0 12\n"
        "3 c 1 LIG C3 3 0 12\n4 c 1 LIG D1 4 0 12\n"
    )
    hybrid_atomtypes = tmp_path / "ffmerged.itp"
    hybrid_atomtypes.write_text("[ atomtypes ]\nhyb 6 12 0 A 0.3 0.2\n")
    hybrid_structure = tmp_path / "merged.pdb"
    hybrid_structure.write_text(
        "HETATM    1  C1  LIG A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  LIG A   1       0.000   1.000   0.000  1.00  0.00           C\n"
        "HETATM    3  C3  LIG A   1       0.000   0.000   1.000  1.00  0.00           C\n"
        "HETATM    4  D1  LIG A   1       2.000   0.000   0.000  1.00  0.00           C\n"
    )
    outdir = tmp_path / "output"

    assemble_hybrid_system.run(
        SimpleNamespace(
            topology=str(topology),
            structure=str(structure),
            ligand_moltype="OLD",
            ligand_resname="LIG",
            hybrid_topology=str(hybrid_topology),
            hybrid_atomtypes=str(hybrid_atomtypes),
            hybrid_structure=str(hybrid_structure),
            outdir=str(outdir),
        )
    )

    output_topology = (outdir / "dual.top").read_text()
    assert '#include "hybrid_atomtypes.itp"' in output_topology
    assert '#include "hybrid.itp"' in output_topology
    assert "HYB                  1" in output_topology
    output_lines = (outdir / "dual.gro").read_text().splitlines()
    assert output_lines[1] == "5"
    dummy_coordinates = np.asarray(
        [
            float(output_lines[6][start:end])
            for start, end in ((20, 28), (28, 36), (36, 44))
        ]
    )
    assert dummy_coordinates == pytest.approx([0.6, 0.9, 0.8])


def test_replace_gro_rejects_different_ligand_pose(tmp_path):
    structure = tmp_path / "gmx.gro"
    structure.write_text(
        "Test\n3\n"
        "    1LIG     C1    1   0.000   0.000   0.000\n"
        "    1LIG     C2    2   0.100   0.000   0.000\n"
        "    1LIG     C3    3   0.800   0.800   0.000\n"
        "1.0 1.0 1.0\n"
    )
    hybrid_structure = tmp_path / "hybrid.pdb"
    hybrid_structure.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  LIG A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  C3  LIG A   1       0.000   1.000   0.000  1.00  0.00           C\n"
    )

    with pytest.raises(ValueError, match="do not describe the same pose"):
        assemble_hybrid_system._replace_gro(structure, hybrid_structure, "LIG")


def test_hybrid_atomtypes_precede_hybrid_when_ligand_is_first_molecule():
    topology = (
        "[ defaults ]\n1 2 yes 0.5 0.8333\n"
        "[ atomtypes ]\nc 6 12 0 A 0.3 0.2\n"
        "[ moleculetype ]\nOLD 3\n[ atoms ]\n1 c 1 LIG C1 1 0 12\n"
        "[ moleculetype ]\nSOL 3\n[ atoms ]\n1 c 1 SOL C1 1 0 12\n"
        "[ system ]\nTest\n[ molecules ]\nOLD 1\nSOL 1\n"
    )

    output = assemble_hybrid_system._replace_molecule_topology(topology, "OLD", "HYB")

    assert output.index('#include "hybrid_atomtypes.itp"') < output.index(
        '#include "hybrid.itp"'
    )
