import pytest
from rdkit import Chem
from rdkit.Geometry import Point3D

from mdtbx.build import pmx_ligand_hybrid


def _molecule():
    molecule = Chem.AddHs(Chem.MolFromSmiles("CO"))
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, Point3D(float(index), 0.0, 0.0))
    molecule.AddConformer(conformer)
    return molecule


def test_atom_order_tracks_parameterization_reordering():
    source = _molecule()
    permutation = list(reversed(range(source.GetNumAtoms())))
    parameterized = Chem.RenumberAtoms(source, permutation)

    order, rmsd = pmx_ligand_hybrid._atom_order(source, parameterized)

    assert [permutation[index] for index in order] == list(range(source.GetNumAtoms()))
    assert rmsd == pytest.approx(0.0)


def test_hybrid_charges_support_common_and_dual_atoms(tmp_path):
    topology = tmp_path / "hybrid.itp"
    topology.write_text(
        "[ moleculetype ]\nLIG 3\n[ atoms ]\n"
        "1 c 1 LIG C1 1 -0.5 12.0 c -0.4 12.0\n"
        "2 h 1 LIG H1 2 -0.5 1.0\n"
    )

    charge_a, charge_b = pmx_ligand_hybrid._hybrid_charges(topology)

    assert charge_a == pytest.approx(-1.0)
    assert charge_b == pytest.approx(-0.9)


def test_gaff_mol2_reader_supports_amber_atom_types(tmp_path):
    mol2 = tmp_path / "ligand.mol2"
    mol2.write_text(
        "@<TRIPOS>MOLECULE\nLIG\n2 1 0 0 0\nSMALL\nbcc\n"
        "@<TRIPOS>ATOM\n"
        "1 CL1 0.0 0.0 0.0 cl 1 LIG -0.1\n"
        "2 H1 1.0 0.0 0.0 ha 1 LIG 0.1\n"
        "@<TRIPOS>BOND\n1 1 2 1\n"
    )

    molecule = pmx_ligand_hybrid._molecule(mol2)

    assert [atom.GetSymbol() for atom in molecule.GetAtoms()] == ["Cl", "H"]
    assert molecule.GetNumBonds() == 1


def test_pdb_elements_reads_chlorine_without_element_column(tmp_path):
    pdb = tmp_path / "ligand.pdb"
    pdb.write_text(
        "ATOM      1  CL1 LIG     1       0.000   0.000   0.000  1.00  0.00\n"
        "ATOM      2  C1  LIG     1       1.000   0.000   0.000  1.00  0.00\n"
    )

    elements = pmx_ligand_hybrid._pdb_elements(pdb)

    assert elements == [17, 6]
