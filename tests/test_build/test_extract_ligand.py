import types

from rdkit import Chem
from rdkit.Chem import AllChem


def _molecule(name, smiles):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(molecule, randomSeed=2026)
    molecule.SetProp("_Name", name)
    return molecule


def test_extract_ligand_selects_named_record(tmp_path):
    from mdtbx.build import extract_ligand

    source = tmp_path / "ligands.sdf"
    with Chem.SDWriter(str(source)) as writer:
        writer.write(_molecule("lig_a", "CC"))
        writer.write(_molecule("lig_b", "C[NH3+]"))

    outdir = tmp_path / "selected"
    extract_ligand.run(
        types.SimpleNamespace(input=str(source), name="lig_b", outdir=str(outdir))
    )

    selected = next(
        molecule
        for molecule in Chem.SDMolSupplier(str(outdir / "source.sdf"), removeHs=False)
        if molecule is not None
    )
    assert selected.GetProp("_Name") == "lig_b"
    assert Chem.GetFormalCharge(selected) == 1
    assert (outdir / "source.mol").is_file()
    assert (outdir / "extraction_manifest.json").is_file()


def test_extract_ligand_rejects_missing_name(tmp_path):
    from mdtbx.build import extract_ligand

    source = tmp_path / "ligands.sdf"
    with Chem.SDWriter(str(source)) as writer:
        writer.write(_molecule("lig_a", "CC"))

    try:
        extract_ligand.run(
            types.SimpleNamespace(
                input=str(source),
                name="missing",
                outdir=str(tmp_path / "selected"),
            )
        )
    except ValueError as error:
        assert "found 0" in str(error)
    else:
        raise AssertionError("Missing ligand name was accepted")
