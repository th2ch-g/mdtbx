"""Extract one named molecule from a multi-record SDF."""

import argparse
import json
from pathlib import Path

from rdkit import Chem

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "extract_ligand",
        help="Extract a named ligand from an SDF for parameterization",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Multi-record SDF")
    parser.add_argument("--name", required=True, help="SDF _Name value")
    parser.add_argument("-o", "--outdir", default="ligand", help="Output directory")
    parser.set_defaults(func=run)


def _select(path, name):
    molecules = [
        molecule
        for molecule in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        if molecule is not None
        and molecule.HasProp("_Name")
        and molecule.GetProp("_Name") == name
    ]
    if len(molecules) != 1:
        raise ValueError(
            f"Expected one SDF record named {name}; found {len(molecules)}"
        )
    molecule = molecules[0]
    if molecule.GetNumConformers() != 1:
        raise ValueError(f"Ligand {name} must contain exactly one conformer")
    return molecule


def run(args):
    source = Path(args.input).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    outdir = Path(args.outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    molecule = _select(source, args.name)
    sdf = outdir / "source.sdf"
    mol = outdir / "source.mol"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(molecule)
    Chem.MolToMolFile(molecule, str(mol))

    manifest = {
        "schema_version": 1,
        "workflow": "extract-ligand",
        "source": str(source),
        "name": args.name,
        "atoms": molecule.GetNumAtoms(),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "sdf": str(sdf),
        "mol": str(mol),
    }
    (outdir / "extraction_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    LOGGER.info("Extracted %s to %s", args.name, outdir)
