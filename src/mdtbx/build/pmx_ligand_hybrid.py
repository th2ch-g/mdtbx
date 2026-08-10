"""Generate a validated PMX ligand hybrid topology."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Geometry import Point3D

from ..logger import generate_logger
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "pmx_ligand_hybrid",
        help="Generate a PMX ligand hybrid topology from an OpenFE atom mapping",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--source-a", required=True, help="Original state-A SDF")
    parser.add_argument("--source-b", required=True, help="Original state-B SDF")
    parser.add_argument("--parameterized-a", required=True, help="State-A GAFF2 MOL2")
    parser.add_argument("--parameterized-b", required=True, help="State-B GAFF2 MOL2")
    parser.add_argument("--structure-a", required=True, help="State-A PMX input PDB")
    parser.add_argument("--structure-b", required=True, help="State-B PMX input PDB")
    parser.add_argument("--topology-a", required=True, help="State-A GROMACS ITP")
    parser.add_argument("--topology-b", required=True, help="State-B GROMACS ITP")
    parser.add_argument("--mapping", required=True, help="OpenFE edge YAML")
    parser.add_argument("--edge", required=True, help="Edge key in the mapping YAML")
    parser.add_argument("--charge-a", type=float, required=True)
    parser.add_argument("--charge-b", type=float, required=True)
    parser.add_argument("--pmx", default="pmx", help="PMX executable")
    parser.add_argument("-o", "--outdir", default="pmx_hybrid")
    parser.set_defaults(func=run)


def _molecule(path):
    path = Path(path).expanduser().resolve()
    suffix = path.suffix.lower()
    if suffix in {".sdf", ".mol"}:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        molecules = [item for item in supplier if item is not None]
        if len(molecules) != 1:
            raise ValueError(f"Expected exactly one molecule in {path}")
        return molecules[0]
    if suffix == ".mol2":
        molecule = _gaff_mol2_molecule(path)
    elif suffix == ".pdb":
        molecule = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)
    else:
        raise ValueError(f"Unsupported molecule format: {path}")
    if molecule is None:
        raise ValueError(f"Could not read molecule: {path}")
    return molecule


def _gaff_element(atom_name, atom_type):
    token = atom_type.partition(".")[0].lower()
    if token.startswith("cl"):
        return "Cl"
    if token.startswith("br"):
        return "Br"
    symbols = {
        "c": "C",
        "h": "H",
        "n": "N",
        "o": "O",
        "s": "S",
        "p": "P",
        "f": "F",
        "i": "I",
    }
    if token[:1] in symbols:
        return symbols[token[0]]
    name = "".join(character for character in atom_name if character.isalpha())
    if name[:2].capitalize() in {"Cl", "Br"}:
        return name[:2].capitalize()
    if name[:1].upper() in symbols.values():
        return name[:1].upper()
    raise ValueError(f"Could not infer element for MOL2 atom {atom_name} ({atom_type})")


def _gaff_mol2_molecule(path):
    """Read Amber/GAFF MOL2 files whose atom types are not SYBYL types."""
    atoms = []
    bonds = []
    section = None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if line.startswith("@<TRIPOS>"):
            section = line.removeprefix("@<TRIPOS>")
            continue
        if not line:
            continue
        fields = line.split()
        if section == "ATOM":
            if len(fields) < 6:
                raise ValueError(f"Malformed MOL2 atom line: {raw}")
            atoms.append(
                (
                    int(fields[0]),
                    _gaff_element(fields[1], fields[5]),
                    tuple(float(value) for value in fields[2:5]),
                )
            )
        elif section == "BOND":
            if len(fields) < 4:
                raise ValueError(f"Malformed MOL2 bond line: {raw}")
            bonds.append((int(fields[1]), int(fields[2])))
    if not atoms:
        raise ValueError(f"No atoms found in MOL2 file: {path}")

    atom_numbers = [number for number, _element, _xyz in atoms]
    if atom_numbers != list(range(1, len(atoms) + 1)):
        raise ValueError("MOL2 atom numbers must be consecutive and one-based")
    editable = Chem.RWMol()
    conformer = Chem.Conformer(len(atoms))
    for index, (_number, element, xyz) in enumerate(atoms):
        atom = Chem.Atom(element)
        atom.SetNoImplicit(True)
        editable.AddAtom(atom)
        conformer.SetAtomPosition(index, Point3D(*xyz))
    for index_a, index_b in bonds:
        editable.AddBond(index_a - 1, index_b - 1, Chem.BondType.SINGLE)
    molecule = editable.GetMol()
    molecule.AddConformer(conformer)
    return molecule


def _atom_order(source, parameterized):
    if source.GetNumAtoms() != parameterized.GetNumAtoms():
        raise ValueError("Parameterization changed the ligand atom count")
    matches = parameterized.GetSubstructMatches(
        source,
        uniquify=False,
        useChirality=True,
        maxMatches=10000,
    )
    if not matches:
        mcs = rdFMCS.FindMCS(
            [source, parameterized],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            timeout=30,
        )
        if mcs.canceled or mcs.numAtoms != source.GetNumAtoms():
            raise ValueError("Could not map the source ligand onto the GAFF2 ligand")
        query = Chem.MolFromSmarts(mcs.smartsString)
        source_matches = source.GetSubstructMatches(
            query, uniquify=False, maxMatches=10000
        )
        target_matches = parameterized.GetSubstructMatches(
            query, uniquify=False, maxMatches=10000
        )
        matches = []
        for source_match in source_matches:
            for target_match in target_matches:
                mapping = [None] * source.GetNumAtoms()
                for query_index, source_index in enumerate(source_match):
                    mapping[source_index] = target_match[query_index]
                if all(index is not None for index in mapping):
                    matches.append(tuple(mapping))
        if not matches:
            raise ValueError("Could not map the source ligand onto the GAFF2 ligand")
    matches = sorted(set(matches))
    source_xyz = np.asarray(source.GetConformer().GetPositions())
    target_xyz = np.asarray(parameterized.GetConformer().GetPositions())
    scored = []
    for match in matches:
        delta = source_xyz - target_xyz[np.asarray(match)]
        scored.append((float(np.mean(np.sum(delta * delta, axis=1))), tuple(match)))
    scored.sort()
    if len(scored) > 1 and math.isclose(
        scored[0][0], scored[1][0], rel_tol=0.0, abs_tol=1e-10
    ):
        raise ValueError("Atom-order mapping is ambiguous after parameterization")
    return scored[0][1], math.sqrt(scored[0][0])


def _openfe_mapping(path, edge):
    data = yaml.safe_load(Path(path).read_text())
    try:
        entry = data["edges"][edge]
        mapping = entry["atom mapping"]
    except (KeyError, TypeError) as error:
        raise ValueError(f"OpenFE edge not found: {edge}") from error
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError("OpenFE atom mapping must be a non-empty mapping")
    return entry, {int(key): int(value) for key, value in mapping.items()}


def _itp_atoms(path):
    atoms = []
    section = None
    for raw in Path(path).read_text().splitlines():
        code = raw.partition(";")[0].strip()
        if not code:
            continue
        if code.startswith("[") and code.endswith("]"):
            section = code.strip("[] ").lower()
            continue
        if section == "atoms" and not code.startswith("#"):
            atoms.append(code.split())
    if not atoms:
        raise ValueError(f"No [ atoms ] entries found in {path}")
    return atoms


def _pdb_elements(path):
    elements = []
    for raw in Path(path).read_text().splitlines():
        if not raw.startswith(("ATOM  ", "HETATM")):
            continue
        element = raw[76:78].strip().capitalize()
        if not element:
            atom_name = "".join(
                character for character in raw[12:16] if character.isalpha()
            )
            upper = atom_name.upper()
            if upper.startswith("CL"):
                element = "Cl"
            elif upper.startswith("BR"):
                element = "Br"
            elif upper:
                element = upper[0]
        try:
            elements.append(Chem.GetPeriodicTable().GetAtomicNumber(element))
        except RuntimeError as error:
            raise ValueError(f"Could not infer PDB element from line: {raw}") from error
    if not elements:
        raise ValueError(f"No atoms found in PDB file: {path}")
    return elements


def _validate_endpoint(structure, topology, parameterized):
    pdb_elements = _pdb_elements(structure)
    if len(pdb_elements) != parameterized.GetNumAtoms():
        raise ValueError("PMX input structure atom count differs from GAFF2 MOL2")
    parameterized_elements = [atom.GetAtomicNum() for atom in parameterized.GetAtoms()]
    if pdb_elements != parameterized_elements:
        raise ValueError("PMX input structure atom order differs from GAFF2 MOL2")
    if len(_itp_atoms(topology)) != parameterized.GetNumAtoms():
        raise ValueError("PMX input topology atom count differs from GAFF2 MOL2")


def _hybrid_charges(path):
    charge_a = 0.0
    charge_b = 0.0
    for fields in _itp_atoms(path):
        if len(fields) < 8:
            raise ValueError("Malformed PMX [ atoms ] entry")
        charge_a += float(fields[6])
        charge_b += float(fields[9]) if len(fields) >= 11 else float(fields[6])
    return charge_a, charge_b


def run(args):
    outdir = Path(args.outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    source_a = _molecule(args.source_a)
    source_b = _molecule(args.source_b)
    parameterized_a = _molecule(args.parameterized_a)
    parameterized_b = _molecule(args.parameterized_b)
    order_a, rmsd_a = _atom_order(source_a, parameterized_a)
    order_b, rmsd_b = _atom_order(source_b, parameterized_b)
    _validate_endpoint(args.structure_a, args.topology_a, parameterized_a)
    _validate_endpoint(args.structure_b, args.topology_b, parameterized_b)

    edge, mapping = _openfe_mapping(args.mapping, args.edge)
    if any(index >= len(order_a) for index in mapping):
        raise ValueError("State-A mapping index is outside the ligand")
    if any(index >= len(order_b) for index in mapping.values()):
        raise ValueError("State-B mapping index is outside the ligand")
    translated = [
        (order_a[index_a] + 1, order_b[index_b] + 1)
        for index_a, index_b in mapping.items()
    ]
    pairs = outdir / "pairs.dat"
    pairs.write_text(
        "".join(f"{index_a} {index_b}\n" for index_a, index_b in translated)
    )

    outputs = {
        "structure_a": outdir / "hybrid_a.pdb",
        "structure_b": outdir / "hybrid_b.pdb",
        "topology": outdir / "hybrid.itp",
        "atomtypes": outdir / "hybrid_atomtypes.itp",
        "log": outdir / "pmx_hybrid.log",
    }
    run_cmd(
        [
            args.pmx,
            "ligandHybrid",
            "-i1",
            str(Path(args.structure_a).resolve()),
            "-i2",
            str(Path(args.structure_b).resolve()),
            "-itp1",
            str(Path(args.topology_a).resolve()),
            "-itp2",
            str(Path(args.topology_b).resolve()),
            "-pairs",
            str(pairs),
            "-oA",
            str(outputs["structure_a"]),
            "-oB",
            str(outputs["structure_b"]),
            "-oitp",
            str(outputs["topology"]),
            "-offitp",
            str(outputs["atomtypes"]),
            "-log",
            str(outputs["log"]),
        ]
    )
    for path in outputs.values():
        if not path.is_file():
            raise FileNotFoundError(f"PMX output not found: {path}")
    charge_a, charge_b = _hybrid_charges(outputs["topology"])
    if not math.isclose(charge_a, args.charge_a, abs_tol=1e-4):
        raise ValueError(f"Unexpected PMX state-A charge: {charge_a}")
    if not math.isclose(charge_b, args.charge_b, abs_tol=1e-4):
        raise ValueError(f"Unexpected PMX state-B charge: {charge_b}")
    manifest = {
        "schema_version": 1,
        "workflow": "pmx-ligand-hybrid",
        "edge": args.edge,
        "ligand_a": edge.get("ligand_a"),
        "ligand_b": edge.get("ligand_b"),
        "source_mapping_count": len(mapping),
        "parameterization_rmsd_a_angstrom": rmsd_a,
        "parameterization_rmsd_b_angstrom": rmsd_b,
        "charge_a": charge_a,
        "charge_b": charge_b,
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    (outdir / "hybrid_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Generated validated PMX hybrid topology in %s", outdir)
