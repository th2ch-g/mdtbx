"""Replace an endpoint ligand with a PMX hybrid molecule."""

import argparse
import json
import re
import shutil
from pathlib import Path

import numpy as np

from ..logger import generate_logger

LOGGER = generate_logger(__name__)
_INCLUDE = re.compile(r'^\s*#include\s+"([^"]+)"')


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "assemble_hybrid_system",
        help="Insert a PMX hybrid ligand into an endpoint GROMACS system",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-p", "--topology", required=True)
    parser.add_argument("-c", "--structure", required=True)
    parser.add_argument("--ligand-moltype", required=True)
    parser.add_argument("--ligand-resname", required=True)
    parser.add_argument("--hybrid-topology", required=True)
    parser.add_argument("--hybrid-atomtypes", required=True)
    parser.add_argument("--hybrid-structure", required=True)
    parser.add_argument("-o", "--outdir", default="hybrid_system")
    parser.set_defaults(func=run)


def _section_name(line):
    code = line.partition(";")[0].strip()
    if code.startswith("[") and code.endswith("]"):
        return code.strip("[] ").lower()
    return None


def _molecule_type(path):
    section = None
    for raw in Path(path).read_text().splitlines():
        name = _section_name(raw)
        if name is not None:
            section = name
            continue
        code = raw.partition(";")[0].strip()
        if section == "moleculetype" and code and not code.startswith("#"):
            return code.split()[0]
    raise ValueError(f"No [ moleculetype ] found in {path}")


def _replace_molecule_topology(text, old_name, hybrid_name):
    lines = text.splitlines(keepends=True)
    starts = []
    for index, line in enumerate(lines):
        if _section_name(line) != "moleculetype":
            continue
        for entry in lines[index + 1 :]:
            code = entry.partition(";")[0].strip()
            if code and not code.startswith("#"):
                starts.append((index, code.split()[0]))
                break
    matches = [item for item in starts if item[1] == old_name]
    if len(matches) != 1:
        raise ValueError(f"Expected one ligand moleculetype named {old_name}")
    start = matches[0][0]
    end = len(lines)
    for index in range(start + 1, len(lines)):
        if _section_name(lines[index]) in {"moleculetype", "system"}:
            end = index
            break
    lines[start:end] = ['#include "hybrid.itp"\n\n']

    first_molecule_definition = next(
        index
        for index, line in enumerate(lines)
        if _section_name(line) == "moleculetype"
        or line.strip() == '#include "hybrid.itp"'
    )
    lines.insert(
        first_molecule_definition,
        '#include "hybrid_atomtypes.itp"\n\n',
    )

    section = None
    replacements = 0
    for index, line in enumerate(lines):
        name = _section_name(line)
        if name is not None:
            section = name
            continue
        code, separator, comment = line.partition(";")
        fields = code.split()
        if section == "molecules" and fields and fields[0] == old_name:
            fields[0] = hybrid_name
            suffix = f";{comment}" if separator else ""
            lines[index] = f"{fields[0]:<20s} {fields[1]}{suffix}"
            if line.endswith("\n") and not lines[index].endswith("\n"):
                lines[index] += "\n"
            replacements += 1
    if replacements != 1:
        raise ValueError(f"Expected one [ molecules ] entry for {old_name}")
    return "".join(lines)


def _pdb_atoms(path):
    atoms = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            {
                "resname": line[17:20].strip(),
                "atomname": line[12:16].strip(),
                "x": float(line[30:38]) / 10.0,
                "y": float(line[38:46]) / 10.0,
                "z": float(line[46:54]) / 10.0,
            }
        )
    if not atoms:
        raise ValueError(f"No atoms found in {path}")
    return atoms


def _gro_atoms(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Malformed GRO file: {path}")
    count = int(lines[1].strip())
    if len(lines) != count + 3:
        raise ValueError(f"GRO atom count does not match file length: {path}")
    atoms = []
    for line in lines[2:-1]:
        atoms.append(
            {
                "resid": int(line[:5]),
                "resname": line[5:10].strip(),
                "atomname": line[10:15].strip(),
                "x": float(line[20:28]),
                "y": float(line[28:36]),
                "z": float(line[36:44]),
                "tail": line[44:],
            }
        )
    return lines[0], atoms, lines[-1]


def _align_hybrid_to_endpoint(endpoint, hybrid, max_rmsd_nm=0.05):
    endpoint_by_name = {}
    for atom in endpoint:
        name = atom["atomname"].upper()
        if name in endpoint_by_name:
            raise ValueError(f"Endpoint ligand atom name is not unique: {name}")
        endpoint_by_name[name] = atom

    hybrid_by_name = {}
    for atom in hybrid:
        name = atom["atomname"].upper()
        if name in hybrid_by_name:
            raise ValueError(f"Hybrid ligand atom name is not unique: {name}")
        hybrid_by_name[name] = atom

    missing = sorted(set(endpoint_by_name) - set(hybrid_by_name))
    if missing:
        raise ValueError(
            "Hybrid structure is missing endpoint ligand atoms: " + ", ".join(missing)
        )

    names = sorted(endpoint_by_name)
    source = np.asarray(
        [[hybrid_by_name[name][axis] for axis in ("x", "y", "z")] for name in names]
    )
    target = np.asarray(
        [[endpoint_by_name[name][axis] for axis in ("x", "y", "z")] for name in names]
    )
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    left, _singular_values, right_transpose = np.linalg.svd(covariance)
    rotation = left @ right_transpose
    if np.linalg.det(rotation) < 0.0:
        left[:, -1] *= -1.0
        rotation = left @ right_transpose

    fitted = (source - source_center) @ rotation + target_center
    rmsd_nm = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    if rmsd_nm > max_rmsd_nm:
        raise ValueError(
            "Hybrid and endpoint ligand coordinates do not describe the same pose: "
            f"alignment RMSD {rmsd_nm:.4f} nm exceeds {max_rmsd_nm:.4f} nm"
        )

    aligned = []
    for atom in hybrid:
        coordinates = np.asarray([atom[axis] for axis in ("x", "y", "z")])
        coordinates = (coordinates - source_center) @ rotation + target_center
        aligned.append(
            {
                **atom,
                "x": float(coordinates[0]),
                "y": float(coordinates[1]),
                "z": float(coordinates[2]),
            }
        )
    return aligned, rmsd_nm


def _replace_gro(path, hybrid_structure, resname):
    title, atoms, box = _gro_atoms(path)
    selected = [index for index, atom in enumerate(atoms) if atom["resname"] == resname]
    if not selected:
        raise ValueError(f"Ligand residue {resname} not found in GRO structure")
    residue_ids = {atoms[index]["resid"] for index in selected}
    if len(residue_ids) != 1 or selected != list(range(selected[0], selected[-1] + 1)):
        raise ValueError("Ligand selection must be one contiguous residue")
    hybrid = _pdb_atoms(hybrid_structure)
    endpoint = [atoms[index] for index in selected]
    hybrid, alignment_rmsd_nm = _align_hybrid_to_endpoint(endpoint, hybrid)
    ligand_resid = atoms[selected[0]]["resid"]
    replacement = [
        {
            **atom,
            "resid": ligand_resid,
            "tail": "",
        }
        for atom in hybrid
    ]
    merged = atoms[: selected[0]] + replacement + atoms[selected[-1] + 1 :]
    # A GRO file must carry velocities on every atom line or on none; the
    # replacement atoms have none, so drop the velocity columns everywhere.
    if any(atom["tail"].strip() for atom in merged):
        LOGGER.info("Dropped input velocities to keep the merged GRO consistent")
        merged = [{**atom, "tail": ""} for atom in merged]
    output = [title, str(len(merged))]
    for index, atom in enumerate(merged, start=1):
        output.append(
            f"{atom['resid'] % 100000:5d}{atom['resname'][:5]:<5s}"
            f"{atom['atomname'][:5]:>5s}{index % 100000:5d}"
            f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}{atom['tail']}"
        )
    output.append(box)
    return (
        "\n".join(output) + "\n",
        len(selected),
        len(hybrid),
        alignment_rmsd_nm,
    )


def _copy_includes(text, source_dir, outdir):
    for line in text.splitlines():
        match = _INCLUDE.match(line)
        if not match:
            continue
        relative = Path(match.group(1))
        if relative.is_absolute() or relative.name in {
            "hybrid.itp",
            "hybrid_atomtypes.itp",
        }:
            continue
        source = source_dir / relative
        if not source.is_file():
            raise FileNotFoundError(f"Topology include not found: {source}")
        target = outdir / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def run(args):
    topology = Path(args.topology).expanduser().resolve()
    structure = Path(args.structure).expanduser().resolve()
    hybrid_topology = Path(args.hybrid_topology).expanduser().resolve()
    hybrid_atomtypes = Path(args.hybrid_atomtypes).expanduser().resolve()
    hybrid_structure = Path(args.hybrid_structure).expanduser().resolve()
    for path in (
        topology,
        structure,
        hybrid_topology,
        hybrid_atomtypes,
        hybrid_structure,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    outdir = Path(args.outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    hybrid_name = _molecule_type(hybrid_topology)
    output_topology = _replace_molecule_topology(
        topology.read_text(), args.ligand_moltype, hybrid_name
    )
    output_gro, endpoint_atoms, hybrid_atoms, alignment_rmsd_nm = _replace_gro(
        structure, hybrid_structure, args.ligand_resname
    )
    _copy_includes(output_topology, topology.parent, outdir)
    shutil.copy2(hybrid_topology, outdir / "hybrid.itp")
    shutil.copy2(hybrid_atomtypes, outdir / "hybrid_atomtypes.itp")
    (outdir / "dual.top").write_text(output_topology)
    (outdir / "dual.gro").write_text(output_gro)
    manifest = {
        "schema_version": 1,
        "workflow": "assemble-hybrid-system",
        "source_moltype": args.ligand_moltype,
        "hybrid_moltype": hybrid_name,
        "endpoint_ligand_atoms": endpoint_atoms,
        "hybrid_ligand_atoms": hybrid_atoms,
        "coordinate_alignment_rmsd_nm": alignment_rmsd_nm,
        "topology": str(outdir / "dual.top"),
        "structure": str(outdir / "dual.gro"),
    }
    (outdir / "assembly_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    LOGGER.info("Assembled hybrid system in %s", outdir)
