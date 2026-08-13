"""Combine two PDB coordinate sets without changing atom names or positions."""

import argparse
from pathlib import Path

from ..logger import generate_logger
from ..utils.paths import ensure_output_parent

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "combine_pdb",
        help="Combine two PDB files while preserving coordinates and atom names",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-1", "--pdb1", required=True, help="First PDB file")
    parser.add_argument("-2", "--pdb2", required=True, help="Second PDB file")
    parser.add_argument("-o", "--output", default="combined.pdb")
    parser.add_argument(
        "--amber-histidines",
        action="store_true",
        help="Rename HIS to HID/HIE/HIP from explicit HD1/HE2 atoms",
    )
    parser.set_defaults(func=run)


def _atoms(path):
    records = [
        line
        for line in Path(path).read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    if not records:
        raise ValueError(f"No coordinate records found in {path}")
    return records


def _renumber(records, start):
    output = []
    serial = start
    for line in records:
        # The serial column is fixed at 5 characters; wrap to the 1..99999 PDB
        # range so combined systems past 99999 atoms keep the columns aligned.
        output.append(f"{line[:6]}{(serial - 1) % 99999 + 1:5d}{line[11:]}")
        serial += 1
    return output, serial


def _amber_histidines(records):
    residue_atoms = {}
    for line in records:
        if line[17:20].strip() != "HIS":
            continue
        key = (line[21], line[22:26], line[26])
        residue_atoms.setdefault(key, set()).add(line[12:16].strip())

    names = {}
    for key, atoms in residue_atoms.items():
        if "HD1" in atoms and "HE2" in atoms:
            names[key] = "HIP"
        elif "HD1" in atoms:
            names[key] = "HID"
        elif "HE2" in atoms:
            names[key] = "HIE"

    output = []
    for line in records:
        key = (line[21], line[22:26], line[26])
        resname = names.get(key)
        output.append(f"{line[:17]}{resname:>3}{line[20:]}" if resname else line)
    return output


def run(args):
    first = Path(args.pdb1).expanduser().resolve()
    second = Path(args.pdb2).expanduser().resolve()
    for path in (first, second):
        if not path.is_file():
            raise FileNotFoundError(path)
    first_atoms = _atoms(first)
    second_atoms = _atoms(second)
    if getattr(args, "amber_histidines", False):
        first_atoms = _amber_histidines(first_atoms)
        second_atoms = _amber_histidines(second_atoms)
    first_records, serial = _renumber(first_atoms, 1)
    second_records, _ = _renumber(second_atoms, serial)
    output = ensure_output_parent(args.output)
    output.write_text(
        "\n".join([*first_records, "TER", *second_records, "TER", "END"]) + "\n"
    )
    LOGGER.info("Combined %s and %s into %s", first, second, output)
