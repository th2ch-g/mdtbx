import argparse
import re
from pathlib import Path

from ..logger import generate_logger
from .atom_selection_parser import AtomSelector
from .parse_top import GromacsTopologyParser

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx partial_tempering -s $SELECTION_TEMPERING -p rep${rep}/gmx_pre.top -o rep${rep}/gmx_pre2.top
    """
    parser = subparsers.add_parser(
        "partial_tempering",
        help="Add _ to selected atom groups in gromacs topology file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-p", "--topology", type=str, required=True, help="GROMACS topology file (.top)"
    )
    parser.add_argument(
        "-s",
        "--selection",
        type=str,
        required=True,
        help="Selection (Custom atom selection language)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.top",
        help="Output topology file (.top)",
    )

    parser.set_defaults(func=run)


def run(args):
    selector = AtomSelector(args.selection)
    parser = GromacsTopologyParser(args.topology)

    all_moleculetypes = parser.get_all_moleculetypes()
    selected_atoms: list[dict[str, str | int]] = []

    for moleculetype in all_moleculetypes:
        atoms = parser.get_atoms_in(moleculetype)
        for atom in atoms:
            if selector.eval(atom):
                selected_atoms.append(atom)

    topology_path = Path(args.topology)
    output_path = Path(args.output)
    with topology_path.open() as f:
        lines = f.readlines()

    updated_count = 0
    for atom in selected_atoms:
        atom_type = atom["atom_type"]
        if not isinstance(atom_type, str) or atom_type.endswith("_"):
            continue
        pattern = rf"^(\s*(?:\S+\s+){{1}}){re.escape(atom_type)}(?=\s|$)"
        replacement = rf"\g<1>{atom_type}_"
        atom_linenumber = atom["linenumber"]
        if not isinstance(atom_linenumber, int):
            raise TypeError("Atom line number must be an integer")
        lines[atom_linenumber], replacements = re.subn(
            pattern, replacement, lines[atom_linenumber], count=1
        )
        updated_count += replacements

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        f.writelines(lines)

    LOGGER.info("Updated %d selected atoms in %s", updated_count, output_path)
