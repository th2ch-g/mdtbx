import argparse
import re

from .atom_selection_parser import AtomSelector
from .parse_top import GromacsTopologyParser
from ..config import *  # NOQA
from ..logger import generate_logger

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
        "-p", "--topology", type=str, required=True, help="Topology file (.gro, .pdb)"
    )
    parser.add_argument(
        "-s",
        "--selection",
        type=str,
        required=True,
        help="Selection (Custom atom selection language)",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="output.top", help="Output struxture file"
    )


def run(args):
    selector = AtomSelector(args.selection)
    parser = GromacsTopologyParser(args.topology)

    LOGGER.info(parser.all_moleculetypes)

    all_moleculetypes = parser.get_all_moleculetypes()
    LOGGER.info(f"all_moleculetypes: {all_moleculetypes}")

    selected_atoms = []

    for moleculetype in all_moleculetypes:
        LOGGER.info(f"moleculetype: {moleculetype}")
        atoms = parser.get_atoms_in(moleculetype)
        for atom in atoms:
            if selector.eval(atom):
                LOGGER.info(f"atom: {atom}")
                # rename to ATOM_ from ATOM
                selected_atoms.append(atom)

    with open(args.topology) as f:
        lines = f.readlines()

    for atom in selected_atoms:
        LOGGER.info(f"atom: {atom}")
        atom_name = atom["name"]
        # pattern = fr'^(\s*(?:\S+\s+){4}){atom_name}(\s+)'
        # replacement = fr'\1{atom_name}_\2'
        pattern = rf"^(\s*(?:\S+\s+){{4}}){re.escape(atom_name)}(?=\s|$)"
        replacement = rf"\1{atom_name}_"
        atom_linenumber = atom["linenumber"]
        lines[atom_linenumber] = re.sub(pattern, replacement, lines[atom_linenumber])

    with open(args.output, "w") as f:
        f.writelines(lines)
