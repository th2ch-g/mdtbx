import argparse
import os
import re
from pathlib import Path

from ..utils.atom_selection_parser import AtomSelector
from ..utils.parse_top import GromacsTopologyParser
from ..utils.common_args import add_selection_arg, add_topology_arg
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def _preprocessor_symbol(output_prefix: str) -> str:
    symbol = re.sub(r"[^A-Z0-9_]", "_", Path(output_prefix).name.upper())
    if not symbol or symbol[0].isdigit():
        symbol = f"_{symbol}"
    return symbol


def add_subcmd(subparsers):
    """
    mdtbx gen_posres -g structure.gro -p topology.top
    """
    parser = subparsers.add_parser(
        "gen_posres",
        help="Generate POSRES",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_topology_arg(parser, help="Topology file (.top)")

    add_selection_arg(
        parser,
        help="Selection for positional restraints. (Custom atom selection language like MDtraj)",
    )

    parser.add_argument(
        "-o",
        "--output_prefix",
        default="posres",
        type=str,
        help="Output file prefix (This also will be constant name)",
    )

    parser.set_defaults(func=run)


def run(args):
    selector = AtomSelector(args.selection)
    parser = GromacsTopologyParser(args.topology)

    with open(args.topology) as f:
        lines = f.readlines()

    topology_path = Path(args.topology).resolve()
    const = _preprocessor_symbol(args.output_prefix)
    force_const = const + "_FC"
    all_moleculetypes = parser.get_all_moleculetypes()
    LOGGER.info(f"all_moleculetypes: {all_moleculetypes}")
    pad_count = 0
    for moleculetypes in all_moleculetypes:
        atoms = parser.get_atoms_in(moleculetypes)
        target_atom_indices = []
        for atom in atoms:
            if selector.eval(atom):
                target_atom_indices.append(atom["index"])

        if len(target_atom_indices) == 0:
            LOGGER.info(f"No target atoms in {moleculetypes}")
            LOGGER.info(f"Skip {moleculetypes}")
            continue

        itp_path = Path(f"{args.output_prefix}_{moleculetypes}.itp")
        itp_path.parent.mkdir(parents=True, exist_ok=True)
        with itp_path.open("w") as f:
            f.write(f"#ifdef {const}\n")
            f.write("[ position_restraints ]\n")
            f.write(f"; {args.selection}\n")
            f.write(";  i funct       fcx        fcy        fcz\n")
            for atom_index in target_atom_indices:
                f.write(f"{atom_index} 1 {force_const} {force_const} {force_const}\n")
            f.write("#endif\n")
        LOGGER.info(f"{itp_path} generated")

        target_insert_linenumber = parser.get_insert_linenumber_in(moleculetypes)
        include_path = os.path.relpath(itp_path.resolve(), topology_path.parent)
        include_directive = f'#include "{Path(include_path).as_posix()}"'
        if not any(line.strip() == include_directive for line in lines):
            lines.insert(
                target_insert_linenumber + 1 + pad_count,
                f"\n{include_directive}\n",
            )
            LOGGER.info(f"{itp_path} inserted")
            pad_count += 1

    # update
    with open(args.topology, "w") as f:
        f.writelines(lines)
