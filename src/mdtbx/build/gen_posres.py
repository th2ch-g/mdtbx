import argparse

from ..utils.atom_selection_parser import AtomSelector
from ..utils.parse_top import GromacsTopologyParser
from ..utils.common_args import add_selection_arg, add_topology_arg
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


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

    const = args.output_prefix.upper()
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

        with open(f"{args.output_prefix}_{moleculetypes}.itp", "w") as f:
            f.write(f"#ifdef {const}\n")
            f.write("[ position_restraints ]\n")
            f.write(f"; {args.selection}\n")
            f.write(";  i funct       fcx        fcy        fcz\n")
            for atom_index in target_atom_indices:
                f.write(f"{atom_index} 1 {force_const} {force_const} {force_const}\n")
            f.write("#endif\n")
        LOGGER.info(f"{args.output_prefix}_{moleculetypes}.itp generated")

        target_insert_linenumber = parser.get_insert_linenumber_in(moleculetypes)
        lines.insert(
            target_insert_linenumber + pad_count,
            f'\n#include "{args.output_prefix}_{moleculetypes}.itp"\n',
        )
        LOGGER.info(f"{args.output_prefix}_{moleculetypes}.itp inserted")

        pad_count += 1

    # update
    with open(args.topology, "w") as f:
        f.writelines(lines)
