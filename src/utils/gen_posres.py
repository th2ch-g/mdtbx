import argparse

from .atom_selection_parser import AtomSelector
from .parse_top import GromacsTopologyParser
from ..config import *  # NOQA
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

    # parser.add_argument("-g", "--gro", required=True, type=str, help="GRO file (.gro)")

    parser.add_argument(
        "-p", "--topology", required=True, type=str, help="Topology file (.top)"
    )

    parser.add_argument(
        "-s",
        "--selection",
        required=True,
        type=str,
        help="Selection for positional restraints. (Custom atom selection language like MDtraj)",
    )

    parser.add_argument(
        "-o",
        "--output_prefix",
        default="posres",
        type=str,
        help="Output file prefix (This also will be constant name)",
    )


def run(args):
    selector = AtomSelector(args.selection)
    parser = GromacsTopologyParser(args.topology)

    print(parser.all_moleculetypes)

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

    # global atom id cannot be used
    # # generate posres.itp
    # gro = md.load(args.gro)
    # target_atom_indices = gro.top.select(args.selection)
    # target_atom_indices = [i + 1 for i in target_atom_indices]  # index start from 0
    #
    # const = args.output_prefix.upper()
    # force_const = const + "_FC"
    # with open(args.output_prefix + ".itp", "w") as f:
    #     f.write(f"#ifdef {const}\n")
    #     f.write("[ position_restraints ]\n")
    #     f.write(";  i funct       fcx        fcy        fcz\n")
    #     for atom_index in target_atom_indices:
    #         f.write(f"{atom_index} 1 {force_const} {force_const} {force_const} \n")
    #     f.write("#endif\n")
    #
    # # insert posres.itp into topology.top
    # # system section treats as global
    # with open(args.topology, "a") as f:
    #     f.write(f'\n#include "{args.output_prefix}.itp"\n')
