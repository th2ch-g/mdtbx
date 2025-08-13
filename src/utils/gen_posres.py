import argparse
import mdtraj as md

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

    parser.add_argument("-g", "--gro", required=True, type=str, help="GRO file (.gro)")

    parser.add_argument(
        "-p", "--topology", required=True, type=str, help="Topology file (.top)"
    )

    parser.add_argument(
        "-s",
        "--selection",
        required=True,
        type=str,
        help="Selection for positional restraints. (MDtraj atom selection language)",
    )

    parser.add_argument(
        "-o",
        "--output_prefix",
        default="posres",
        type=str,
        help="Output file prefix (This also will be constant name)",
    )


def run(args):
    # generate posres.itp
    gro = md.load(args.gro)
    target_atom_indices = gro.top.select(args.selection)
    target_atom_indices = [i + 1 for i in target_atom_indices]  # index start from 0

    const = args.output_prefix.upper()
    force_const = const + "_FC"
    with open(args.output_prefix + ".itp", "w") as f:
        f.write(f"#ifdef {const}\n")
        f.write("[ position_restraints ]\n")
        f.write(";  i funct       fcx        fcy        fcz\n")
        for atom_index in target_atom_indices:
            f.write(f"{atom_index} 1 {force_const} {force_const} {force_const} \n")
        f.write("#endif\n")

    # insert posres.itp into topology.top
    # system section treats as global
    with open(args.topology, "a") as f:
        f.write(f'\n#include "{args.output_prefix}.itp"\n')
