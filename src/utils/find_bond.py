import argparse
from pymol import cmd

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx find_bond -s structure -s1 selection -s2 selection -c cutoff -o output
    """
    parser = subparsers.add_parser(
        "find_bond",
        help="Find bond (e.g. CYS-CYS)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Structure file (can be parsed in PyMOL)",
    )

    parser.add_argument(
        "-s1",
        "--selection1",
        type=str,
        default="resn CYS and name SG",
        help="Selection 1 (can be parsed in PyMOL)",
    )

    parser.add_argument(
        "-s2",
        "--selection2",
        default="resn CYS and name SG",
        type=str,
        help="Selection 2 (can be parsed in PyMOL)",
    )

    parser.add_argument(
        "-c",
        "--cutoff",
        default=3.0,
        type=float,
        help="Cutoff distance [angstrom]",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name",
    )


def run(args):
    cmd.load(args.structure, "target")
    pairs = cmd.find_pairs(args.selection1, args.selection2, cutoff=args.cutoff)
    bonds_set = set()
    bonds = []
    for idx1, idx2 in pairs:
        if idx1[1] in bonds_set or idx2[1] in bonds_set:
            continue
        bonds_set.add(idx1[1])
        bonds_set.add(idx2[1])
        tmp = []
        cmd.iterate_state(1, f"index {idx1[1]}", "tmp.append(resi)", space=locals())
        cmd.iterate_state(1, f"index {idx2[1]}", "tmp.append(resi)", space=locals())
        bonds.append(tmp)
    bonds_str_lines = []
    for res1, res2 in bonds:
        bonds_str_lines.append(f"bond {SYSTEM_NAME}.{res1}.SG {SYSTEM_NAME}.{res2}.SG")  # NOQA
    bonds_str = "\n".join(bonds_str_lines)

    if args.output is not None:
        with open(args.output, "w") as f:
            f.write(bonds_str)
        LOGGER.info(f"{args.output} generated")
    else:
        if len(bonds) == 0:
            LOGGER.info("No bond found")
        else:
            LOGGER.info(f"{len(bonds)} bonds found")
            print(bonds_str)
