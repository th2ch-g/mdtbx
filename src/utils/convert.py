import argparse
import mdtraj as md
from pymol import cmd
from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx convert -f input_structure -o output_structure
    """
    parser = subparsers.add_parser(
        "convert",
        help="Convert files to other format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-f",
        "--structure",
        required=True,
        type=str,
        help="Structure file (cat be parsed in PyMOL or MDtraj)",
    )

    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Output file path"
    )

    parser.add_argument(
        "--type",
        default="pymol",
        type=str,
        help="Type of conversion",
        choices=["pymol", "mdtraj"],
    )


def run(args):
    if args.type == "pymol":
        cmd.load(args.structure)
        cmd.save(args.output)
    elif args.type == "mdtraj":
        traj = md.load(args.structure)
        traj.save(args.output)
    LOGGER.info(f"{args.output} generated")
