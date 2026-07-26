import argparse

import mdtraj as md
from pymol import cmd

from ..logger import generate_logger
from .paths import ensure_output_parent
from .pymol_session import pymol_session

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
        help="Structure file (can be parsed in PyMOL or MDtraj)",
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

    parser.set_defaults(func=run)


def run(args):
    if args.type not in {"pymol", "mdtraj"}:
        raise ValueError(f"Unsupported conversion type: {args.type}")

    output_path = ensure_output_parent(args.output)

    if args.type == "pymol":
        with pymol_session(cmd, args.structure):
            cmd.save(str(output_path), "target")
    elif args.type == "mdtraj":
        traj = md.load(args.structure)
        traj.save(str(output_path))
    LOGGER.info(f"{output_path} generated")
