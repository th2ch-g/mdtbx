import argparse

from pymol import cmd, editor

from ..logger import generate_logger
from ..utils.pymol_session import pymol_session
from .pdb_caps import normalize_methyl_hydrogen_names

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "addace",
        help="Add ACE to protein",
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
        "-o", "--output_prefix", default="out_ace", type=str, help="Output file prefix"
    )

    parser.set_defaults(func=run)


def run(args):
    with pymol_session(cmd, args.structure):
        for chain in cmd.get_chains("target and polymer.protein"):
            if chain:
                selection = f"first (chain {chain}) and name N"
            else:
                selection = "first polymer.protein and name N"
            editor.attach_amino_acid(selection, "ace")
            LOGGER.info(f"ACE added to chain '{chain}'")
        cmd.set("retain_order", 0)
        cmd.sort()
        cmd.save(f"{args.output_prefix}.pdb")

    with open(f"{args.output_prefix}.pdb") as ref:
        lines = [
            normalize_methyl_hydrogen_names(line) if "ACE" in line else line
            for line in ref
        ]

    with open(f"{args.output_prefix}.pdb", "w") as f:
        f.writelines(lines)
