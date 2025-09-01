import argparse
from pymol import cmd, editor

from ..config import *  # NOQA
from ..logger import generate_logger

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


def run(args):
    cmd.load(args.structure, "target")
    for chain in cmd.get_chains("target and polymer.protein"):
        editor.attach_amino_acid(f"first (chain {chain}) and name N", "ace")
        LOGGER.info(f"ACE added to {chain}")
    cmd.set("retain_order", 0)
    cmd.sort()
    cmd.save(f"{args.output_prefix}.pdb")

    with open(f"{args.output_prefix}.pdb") as ref:
        lines = ref.readlines()

    with open(f"{args.output_prefix}.pdb") as f:
        for idx, line in enumerate(f):
            line = line.rstrip()

            if "ACE" in line:
                for target_idx in range(1, 3 + 1):
                    lines[idx] = lines[idx].replace(
                        f"HH3{target_idx}", f" H{target_idx} ", 1
                    )

    with open(f"{args.output_prefix}.pdb", "w") as f:
        f.writelines(lines)
