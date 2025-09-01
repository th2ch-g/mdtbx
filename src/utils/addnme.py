import argparse
from pymol import cmd, editor

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "addnme",
        help="Add NME to protein",
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
        "-o", "--output_prefix", default="out_nme", type=str, help="Output file prefix"
    )


def run(args):
    cmd.load(args.structure, "target")
    cmd.select("oxts", "name OXT")
    cmd.remove("oxts")
    LOGGER.info("OXTs removed")
    for chain in cmd.get_chains("target and polymer.protein"):
        editor.attach_amino_acid(f"last (chain {chain}) and name C", "nme")
        LOGGER.info(f"NME added to {chain}")
    cmd.save(f"{args.output_prefix}.pdb")

    with open(f"{args.output_prefix}.pdb") as ref:
        lines = ref.readlines()

    with open(f"{args.output_prefix}.pdb") as f:
        for idx, line in enumerate(f):
            line = line.rstrip()

            if "NME" in line:
                for target_idx in range(1, 3 + 1):
                    lines[idx] = lines[idx].replace(
                        f"HH3{target_idx}", f" H{target_idx} ", 1
                    )

                lines[idx] = lines[idx].replace("CH3", "C  ", 1)

    with open(f"{args.output_prefix}.pdb", "w") as f:
        f.writelines(lines)
