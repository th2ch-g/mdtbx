import argparse
import subprocess

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx addh -s input.pdb -o input_h --method reduce
    """
    parser = subparsers.add_parser(
        "addh",
        help="Add hydrogen atoms to structure",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Input structure file (PDB etc.)",
    )

    parser.add_argument(
        "-o",
        "--output_prefix",
        default="out_h",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "--method",
        default="reduce",
        choices=["reduce", "pymol"],
        help="Method for hydrogen addition",
    )

    parser.set_defaults(func=run)


def run(args):
    output_pdb = f"{args.output_prefix}.pdb"

    if args.method == "reduce":
        cmd = f"reduce -build {args.structure} > {output_pdb}"
        subprocess.run(cmd, shell=True, check=True)
    else:
        from pymol import cmd as pymol_cmd

        pymol_cmd.load(args.structure, "target")
        pymol_cmd.h_add("target")
        pymol_cmd.save(output_pdb, "target")

    LOGGER.info(f"{output_pdb} generated")
