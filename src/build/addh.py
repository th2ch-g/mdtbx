import argparse

from ..logger import generate_logger
from ..utils.proc import run_cmd

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
        # Avoid shell=True: write reduce's stdout straight to the output file so
        # paths with spaces or shell metacharacters cannot break or inject.
        with open(output_pdb, "w") as out:
            run_cmd(["reduce", "-build", args.structure], stdout=out)
    else:
        from pymol import cmd as pymol_cmd

        from ..utils.pymol_session import pymol_session

        with pymol_session(pymol_cmd, args.structure) as pymol_cmd:
            pymol_cmd.h_add("target")
            pymol_cmd.save(output_pdb, "target")

    LOGGER.info(f"{output_pdb} generated")
