import argparse
from pathlib import Path

from ..logger import generate_logger
from ..utils.gmx import remove_gmx_backups
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx trjcat --skip <int> --keep_selection <str> --centering_selection <str> --num_of_step <int>
    """
    parser = subparsers.add_parser(
        "trjcat",
        help="Concatenate trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-k", "--keep_selection", default="non-Water", type=str, help="Keep selection"
    )

    parser.add_argument(
        "-c",
        "--centering_selection",
        default="Protein",
        type=str,
        help="Centering selection",
    )

    parser.add_argument(
        "-n", "--num_of_step", required=True, type=int, help="Number of steps"
    )

    parser.add_argument(
        "-idx", "--index", default="index.ndx", type=str, help="Index file"
    )

    parser.add_argument(
        "--pbc",
        default="mol",
        type=str,
        help="PBC option for gmx trjconv",
        choices=["none", "mol", "res", "atom", "nojump", "cluster", "whole"],
    )

    parser.add_argument(
        "--prefix", default="prd", type=str, help="Prefix of trajectory files"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output trajectory path (default: <prefix>_all_rmmol.xtc)",
    )

    parser.add_argument("--skip", default=1, type=int, help="Number of frames to skip")

    parser.add_argument(
        "--no-resnr", action="store_true", help="Do not run gmx editconf -resnr 1"
    )
    parser.add_argument(
        "--keep-concatenated",
        action="store_true",
        help="Keep the intermediate concatenated trajectory before PBC processing",
    )

    parser.set_defaults(func=run)


def run(args):
    if args.num_of_step < 1:
        raise ValueError("--num_of_step must be positive")
    if args.skip < 1:
        raise ValueError("--skip must be positive")

    # ref: https://zenn.dev/kh01734/articles/012380a58949d1
    # make new represent topology
    topology = f"{args.prefix}1.tpr"
    concatenated = f"{args.prefix}_all.xtc"
    output = getattr(args, "output", None) or f"{args.prefix}_all_rmmol.xtc"
    cmd = [
        "gmx",
        "convert-tpr",
        "-s",
        topology,
        "-n",
        args.index,
        "-o",
        "rmmol_top.tpr",
    ]
    run_cmd(cmd, input=f"{args.keep_selection}\n", log="rmmol_top.tpr generated")

    if args.no_resnr:
        cmd = ["gmx", "editconf", "-f", "rmmol_top.tpr", "-o", "rmmol_top.gro"]
    else:
        cmd = [
            "gmx",
            "editconf",
            "-f",
            "rmmol_top.tpr",
            "-o",
            "rmmol_top.gro",
            "-resnr",
            "1",
        ]
    run_cmd(cmd, log="rmmol_top.gro generated")

    # trjcat -> trjconv
    c_cmd = "c\n" * args.num_of_step
    trj_files = [f"{args.prefix}{step}.xtc" for step in range(1, args.num_of_step + 1)]

    cmd = ["gmx", "trjcat", "-f", *trj_files, "-o", concatenated, "-settime"]
    run_cmd(cmd, input=c_cmd, log=f"{concatenated} generated")

    if args.pbc == "cluster":
        selection_input = (
            f"{args.centering_selection}\n"
            f"{args.centering_selection}\n"
            f"{args.keep_selection}\n"
        )
    else:
        selection_input = f"{args.centering_selection}\n{args.keep_selection}\n"
    cmd = [
        "gmx",
        "trjconv",
        "-f",
        concatenated,
        "-s",
        topology,
        "-n",
        args.index,
        "-o",
        output,
        "-skip",
        str(args.skip),
        "-pbc",
        args.pbc,
        "-center",
    ]
    run_cmd(cmd, input=selection_input, log=f"{output} generated")

    if not getattr(args, "keep_concatenated", False):
        Path(concatenated).unlink(missing_ok=True)
        remove_gmx_backups(Path(concatenated).parent)
        LOGGER.info(f"{concatenated} and GROMACS backup files removed")
