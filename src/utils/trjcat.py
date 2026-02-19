import argparse
import subprocess

from ..logger import generate_logger

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

    # parser.add_argument(
    #     "-s", "--skip", default=1, type=int, help="Number of frames to skip"
    # )

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

    parser.add_argument("--skip", default=1, type=int, help="Number of frames to skip")

    parser.add_argument(
        "--no-resnr", action="store_true", help="Do not run gmx editconf -resnr 1"
    )


def run(args):
    # ref: https://zenn.dev/kh01734/articles/012380a58949d1
    # make new represent topology
    cmd = f"echo {args.keep_selection} | gmx convert-tpr -s {args.prefix}1.tpr -n {args.index} -o rmmol_top.tpr"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("rmmol_top.tpr generated")

    if not args.no_resnr:
        cmd = "gmx editconf -f rmmol_top.tpr -o rmmol_top.gro"
    else:
        cmd = "gmx editconf -f rmmol_top.tpr -o rmmol_top.gro -resnr 1"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("rmmol_top.gro generated")

    # make gro with trjconv
    # cmd = f"echo {args.centering_selection} {args.keep_selection} | gmx trjconv -f {args.prefix}1.gro -s {args.prefix}1.tpr -n {args.index} -o rmmol_top.gro -center"
    # subprocess.run(cmd, shell=True, check=True)
    # LOGGER.info("rmmol_top.gro generated")
    #
    # if not args.no_editconf:
    #     cmd = "gmx editconf -f rmmol_top.gro -o rmmol_top.gro -resnr 1"
    #     subprocess.run(cmd, shell=True, check=True)
    #     LOGGER.info("gmx editconf -f rmmol_top.gro -o rmmol_top.gro -resnr 1 runned")

    # trjcat -> trjconv
    c_cmd = "c\n" * args.num_of_step
    trj_files = [f"{args.prefix}{step}.xtc" for step in range(1, args.num_of_step + 1)]
    trj_files = " ".join(trj_files)

    cmd = (
        f"echo '{c_cmd}' | gmx trjcat -f {trj_files} -o {args.prefix}_all.xtc -settime"
    )
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.prefix}_all.xtc generated")

    if args.pbc == "cluster":
        cmd = f"echo {args.centering_selection} {args.centering_selection} {args.keep_selection} | gmx trjconv -f {args.prefix}_all.xtc -s {args.prefix}1.tpr -n {args.index} -o {args.prefix}_all_rmmol.xtc -skip {args.skip} -pbc {args.pbc} -center"
    else:
        cmd = f"echo {args.centering_selection} {args.keep_selection} | gmx trjconv -f {args.prefix}_all.xtc -s {args.prefix}1.tpr -n {args.index} -o {args.prefix}_all_rmmol.xtc -skip {args.skip} -pbc {args.pbc} -center"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.prefix}_all_rmmol.xtc generated")

    cmd = f"rm -f {args.prefix}_all.xtc \#*"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.prefix}_all.xtc and backup files removed")
