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

    parser.add_argument(
        "-s", "--skip", default=1, type=int, help="Number of frames to skip"
    )

    parser.add_argument(
        "-k", "--keep_selection", default="non-Water", type=str, help="Keep selection"
    )

    parser.add_argument(
        "-c", "--centering_selection", default="Protein", type=str, help="Centering selection"
    )

    parser.add_argument(
        "-n", "--num_of_step", required=True, type=int, help="Number of steps"
    )

    parser.add_argument(
        "-g", "--gmx", default="gmx", type=str, help="gmx command"
    )

    parser.add_argument(
        "-idx", "--index", default="index.ndx", type=str, help="Index file"
    )


def run(args):
    # make new represent topology
    cmd = f"echo {args.keep_selection} | {args.gmx} trjconv -f prd_1.gro -s prd_1.tpr -n {args.index} -o rmmol_top.gro"
    subprocess.run(cmd, shell=True)

    for step in range(1, args.num_of_step+1):
        cmd = f"echo {args.centering_selection} {args.keep_selection} | {args.gmx} trjconv -f prd_{step}.xtc -s prd_{step}.tpr -n {args.index} -o prd_{step}_skip{args.skip}_rmmol.xtc -pbc mol -center -skip {args.skip}"
        subprocess.run(cmd, shell=True)

    c_cmd = "c\n" * 100
    trj_files = [f"prd_{step}_skip{args.skip}_rmmol.xtc" for step in range(1, args.num_of_step+1)]
    trj_files = " ".join(trj_files)

    cmd = f"echo {c_cmd} | {args.gmx} trjcat -f {trj_files} -o prd_all_skip{args.skip}_rmmol.xtc -settime"

    subprocess.run(cmd, shell=True)

    cmd = f"rm -f {trj_files}"
    subprocess.run(cmd, shell=True)

