import argparse
import subprocess

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx extract_str --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10" -o target.pdb --time 123
    """
    parser = subparsers.add_parser(
        "extract_str",
        help="Extract structure in specific time point",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p", "--topology", type=str, required=True, help="Topology file (.gro, .pdb)"
    )
    parser.add_argument(
        "-t",
        "--trajectory",
        type=str,
        required=True,
        help="Trajectory file (.xtc, .trr)",
    )
    parser.add_argument(
        "-s",
        "--selection",
        type=str,
        required=True,
        help="Selection (MDtraj atom selection language)",
    )
    parser.add_argument(
        "--time",
        type=int,
        required=True,
        help="Time to extract structure",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="target.pdb", help="Output struxture file"
    )
    parser.add_argument(
        "--gmx",
        action="store_true",
        help="Use gmx",
    )
    parser.add_argument(
        "-idx",
        "--index",
        type=str,
        help="Index file (.ndx)",
    )


def run(args):
    if args.gmx:
        if args.index is not None:
            INDEX_OPTION = f"-n {args.index}"
        else:
            INDEX_OPTION = ""
        cmd = f"gmx trjconv -s {args.topology} -f {args.trajectory} -o {args.output} -b {args.time} -e {args.time} {INDEX_OPTION}"
        subprocess.run(cmd, input=f"{args.selection}\n", shell=True, check=True)
    else:
        # mdtraj
        import mdtraj as md

        trj = md.load(args.trajectory, top=args.topology)
        trj = trj[args.time - 1]
        atom_indices = trj.top.select(args.selection)
        final_trj = trj.atom_slice(atom_indices)
        final_trj.save_pdb(args.output)
    LOGGER.info("Done")
