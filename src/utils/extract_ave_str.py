import argparse
import subprocess

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx comdist --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10" -o ave.pdb
    """
    parser = subparsers.add_parser(
        "extract_ave_str",
        help="Extract average structures in trajectory",
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
        "-o", "--output", type=str, default="ave.pdb", help="Output struxture file"
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
        cmd = f"gmx covar -s {args.topology} -f {args.trajectory} -o eigenval.xvg -v eigenvec.trr -av {args.output} {INDEX_OPTION}"
        subprocess.run(cmd, input=f"{args.selection}\n", shell=True, check=True)
    else:
        # mdtraj
        import mdtraj as md

        trj = md.load(args.trajectory, top=args.topology)
        ave_crds = trj.xyz.mean(axis=0)
        atom_indices = trj.top.select(args.selection)
        final_trj = ave_crds[atom_indices]
        final_trj.save_pdb(args.output)
    LOGGER.info("Done")
