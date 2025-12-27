import argparse
import numpy as np
import subprocess

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx densmap --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10"  -o desmap.npy
    """
    parser = subparsers.add_parser(
        "densmap",
        help="Extract densmap",
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
        help="Selection (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="comdist.npy", help="Output file (.npy)"
    )

    parser.add_argument(
        "--gmx", action="store_true", help="Use Gromacs instead of MDtraj"
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
            INDEX_FILE = f"-n {args.index}"
        else:
            INDEX_FILE = ""
        cmd = f"gmx densmap -f {args.trajectory} -s {args.topology} {INDEX_FILE} -od densmap.dat"
        subprocess.run(cmd, input=f"{args.selection}\n", shell=True, check=True)
        densmap = np.loadtxt("densmap.dat")
        # X = a[0, 1:]
        # Y = a[1:, 0]
        # dens = a[1:, 1]
        # c = ax.pcolormesh(X, Y, dens)
    else:
        # trj = md.load(args.trajectory, top=args.topology)
        raise NotImplementedError

    np.save(args.output, densmap)
    LOGGER.info(f"Saved to {args.output}")
