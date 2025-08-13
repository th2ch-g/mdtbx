import argparse
import mdtraj as md
import numpy as np

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx xyz --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10" -o cv.npy
    """
    parser = subparsers.add_parser(
        "xyz",
        help="Extract XYZ coordinates",
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


def run(args):
    trj = md.load(args.trajectory, top=args.topology)
    xyz = trj.atom_slice(trj.topology.select(args.selection)).xyz
    np.save(args.output, xyz)
    LOGGER.info(f"Saved to {args.output}")
