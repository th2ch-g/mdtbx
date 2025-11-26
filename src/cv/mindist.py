import argparse
import mdtraj as md
import numpy as np
import itertools

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx mindist --topology structure.pdb --trajectory trajectory.xtc --selection1 "resid 1 to 10" --selection2 "resid 11 to 20" -o cv.npy
    """
    parser = subparsers.add_parser(
        "mindist",
        help="Extract minimum distance between two sets of atoms",
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
        "-s1",
        "--selection1",
        type=str,
        required=True,
        help="Selection 1 (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-s2",
        "--selection2",
        type=str,
        required=True,
        help="Selection 2 (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="mindist.npy", help="Output file (.npy)"
    )


def run(args):
    trj = md.load(args.trajectory, top=args.topology)

    selection1_indices = trj.top.select(args.selection1)
    selection2_indices = trj.top.select(args.selection2)

    atom_pairs = list(itertools.product(selection1_indices, selection2_indices))

    distances = md.compute_distances(trj, atom_pairs)

    min_dist_per_frame = np.min(distances, axis=1)

    np.save(args.output, min_dist_per_frame)
    LOGGER.info(f"Saved to {args.output}")
