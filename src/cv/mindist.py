import argparse
import itertools

import mdtraj as md
import numpy as np

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.mdtraj import select_atom_indices
from ..utils.numpy_io import save_npy

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

    add_topology_arg(parser)
    add_trajectory_arg(parser)
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
    add_output_arg(parser, default="mindist.npy")

    parser.set_defaults(func=run)


def run(args):
    trj = md.load(args.trajectory, top=args.topology)

    selection1_indices = select_atom_indices(
        trj.topology, args.selection1, label="selection1"
    )
    selection2_indices = select_atom_indices(
        trj.topology, args.selection2, label="selection2"
    )

    atom_pairs = list(itertools.product(selection1_indices, selection2_indices))

    distances = md.compute_distances(trj, atom_pairs)

    min_dist_per_frame = np.min(distances, axis=1)

    output_path = save_npy(args.output, min_dist_per_frame)
    LOGGER.info(f"Saved to {output_path}")
