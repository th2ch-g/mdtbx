import argparse
import mdtraj as md
import numpy as np

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_selection_arg,
    add_topology_arg,
    add_trajectory_arg,
)

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

    add_topology_arg(parser)
    add_trajectory_arg(parser)
    add_selection_arg(parser, help="Selection (MDtraj Atom selection language)")
    add_output_arg(parser, default="xyz.npy")

    parser.set_defaults(func=run)


def run(args):
    trj = md.load(args.trajectory, top=args.topology)
    xyz = trj.atom_slice(trj.topology.select(args.selection)).xyz
    np.save(args.output, xyz)
    LOGGER.info(f"Saved to {args.output}")
