import argparse
import mdtraj as md

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_selection_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.mdtraj import select_atom_indices
from ..utils.numpy_io import save_npy

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
    atom_indices = select_atom_indices(trj.topology, args.selection)
    xyz = trj.atom_slice(atom_indices).xyz
    output_path = save_npy(args.output, xyz)
    LOGGER.info(f"Saved to {output_path}")
