import argparse
import numpy as np

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_selection_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.gmx import gmx_index_args, gmx_tempfile
from ..utils.mdtraj import select_atom_indices
from ..utils.numpy_io import save_npy
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)

_AXIS_MAP = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}


def add_subcmd(subparsers):
    """
    mdtbx densmap --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10"  -o densmap.npy
    """
    parser = subparsers.add_parser(
        "densmap",
        help="Extract 2D density map",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_topology_arg(parser)
    add_trajectory_arg(parser)
    add_selection_arg(
        parser,
        help=(
            "Selection. MDtraj Atom selection language for the MDtraj path; "
            "with --gmx this must be a Gromacs index group name or number"
        ),
    )
    add_output_arg(parser, default="densmap.npy")
    parser.add_argument(
        "--bins", type=int, default=100, help="Number of bins along each axis"
    )
    parser.add_argument(
        "--axis",
        type=str,
        default="xy",
        choices=["xy", "xz", "yz"],
        help="Projection plane (MDtraj path only)",
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

    parser.set_defaults(func=run)


def run(args):
    if args.bins <= 0:
        raise ValueError("--bins must be positive")

    if args.gmx:
        with gmx_tempfile(".dat") as density_path:
            cmd = [
                "gmx",
                "densmap",
                "-f",
                args.trajectory,
                "-s",
                args.topology,
                *gmx_index_args(args.index),
                "-od",
                density_path,
            ]
            run_cmd(cmd, input=f"{args.selection}\n")
            # Layout: first row is X, first column is Y, and the body is density.
            a = np.loadtxt(density_path)
        # Save in the same [counts, edges0, edges1] object-array form as the MDtraj path
        densmap = np.empty(3, dtype=object)
        densmap[0] = a[1:, 1:]
        densmap[1] = a[0, 1:]
        densmap[2] = a[1:, 0]
    else:
        import mdtraj as md

        trj = md.load(args.trajectory, top=args.topology)
        atom_indices = select_atom_indices(trj.topology, args.selection)

        xyz = trj.xyz[:, atom_indices, :]  # (n_frames, n_atoms, 3) [nm]
        ax0, ax1 = _AXIS_MAP[args.axis]
        pos0 = xyz[:, :, ax0].ravel()
        pos1 = xyz[:, :, ax1].ravel()

        counts, edges0, edges1 = np.histogram2d(pos0, pos1, bins=args.bins)
        # Save as dict-like structured array: [counts, edges0, edges1]
        # Use object array to preserve shapes
        densmap = np.empty(3, dtype=object)
        densmap[0] = counts
        densmap[1] = edges0
        densmap[2] = edges1
        LOGGER.info(f"Density map shape: {counts.shape}, axis: {args.axis}")

    output_path = save_npy(args.output, densmap)
    LOGGER.info(f"Saved to {output_path}")
