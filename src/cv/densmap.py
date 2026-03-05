import argparse
import numpy as np
import subprocess

from ..config import *  # NOQA
from ..logger import generate_logger

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
        "-o", "--output", type=str, default="densmap.npy", help="Output file (.npy)"
    )
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
        import mdtraj as md

        trj = md.load(args.trajectory, top=args.topology)
        atom_indices = trj.topology.select(args.selection)
        if len(atom_indices) == 0:
            LOGGER.error(f"No atoms selected by: {args.selection}")
            return

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

    np.save(args.output, densmap)
    LOGGER.info(f"Saved to {args.output}")
