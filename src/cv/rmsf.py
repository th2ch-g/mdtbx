import argparse
import mdtraj as md
import numpy as np
import subprocess

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx rmsf --topology structure.pdb --trajectory trajectory.xtc --selection "Protein" -o cv.npy --gmx --resolution atom
    """
    parser = subparsers.add_parser(
        "rmsf",
        help="Extract RMSF",
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
        "--selection",
        type=str,
        default="Protein",
        help="Selection",
    )
    parser.add_argument(
        "--resolution",
        type=str,
        choices=["atom", "residue"],
        default="residue",
        help="Resolution for RMSF",
    )
    parser.add_argument(
        "--gmx",
        action="store_true",
        help="Use gmx rmsf",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="rmsf.npy", help="Output file (.npy)"
    )


def run(args):
    if args.gmx:
        # gmx
        if args.resolution == "atom":
            RESOLUTION = ""
        else:
            RESOLUTION = "-res"
        cmd = f"echo {args.selection} | gmx rmsf -f {args.trajectory} -s {args.topology} -o {args.output} -xvg none {RESOLUTION}"
        subprocess.run(cmd, shell=True, check=True)
        rmsf = np.loadtxt(args.output)
    else:
        # mdtraj
        trj = md.load(args.trajectory, top=args.topology)
        atom_indices = trj.top.select(args.selection)
        rmsf = md.rmsf(trj, atom_indices)
    np.save(args.output, rmsf)
    LOGGER.info(f"Saved to {args.output}")
