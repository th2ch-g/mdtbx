import argparse
import mdtraj as md
import numpy as np

from pathlib import Path

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.proc import run_cmd

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

    add_topology_arg(parser)
    add_trajectory_arg(parser)
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
    add_output_arg(parser, default="rmsf.npy")

    parser.set_defaults(func=run)


def run(args):
    if args.gmx:
        # gmx
        if args.resolution == "atom":
            RESOLUTION = ""
        else:
            RESOLUTION = "-res"
        cmd = f"echo {args.selection} | gmx rmsf -f {args.trajectory} -s {args.topology} -o tmp_rmsf.xvg -xvg none {RESOLUTION}"
        run_cmd(cmd, log="Saved to tmp_rmsf.xvg")
        try:
            rmsf = np.loadtxt("tmp_rmsf.xvg")[:, 1]
        finally:
            Path("tmp_rmsf.xvg").unlink(missing_ok=True)
            LOGGER.info("Removed tmp_rmsf.xvg")
    else:
        # mdtraj
        trj = md.load(args.trajectory, top=args.topology)
        atom_indices = trj.top.select(args.selection)
        rmsf = md.rmsf(trj, trj, 0, atom_indices=atom_indices)
    np.save(args.output, rmsf)
    LOGGER.info(f"Saved to {args.output}")
