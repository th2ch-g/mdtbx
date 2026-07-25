import argparse
import numpy as np

import mdtraj as md

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.gmx import gmx_tempfile
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
        with gmx_tempfile(".xvg") as rmsf_path:
            cmd = [
                "gmx",
                "rmsf",
                "-f",
                args.trajectory,
                "-s",
                args.topology,
                "-o",
                rmsf_path,
                "-xvg",
                "none",
            ]
            if args.resolution == "residue":
                cmd.append("-res")
            run_cmd(cmd, input=f"{args.selection}\n", log=f"Saved to {rmsf_path}")
            rmsf = np.loadtxt(rmsf_path, ndmin=2)[:, 1]
    else:
        # mdtraj
        trj = md.load(args.trajectory, top=args.topology)
        atom_indices = trj.top.select(args.selection)
        rmsf = md.rmsf(trj, trj, 0, atom_indices=atom_indices)
    np.save(args.output, rmsf)
    LOGGER.info(f"Saved to {args.output}")
