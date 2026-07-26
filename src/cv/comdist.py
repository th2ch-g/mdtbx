import argparse
import numpy as np

import mdtraj as md

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.gmx import gmx_index_args, gmx_tempfile
from ..utils.mdtraj import select_atom_indices
from ..utils.numpy_io import save_npy
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx comdist --topology structure.pdb --trajectory trajectory.xtc --selection1 "resid 1 to 10" --selection2 "resid 11 to 20" -o cv.npy
    """
    parser = subparsers.add_parser(
        "comdist",
        help="Extract center of mass distance between two sets of atoms",
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
    add_output_arg(parser, default="comdist.npy")

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
        with gmx_tempfile(".xvg") as displacement_path:
            cmd = [
                "gmx",
                "distance",
                "-f",
                args.trajectory,
                "-s",
                args.topology,
                *gmx_index_args(args.index),
                "-oxyz",
                displacement_path,
                "-xvg",
                "none",
                "-pbc",
                "no",
                "-select",
                f"com of group {args.selection1} plus com of group {args.selection2}",
            ]
            run_cmd(cmd, log=f"Saved to {displacement_path}")
            inter_com_xyz = np.loadtxt(displacement_path, ndmin=2)
            comdist = np.linalg.norm(inter_com_xyz[:, [1, 2, 3]], axis=1)
    else:
        trj = md.load(args.trajectory, top=args.topology)
        selection1 = select_atom_indices(
            trj.topology, args.selection1, label="selection1"
        )
        selection2 = select_atom_indices(
            trj.topology, args.selection2, label="selection2"
        )
        com1 = md.compute_center_of_mass(trj.atom_slice(selection1))
        com2 = md.compute_center_of_mass(trj.atom_slice(selection2))
        comdist = np.linalg.norm(com1 - com2, axis=1)
    LOGGER.info(f"COM distance: {np.mean(comdist):.3f} +/- {np.std(comdist):.3f}")
    output_path = save_npy(args.output, comdist)
    LOGGER.info(f"Saved to {output_path}")
