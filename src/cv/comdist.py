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
from ..utils.gmx import gmx_index_flag
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
        INDEX_FILE = gmx_index_flag(args.index)
        cmd = f"gmx distance -f {args.trajectory} -s {args.topology} {INDEX_FILE} -oxyz tmp_interCOM_xyz.xvg -xvg none -pbc no -select 'com of group {args.selection1} plus com of group {args.selection2}'"
        run_cmd(cmd, log="Saved to tmp_interCOM_xyz.xvg")
        try:
            inter_com_xyz = np.loadtxt("tmp_interCOM_xyz.xvg")
            comdist = np.linalg.norm(inter_com_xyz[:, [1, 2, 3]], axis=1)
        finally:
            Path("tmp_interCOM_xyz.xvg").unlink(missing_ok=True)
            LOGGER.info("Removed tmp_interCOM_xyz.xvg")
    else:
        trj = md.load(args.trajectory, top=args.topology)
        com1 = md.compute_center_of_mass(
            trj.atom_slice(trj.topology.select(args.selection1))
        )
        com2 = md.compute_center_of_mass(
            trj.atom_slice(trj.topology.select(args.selection2))
        )
        comdist = np.linalg.norm(com1 - com2, axis=1)
    LOGGER.info(f"COM distance: {np.mean(comdist):.3f} +/- {np.std(comdist):.3f}")
    np.save(args.output, comdist)
    LOGGER.info(f"Saved to {args.output}")
