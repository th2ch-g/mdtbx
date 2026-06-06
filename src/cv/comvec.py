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
    mdtbx comvec --topology structure.pdb --trajectory trajectory.xtc --selection1 "resid 1 to 10" --selection2 "resid 11 to 20" -o cv.npy
    """
    parser = subparsers.add_parser(
        "comvec",
        help="Extract center of mass vector between two sets of atoms",
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
    add_output_arg(parser, default="comvec.npy")
    parser.add_argument("--gmx", action="store_true", help="Use gmx instead of MDtraj")
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
            # gmx -oxyz yields pos2 - pos1 (com2 - com1); negate to match the
            # mdtraj branch which returns com1 - com2, so output direction is
            # backend-independent.
            comvec = -inter_com_xyz[:, [1, 2, 3]]
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
        comvec = com1 - com2
    np.save(args.output, comvec)
    LOGGER.info(f"Saved to {args.output}")
