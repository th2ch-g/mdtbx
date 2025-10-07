import argparse
import mdtraj as md
import numpy as np
import subprocess

from ..config import *  # NOQA
from ..logger import generate_logger

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
    parser.add_argument(
        "-o", "--output", type=str, default="comvec.npy", help="Output file (.npy)"
    )
    parser.add_argument("--gmx", action="store_true", help="Use gmx instead of MDtraj")
    parser.add_argument(
        "-idx",
        "--index",
        type=str,
        help="Index file (.ndx)",
    )


def run(args):
    if args.gmx:
        if args.index is not None:
            INDEX_FILE = f"-n {args.index}"
        else:
            INDEX_FILE = ""
        cmd = f"gmx distance -f {args.trajectory} -s {args.topology} {INDEX_FILE} -oxyz tmp_interCOM_xyz.xvg -xvg none -pbc no -select 'com of group {args.selection1} plus com of group {args.selection2}'"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info("Saved to tmp_interCOM_xyz.xvg")
        inter_com_xyz = np.loadtxt("tmp_interCOM_xyz.xvg")
        comvec = inter_com_xyz[:, [1, 2, 3]]
        cmd = "rm -f tmp_interCOM_xyz.xvg"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info("Removed tmp_interCOM_xyz.xvg")
    else:
        trj = md.load(args.trajectory, top=args.topology)
        com1 = md.compute_center_of_mass(
            trj.atom_slice(trj.topology.select(args.selection1))
        )
        com2 = md.compute_center_of_mass(
            trj.atom_slice(trj.topology.select(args.selection2))
        )
        comvec = np.array([com1[i] - com2[i] for i in range(len(com1))])
    np.save(args.output, comvec)
    LOGGER.info(f"Saved to {args.output}")
