import argparse
import mdtraj as md
import numpy as np

from ..config import *  # NOQA
from ..logger import generate_logger

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

    parser.add_argument("-p", "--topology", type=str, required=True, help="Topology file (.gro, .pdb)")
    parser.add_argument("-t", "--trajectory", type=str, required=True, help="Trajectory file (.xtc, .trr)")
    parser.add_argument("-s1", "--selection1", type=str, required=True, help="Selection 1 (MDtraj Atom selection language)")
    parser.add_argument("-s2", "--selection2", type=str, required=True, help="Selection 2 (MDtraj Atom selection language)")
    parser.add_argument("-o", "--output", type=str, default="comdist.npy", help="Output file (.npy)")

def run(args):
    trj = md.load(args.trajectory, top=args.topology)
    com1 = md.compute_center_of_mass(trj.atom_slice(trj.topology.select(selection1)))
    com2 = md.compute_center_of_mass(trj.atom_slice(trj.topology.select(selection2)))
    comdist = np.array([np.linalg.norm(com1[i] - com2[i]) for i in range(len(com1))])
    LOGGER.info(f"COM distance: {np.mean(comdist):.3f} +/- {np.std(comdist):.3f}")
    np.save(args.output, comdist)
    LOGGER.info(f"Saved to {args.output}")
