import argparse
import mdtraj as md
import numpy as np

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx rmsd --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10" -o cv.npy
    """
    parser = subparsers.add_parser(
        "rmsd",
        help="Extract RMSD",
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
        "-r",
        "--reference",
        type=str,
        required=True,
        help="Reference trajectory file (.xtc, .trr)",
    )
    parser.add_argument(
        "-sct",
        "--selection_cal_trj",
        type=str,
        required=True,
        help="Selection for trajectory (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-scr",
        "--selection_cal_ref",
        type=str,
        required=True,
        help="Selection for reference (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-sft",
        "--selection_fit_trj",
        type=str,
        required=True,
        help="Selection for trajectory (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-sfr",
        "--selection_fit_ref",
        type=str,
        required=True,
        help="Selection for reference (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="rmsd.npy", help="Output file (.npy)"
    )


def run(args):
    trj = md.load(args.trajectory, top=args.topology)
    ref = md.load(args.reference)

    n_fit_trj = len(trj.top.select(args.selection_fit_trj))
    n_fit_ref = len(ref.top.select(args.selection_fit_ref))
    if n_fit_trj != n_fit_ref:
        LOGGER.error(
            f"Number of atoms in fit selection for trajectory ({n_fit_trj}) and reference ({n_fit_ref}) are different"
        )

    n_cal_trj = len(trj.top.select(args.selection_cal_trj))
    n_cal_ref = len(ref.top.select(args.selection_cal_ref))
    if n_cal_trj != n_cal_ref:
        LOGGER.error(
            f"Number of atoms in cal selection for trajectory ({n_cal_trj}) and reference ({n_cal_ref}) are different"
        )

    trj.superpose(
        ref,
        0,
        atom_indices=trj.top.select(args.selection_fit_trj),
        ref_atom_indices=ref.top.select(args.selection_fit_ref),
    )
    # md.rmsd performs superposition automatically, so we don't use that typeality here
    rmsd = np.sqrt(
        3
        * np.mean(
            np.square(
                trj.xyz[:, trj.top.select(args.selection_cal_trj)]
                - ref.xyz[:, ref.top.select(args.selection_cal_ref)]
            ),
            axis=(1, 2),
        )
    )
    np.save(args.output, rmsd)
    LOGGER.info(f"Saved to {args.output}")
