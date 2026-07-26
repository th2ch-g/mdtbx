import argparse

import mdtraj as md
import numpy as np

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.mdtraj import select_atom_indices
from ..utils.numpy_io import save_npy

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx rmsd --topology structure.pdb --trajectory trajectory.xtc --reference ref.pdb --selection_cal_trj "resid 1 to 10" --selection_cal_ref "resid 1 to 10" --selection_fit_trj "resid 1 to 10" --selection_fit_ref "resid 1 to 10" -o cv.npy
    """
    parser = subparsers.add_parser(
        "rmsd",
        help="Extract RMSD",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_topology_arg(parser)
    add_trajectory_arg(parser)
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=True,
        help="Reference structure file (.gro, .pdb)",
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
    add_output_arg(parser, default="rmsd.npy")

    parser.set_defaults(func=run)


def run(args):
    trj = md.load(args.trajectory, top=args.topology)
    ref = md.load(args.reference)[0]

    fit_trj = select_atom_indices(
        trj.topology, args.selection_fit_trj, label="trajectory fit selection"
    )
    fit_ref = select_atom_indices(
        ref.topology, args.selection_fit_ref, label="reference fit selection"
    )
    if len(fit_trj) != len(fit_ref):
        raise ValueError(
            "Fit selections contain different atom counts: "
            f"trajectory={len(fit_trj)}, reference={len(fit_ref)}"
        )

    cal_trj = select_atom_indices(
        trj.topology, args.selection_cal_trj, label="trajectory calculation selection"
    )
    cal_ref = select_atom_indices(
        ref.topology, args.selection_cal_ref, label="reference calculation selection"
    )
    if len(cal_trj) != len(cal_ref):
        raise ValueError(
            "Calculation selections contain different atom counts: "
            f"trajectory={len(cal_trj)}, reference={len(cal_ref)}"
        )

    trj.superpose(
        ref,
        0,
        atom_indices=fit_trj,
        ref_atom_indices=fit_ref,
    )
    # md.rmsd performs superposition automatically, so we don't use that functionality here
    rmsd = np.sqrt(
        3
        * np.mean(
            np.square(trj.xyz[:, cal_trj] - ref.xyz[0, cal_ref]),
            axis=(1, 2),
        )
    )
    output_path = save_npy(args.output, rmsd)
    LOGGER.info(f"Saved to {output_path}")
