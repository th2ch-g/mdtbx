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
from ..utils.mdtraj import select_atom_indices
from ..utils.numpy_io import save_npy
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)

_MDTRAJ_SELECTION_ALIASES = {
    "Protein": "protein",
    "Water": "water",
    "System": "all",
}


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


def average_by_residue(topology, atom_indices, values):
    """Average per-atom values for each selected residue in topology order."""
    if len(atom_indices) != len(values):
        raise ValueError("atom_indices and values must have the same length")

    grouped_values = {}
    for atom_index, value in zip(atom_indices, values, strict=True):
        residue_index = topology.atom(int(atom_index)).residue.index
        grouped_values.setdefault(residue_index, []).append(value)
    return np.asarray(
        [np.mean(residue_values) for residue_values in grouped_values.values()]
    )


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
        trj = md.load(args.trajectory, top=args.topology)
        selection = _MDTRAJ_SELECTION_ALIASES.get(args.selection, args.selection)
        atom_indices = select_atom_indices(trj.topology, selection)
        rmsf = md.rmsf(trj, trj, 0, atom_indices=atom_indices)
        if args.resolution == "residue":
            rmsf = average_by_residue(trj.topology, atom_indices, rmsf)
    output_path = save_npy(args.output, rmsf)
    LOGGER.info(f"Saved to {output_path}")
