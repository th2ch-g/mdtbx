import argparse
import math

import numpy as np

from ..logger import generate_logger
from ..utils.common_args import add_output_arg, add_topology_arg
from ..utils.numpy_io import save_npy
from .distmap import load_representative_coordinates, pairwise_distances

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx contactmap -p structure.pdb -t trajectory.xtc -o contactmap.npy
    """
    parser = subparsers.add_parser(
        "contactmap",
        help="Extract residue contact matrix from representative atoms",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_topology_arg(parser)
    parser.add_argument(
        "-t",
        "--trajectory",
        type=str,
        help="Trajectory file (.xtc, .trr). If omitted, coordinates are read from topology",
    )
    parser.add_argument(
        "-s",
        "--selection",
        type=str,
        default="protein and name CA",
        help="MDtraj selection for representative atoms (typically one atom per residue)",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=6.0,
        help="Contact cutoff in angstrom",
    )
    add_output_arg(parser, default="contactmap.npy")

    parser.set_defaults(func=run)


def calculate_contact_map(distance_matrices: np.ndarray, cutoff: float) -> np.ndarray:
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be positive and finite")
    if distance_matrices.ndim != 3 or distance_matrices.shape[0] == 0:
        raise ValueError(
            "distance_matrices must have shape (n_frames, n_atoms, n_atoms)"
        )
    if distance_matrices.shape[1] != distance_matrices.shape[2]:
        raise ValueError("distance matrices must be square")
    if not np.all(np.isfinite(distance_matrices)):
        raise ValueError("distance matrices must contain only finite values")
    contact_map = (distance_matrices < cutoff).astype(float).mean(axis=0)
    np.fill_diagonal(contact_map, 0.0)
    return contact_map


def run(args):
    coordinates = load_representative_coordinates(
        args.topology,
        args.trajectory,
        args.selection,
    )
    if coordinates is None:
        return

    distance_matrices = pairwise_distances(coordinates)
    contact_map = calculate_contact_map(distance_matrices, args.cutoff)

    LOGGER.info(f"Contact map shape: {contact_map.shape}")
    output_path = save_npy(args.output, contact_map)
    LOGGER.info(f"Saved to {output_path}")
