import argparse

import numpy as np

from ..logger import generate_logger
from .distmap import load_representative_coordinates

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

    parser.add_argument(
        "-p", "--topology", type=str, required=True, help="Topology file (.gro, .pdb)"
    )
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
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="contactmap.npy",
        help="Output file (.npy)",
    )

    parser.set_defaults(func=run)


def calculate_contact_map(distance_matrices: np.ndarray, cutoff: float) -> np.ndarray:
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

    diff = coordinates[:, :, None, :] - coordinates[:, None, :, :]
    distance_matrices = np.linalg.norm(diff, axis=-1)
    contact_map = calculate_contact_map(distance_matrices, args.cutoff)

    LOGGER.info(f"Contact map shape: {contact_map.shape}")
    np.save(args.output, contact_map)
    LOGGER.info(f"Saved to {args.output}")
