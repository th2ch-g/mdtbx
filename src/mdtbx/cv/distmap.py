import argparse

import mdtraj as md
import numpy as np

from ..logger import generate_logger
from ..utils.common_args import add_output_arg, add_topology_arg
from ..utils.mdtraj import select_atom_indices
from ..utils.numpy_io import save_npy

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx distmap -p structure.pdb -t trajectory.xtc -o distmap.npy
    """
    parser = subparsers.add_parser(
        "distmap",
        help="Extract residue distance matrix from representative atoms",
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
    add_output_arg(parser, default="distmap.npy")

    parser.set_defaults(func=run)


def load_representative_coordinates(
    topology_path: str,
    trajectory_path: str | None,
    selection: str,
) -> np.ndarray:
    if trajectory_path is not None:
        traj = md.load(trajectory_path, top=topology_path)
    else:
        traj = md.load(topology_path)

    atom_indices = select_atom_indices(traj.topology, selection)

    coordinates = traj.xyz[:, atom_indices, :] * 10.0
    LOGGER.info(
        "Selected %s representative atoms across %s frame(s)",
        len(atom_indices),
        traj.n_frames,
    )
    return coordinates


def _validate_coordinates(coordinates: np.ndarray) -> None:
    if coordinates.ndim != 3 or coordinates.shape[0] == 0 or coordinates.shape[1] == 0:
        raise ValueError("coordinates must have shape (n_frames, n_atoms, 3)")
    if coordinates.shape[2] != 3:
        raise ValueError("coordinates must have shape (n_frames, n_atoms, 3)")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("coordinates must contain only finite values")


def pairwise_distances(coordinates: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(coordinates, dtype=float)
    _validate_coordinates(coordinates)
    diff = coordinates[:, :, None, :] - coordinates[:, None, :, :]
    return np.linalg.norm(diff, axis=-1)


def reduce_pairwise_distances(
    coordinates: np.ndarray,
    transform=None,
    chunk_frames: int = 100,
) -> np.ndarray:
    """Average transform(pairwise distances) over frames in bounded memory.

    The broadcast in pairwise_distances materializes an
    (n_frames, n_atoms, n_atoms, 3) tensor, which for long trajectories does
    not fit in memory; both distmap and contactmap only need per-frame
    reductions, so accumulate chunk sums instead.
    """
    if chunk_frames < 1:
        raise ValueError("chunk_frames must be positive")
    coordinates = np.asarray(coordinates, dtype=float)
    _validate_coordinates(coordinates)
    n_frames = coordinates.shape[0]
    total = None
    for start in range(0, n_frames, chunk_frames):
        block = pairwise_distances(coordinates[start : start + chunk_frames])
        if transform is not None:
            block = transform(block)
        block_sum = block.sum(axis=0, dtype=float)
        total = block_sum if total is None else total + block_sum
    return total / n_frames


def calculate_distance_map(coordinates: np.ndarray) -> np.ndarray:
    return reduce_pairwise_distances(coordinates)


def run(args):
    coordinates = load_representative_coordinates(
        args.topology,
        args.trajectory,
        args.selection,
    )
    distance_map = calculate_distance_map(coordinates)
    LOGGER.info(f"Distance map shape: {distance_map.shape}")
    output_path = save_npy(args.output, distance_map)
    LOGGER.info(f"Saved to {output_path}")
