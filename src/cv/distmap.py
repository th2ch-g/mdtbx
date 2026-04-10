import argparse

import mdtraj as md
import numpy as np

from ..logger import generate_logger

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
        "-o", "--output", type=str, default="distmap.npy", help="Output file (.npy)"
    )

    parser.set_defaults(func=run)


def load_representative_coordinates(
    topology_path: str,
    trajectory_path: str | None,
    selection: str,
) -> np.ndarray | None:
    if trajectory_path is not None:
        traj = md.load(trajectory_path, top=topology_path)
    else:
        traj = md.load(topology_path)

    atom_indices = traj.topology.select(selection)
    if len(atom_indices) == 0:
        LOGGER.error(f"No atoms selected by: {selection}")
        return None

    coordinates = traj.xyz[:, atom_indices, :] * 10.0
    LOGGER.info(
        "Selected %s representative atoms across %s frame(s)",
        len(atom_indices),
        traj.n_frames,
    )
    return coordinates


def calculate_distance_map(coordinates: np.ndarray) -> np.ndarray:
    diff = coordinates[:, :, None, :] - coordinates[:, None, :, :]
    return np.linalg.norm(diff, axis=-1).mean(axis=0)


def run(args):
    coordinates = load_representative_coordinates(
        args.topology,
        args.trajectory,
        args.selection,
    )
    if coordinates is None:
        return

    distance_map = calculate_distance_map(coordinates)
    LOGGER.info(f"Distance map shape: {distance_map.shape}")
    np.save(args.output, distance_map)
    LOGGER.info(f"Saved to {args.output}")
