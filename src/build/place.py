import argparse
import sys

import numpy as np

from ..logger import generate_logger

LOGGER = generate_logger(__name__)

# Pairwise atom distance below this is treated as a serious clash.
# Smaller than typical vdW contact (~3.0 A) but larger than covalent bonds (~1.5 A).
DEFAULT_CLASH_CUTOFF = 2.0


def add_subcmd(subparsers):
    """
    mdtbx place -1 chainA.pdb -2 chainB.pdb -d 30.0 --seed 42 -o placed.pdb
    """
    parser = subparsers.add_parser(
        "place",
        help="Place two single-chain PDBs at a given distance with random rotation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-1",
        "--pdb1",
        required=True,
        type=str,
        help="First PDB file (single chain, kept fixed at origin)",
    )

    parser.add_argument(
        "-2",
        "--pdb2",
        required=True,
        type=str,
        help="Second PDB file (single chain, randomly rotated then translated)",
    )

    parser.add_argument(
        "-d",
        "--distance",
        required=True,
        type=float,
        help="Distance between centers of mass of the two chains (Angstrom)",
    )

    parser.add_argument(
        "--seed",
        default=42,
        type=int,
        help="Random seed for rotation of the second chain",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="placed.pdb",
        type=str,
        help="Output PDB file (combined)",
    )

    parser.add_argument(
        "--clash_cutoff",
        default=DEFAULT_CLASH_CUTOFF,
        type=float,
        help="Min allowed inter-chain atom distance (Angstrom)",
    )

    parser.add_argument(
        "--ignore",
        action="store_true",
        help="Warn instead of aborting when a serious clash is detected",
    )

    parser.set_defaults(func=run)


def _spawn_rngs(seed: int, n: int) -> list[np.random.Generator]:
    """
    Derive n independent RNG streams from one user-facing seed.
    Using SeedSequence.spawn avoids correlation between consumers
    that would otherwise share the same seed.
    """
    children = np.random.SeedSequence(seed).spawn(n)
    return [np.random.default_rng(c) for c in children]


def random_rotation_matrix(seed: int) -> np.ndarray:
    """
    Build a uniformly distributed 3x3 rotation matrix from a random seed.

    Strategy: draw three Euler angles (alpha, beta, gamma) and combine
    Rz(alpha) @ Ry(beta) @ Rz(gamma). beta is sampled via arccos(1 - 2u)
    so that the rotation is uniform on SO(3) rather than biased toward poles.
    """
    rng = _spawn_rngs(seed, 2)[0]
    alpha = rng.uniform(0.0, 2.0 * np.pi)
    beta = np.arccos(1.0 - 2.0 * rng.uniform(0.0, 1.0))
    gamma = rng.uniform(0.0, 2.0 * np.pi)

    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)

    rz1 = np.array([[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]])
    ry = np.array([[cb, 0.0, sb], [0.0, 1.0, 0.0], [-sb, 0.0, cb]])
    rz2 = np.array([[cg, -sg, 0.0], [sg, cg, 0.0], [0.0, 0.0, 1.0]])

    return rz1 @ ry @ rz2


def random_unit_vector(seed: int) -> np.ndarray:
    """
    Sample a uniformly distributed unit vector on the 2-sphere from a seed.
    Muller's method: draw 3D standard normal, then normalize.
    """
    rng = _spawn_rngs(seed, 2)[1]
    v = rng.standard_normal(3)
    return v / np.linalg.norm(v)


def min_interchain_distance(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Smallest pairwise distance between two atom coordinate sets.
    coords1: (N1, 3), coords2: (N2, 3) -> scalar in Angstrom.
    """
    diff = coords1[:, None, :] - coords2[None, :, :]
    dists = np.linalg.norm(diff, axis=-1)
    return float(dists.min())


def run(args):
    from pymol import cmd as pymol_cmd

    pymol_cmd.reinitialize()
    pymol_cmd.load(args.pdb1, "mol1")
    pymol_cmd.load(args.pdb2, "mol2")

    com1 = np.array(pymol_cmd.centerofmass("mol1"))
    com2 = np.array(pymol_cmd.centerofmass("mol2"))
    LOGGER.info(f"mol1 COM (input): {com1}")
    LOGGER.info(f"mol2 COM (input): {com2}")

    # 1. center mol1 at origin
    pymol_cmd.translate(list(-com1), "mol1", camera=0)

    # 2. center mol2 at origin (so rotation is around its own COM)
    pymol_cmd.translate(list(-com2), "mol2", camera=0)

    # 3. rotate mol2 with a seeded random rotation (4x4 TTT matrix)
    R = random_rotation_matrix(args.seed)
    ttt = [
        R[0, 0],
        R[0, 1],
        R[0, 2],
        0.0,
        R[1, 0],
        R[1, 1],
        R[1, 2],
        0.0,
        R[2, 0],
        R[2, 1],
        R[2, 2],
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
    ]
    pymol_cmd.transform_selection("mol2", ttt)
    LOGGER.info(f"mol2 rotated with seed {args.seed}")

    # 4. translate mol2 by distance d along a seeded random direction
    direction = random_unit_vector(args.seed)
    translation = direction * float(args.distance)
    pymol_cmd.translate(list(translation), "mol2", camera=0)
    LOGGER.info(
        f"mol2 translated by {translation} (direction={direction}, |t|={args.distance})"
    )

    com1_new = np.array(pymol_cmd.centerofmass("mol1"))
    com2_new = np.array(pymol_cmd.centerofmass("mol2"))
    LOGGER.info(f"mol1 COM (placed): {com1_new}")
    LOGGER.info(f"mol2 COM (placed): {com2_new}")
    LOGGER.info(f"COM-COM distance: {np.linalg.norm(com2_new - com1_new):.3f} A")

    # Clash check on the placed coordinates (heavy + hydrogen, all atoms).
    coords1 = np.asarray(pymol_cmd.get_coords("mol1"))
    coords2 = np.asarray(pymol_cmd.get_coords("mol2"))
    min_dist = min_interchain_distance(coords1, coords2)
    LOGGER.info(f"min inter-chain atom distance: {min_dist:.3f} A")

    if min_dist < args.clash_cutoff:
        msg = (
            f"serious clash detected: min inter-chain distance "
            f"{min_dist:.3f} A < cutoff {args.clash_cutoff} A"
        )
        if args.ignore:
            LOGGER.warning(msg + " (continuing because --ignore is set)")
        else:
            LOGGER.error(msg + " (use --ignore to write the PDB anyway)")
            sys.exit(1)

    pymol_cmd.create("placed", "mol1 or mol2")
    pymol_cmd.save(args.output, "placed")
    LOGGER.info(f"{args.output} generated")
