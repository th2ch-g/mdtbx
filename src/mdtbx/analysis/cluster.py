"""K-means clustering for molecular feature trajectories."""

import argparse

import numpy as np

from ..utils.numpy_io import save_npz
from ._kinetics import concatenate, load_feature_trajectories


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "cluster",
        help="Cluster feature trajectories with deeptime K-means",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input 2D feature .npy files or tICA .npz archives",
    )
    parser.add_argument("--n-clusters", type=int, required=True, help="Cluster count")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--max-iter", type=int, default=500, help="Maximum iterations")
    parser.add_argument("--n-jobs", type=int, default=1, help="Parallel workers")
    parser.add_argument("-o", "--output", default="clusters.npz", help="Output archive")
    parser.set_defaults(func=run)


def run(args):
    from deeptime.clustering import KMeans

    if args.n_clusters < 1:
        raise ValueError("--n-clusters must be positive")
    if args.seed < 0:
        raise ValueError("--seed must be non-negative")
    if args.max_iter < 1:
        raise ValueError("--max-iter must be positive")
    if args.n_jobs < 1:
        raise ValueError("--n-jobs must be positive")
    trajectories, sources = load_feature_trajectories(args.input)
    features, _ = concatenate(trajectories)
    if args.n_clusters > len(features):
        raise ValueError("--n-clusters must not exceed the number of frames")
    estimator = KMeans(
        n_clusters=args.n_clusters,
        max_iter=args.max_iter,
        fixed_seed=args.seed,
        n_jobs=args.n_jobs,
    )
    model = estimator.fit_fetch(features)
    dtrajs = [
        np.asarray(model.transform(item), dtype=np.int64) for item in trajectories
    ]
    discrete, discrete_lengths = concatenate(dtrajs)
    output = save_npz(
        args.output,
        schema_version=np.int64(1),
        kind=np.str_("mdtbx.cluster"),
        dtrajs=discrete,
        lengths=discrete_lengths,
        centers=np.asarray(model.cluster_centers),
        inertias=np.asarray(model.inertias),
        inertia=np.float64(model.inertia),
        converged=np.bool_(model.converged),
        n_clusters=np.int64(args.n_clusters),
        seed=np.int64(args.seed),
        max_iter=np.int64(args.max_iter),
        n_jobs=np.int64(args.n_jobs),
        source_paths=np.asarray(sources, dtype=np.str_),
    )
    return {
        "output": str(output),
        "trajectories": len(trajectories),
        "frames": len(features),
        "converged": bool(model.converged),
    }
