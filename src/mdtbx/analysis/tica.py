"""Time-lagged independent component analysis."""

import argparse

import numpy as np

from ..utils.numpy_io import save_npz
from ._kinetics import concatenate, load_feature_trajectories


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "tica",
        help="Fit tICA to one or more feature trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input 2D feature .npy files or tICA .npz archives",
    )
    parser.add_argument("--lagtime", type=int, required=True, help="Lag time in frames")
    parser.add_argument(
        "-n", "--n-components", type=int, default=2, help="Output dimensions"
    )
    parser.add_argument("--epsilon", type=float, default=1e-6, help="Rank cutoff")
    parser.add_argument(
        "--scaling",
        choices=("kinetic_map", "commute_map", "none"),
        default="kinetic_map",
        help="tICA component scaling",
    )
    parser.add_argument("-o", "--output", default="tica.npz", help="Output archive")
    parser.set_defaults(func=run)


def run(args):
    from deeptime.decomposition import TICA

    if args.lagtime < 1:
        raise ValueError("--lagtime must be positive")
    if args.n_components < 1:
        raise ValueError("--n-components must be positive")
    if not np.isfinite(args.epsilon) or args.epsilon <= 0:
        raise ValueError("--epsilon must be a positive finite number")
    trajectories, sources = load_feature_trajectories(args.input)
    if any(len(trajectory) <= args.lagtime for trajectory in trajectories):
        raise ValueError("Every trajectory must be longer than --lagtime")
    feature_count = trajectories[0].shape[1]
    if args.n_components > feature_count:
        raise ValueError("--n-components must not exceed the feature dimension")
    scaling = None if args.scaling == "none" else args.scaling
    model = TICA(
        lagtime=args.lagtime,
        dim=args.n_components,
        epsilon=args.epsilon,
        scaling=scaling,
    ).fit_fetch(trajectories)
    transformed = [
        np.asarray(model.transform(item), dtype=np.float64) for item in trajectories
    ]
    if any(item.shape[1] != args.n_components for item in transformed):
        raise ValueError("tICA rank is smaller than --n-components")
    if not all(np.isfinite(item).all() for item in transformed):
        raise ValueError("tICA produced non-finite projections")
    scores, lengths = concatenate(transformed)
    components = np.asarray(model.instantaneous_coefficients)[:, : args.n_components]
    timelagged = np.asarray(model.timelagged_coefficients)[:, : args.n_components]
    output = save_npz(
        args.output,
        schema_version=np.int64(1),
        kind=np.str_("mdtbx.tica"),
        scores=scores,
        lengths=lengths,
        components=components,
        timelagged_components=timelagged,
        singular_values=np.asarray(model.singular_values)[: args.n_components],
        timescales=np.asarray(model.timescales())[: args.n_components],
        mean=np.asarray(model.mean_0),
        lagtime=np.int64(args.lagtime),
        n_components=np.int64(args.n_components),
        epsilon=np.float64(args.epsilon),
        scaling=np.str_(args.scaling),
        source_paths=np.asarray(sources, dtype=np.str_),
    )
    return {
        "output": str(output),
        "trajectories": len(trajectories),
        "frames": len(scores),
    }
