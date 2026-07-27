"""Maximum-likelihood Markov state model estimation."""

import argparse

import numpy as np

from ..utils.numpy_io import save_npz
from ._kinetics import load_discrete_trajectories

_COUNT_MODES = ("sample", "sliding", "sliding-effective", "effective")


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "msm",
        help="Estimate a maximum-likelihood MSM from discrete trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input integer .npy files or cluster .npz archives",
    )
    parser.add_argument("--lagtime", type=int, required=True, help="Lag time in frames")
    parser.add_argument(
        "--count-mode",
        choices=_COUNT_MODES,
        default="effective",
        help="Transition counting mode",
    )
    parser.add_argument(
        "--nonreversible",
        action="store_true",
        help="Estimate a non-reversible transition matrix",
    )
    parser.add_argument("-o", "--output", default="msm.npz", help="Output archive")
    parser.set_defaults(func=run)


def run(args):
    from deeptime.markov import TransitionCountEstimator
    from deeptime.markov.msm import MaximumLikelihoodMSM

    if args.lagtime < 1:
        raise ValueError("--lagtime must be positive")
    trajectories, sources = load_discrete_trajectories(args.input)
    if any(len(trajectory) <= args.lagtime for trajectory in trajectories):
        raise ValueError("Every trajectory must be longer than --lagtime")
    count_model = TransitionCountEstimator(
        lagtime=args.lagtime,
        count_mode=args.count_mode,
    ).fit_fetch(trajectories)
    reversible = not args.nonreversible
    model = MaximumLikelihoodMSM(reversible=reversible, use_lcc=True).fit_fetch(
        count_model
    )
    state_symbols = np.asarray(model.state_symbols(), dtype=np.int64)
    output = save_npz(
        args.output,
        schema_version=np.int64(1),
        kind=np.str_("mdtbx.msm"),
        transition_matrix=np.asarray(model.transition_matrix),
        count_matrix=np.asarray(model.count_model.count_matrix),
        full_count_matrix=np.asarray(count_model.count_matrix_full),
        stationary_distribution=np.asarray(model.stationary_distribution),
        eigenvalues=np.asarray(model.eigenvalues),
        timescales=np.asarray(model.timescales()),
        state_symbols=state_symbols,
        state_histogram=np.asarray(model.count_model.state_histogram),
        lagtime=np.int64(args.lagtime),
        count_mode=np.str_(args.count_mode),
        reversible=np.bool_(reversible),
        source_paths=np.asarray(sources, dtype=np.str_),
    )
    return {
        "output": str(output),
        "trajectories": len(trajectories),
        "states": len(state_symbols),
    }
