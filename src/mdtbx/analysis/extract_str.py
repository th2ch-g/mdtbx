import argparse
import sys

from ..logger import generate_logger
from ..utils.common_args import (
    add_gmx_args,
    add_output_arg,
    add_selection_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.gmx import gmx_index_args
from ..utils.mdtraj import select_atom_indices
from ..utils.paths import ensure_output_parent
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx extract_str --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10" -o target.pdb --time 123

    --time is a time point in ps, not a frame index, for both the mdtraj and
    the --gmx backend.
    """
    parser = subparsers.add_parser(
        "extract_str",
        help="Extract structure in specific time point",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_topology_arg(parser)
    add_trajectory_arg(parser)
    add_selection_arg(parser)
    parser.add_argument(
        "--time",
        type=float,
        required=True,
        help="Time point to extract [ps]",
    )
    add_output_arg(parser, default="target.pdb", help="Output struxture file")
    add_gmx_args(parser)

    parser.set_defaults(func=run)


def run(args):
    output_path = ensure_output_parent(args.output)

    if args.gmx:
        cmd = [
            "gmx",
            "trjconv",
            "-s",
            args.topology,
            "-f",
            args.trajectory,
            "-o",
            str(output_path),
            "-b",
            str(args.time),
            "-e",
            str(args.time),
            *gmx_index_args(args.index),
        ]
        run_cmd(cmd, input=f"{args.selection}\n")
    else:
        # mdtraj
        import mdtraj as md

        trj = md.load(args.trajectory, top=args.topology)
        # --time is a time point in ps for both backends: the gmx path feeds it
        # to `trjconv -b/-e`, so the mdtraj path resolves it against trj.time
        # instead of treating it as a frame index.
        times = trj.time
        if not (times.min() <= args.time <= times.max()):
            LOGGER.error(
                f"--time {args.time} ps is out of range "
                f"({times.min()} to {times.max()} ps)"
            )
            sys.exit(1)
        frame = int(abs(times - args.time).argmin())
        if times[frame] != args.time:
            LOGGER.warning(
                f"No frame at {args.time} ps; using the nearest frame "
                f"{frame} at {times[frame]} ps"
            )
        trj = trj[frame]
        atom_indices = select_atom_indices(trj.topology, args.selection)
        final_trj = trj.atom_slice(atom_indices)
        final_trj.save_pdb(str(output_path))
    LOGGER.info("Done")
