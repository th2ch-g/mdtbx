import argparse
import sys

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_selection_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.gmx import gmx_index_flag
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx extract_str --topology structure.pdb --trajectory trajectory.xtc --selection "resid 1 to 10" -o target.pdb --time 123
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
        type=int,
        required=True,
        help="Time to extract structure",
    )
    add_output_arg(parser, default="target.pdb", help="Output struxture file")
    parser.add_argument(
        "--gmx",
        action="store_true",
        help="Use gmx",
    )
    parser.add_argument(
        "-idx",
        "--index",
        type=str,
        help="Index file (.ndx)",
    )

    parser.set_defaults(func=run)


def run(args):
    if args.gmx:
        INDEX_OPTION = gmx_index_flag(args.index)
        cmd = f"gmx trjconv -s {args.topology} -f {args.trajectory} -o {args.output} -b {args.time} -e {args.time} {INDEX_OPTION}"
        run_cmd(cmd, input=f"{args.selection}\n")
    else:
        # mdtraj
        import mdtraj as md

        trj = md.load(args.trajectory, top=args.topology)
        if not (1 <= args.time <= trj.n_frames):
            LOGGER.error(f"--time {args.time} is out of range (1 to {trj.n_frames})")
            sys.exit(1)
        trj = trj[args.time - 1]
        atom_indices = trj.top.select(args.selection)
        final_trj = trj.atom_slice(atom_indices)
        final_trj.save_pdb(args.output)
    LOGGER.info("Done")
