import argparse

from ..logger import generate_logger
from ..utils.common_args import (
    add_output_arg,
    add_selection_arg,
    add_topology_arg,
    add_trajectory_arg,
)
from ..utils.gmx import gmx_index_args, gmx_tempfile
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx extract_ave_str
    """
    parser = subparsers.add_parser(
        "extract_ave_str",
        help="Extract average structures in trajectory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_topology_arg(parser)
    add_trajectory_arg(parser)
    add_selection_arg(parser)
    add_output_arg(parser, default="ave.pdb", help="Output struxture file")
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
        with (
            gmx_tempfile(".xvg") as eigenvalue_path,
            gmx_tempfile(".trr") as eigenvector_path,
        ):
            cmd = [
                "gmx",
                "covar",
                "-s",
                args.topology,
                "-f",
                args.trajectory,
                "-o",
                eigenvalue_path,
                "-v",
                eigenvector_path,
                "-av",
                args.output,
                *gmx_index_args(args.index),
            ]
            run_cmd(cmd, input=f"{args.selection}\n{args.selection}\n")
    else:
        # mdtraj
        import mdtraj as md

        trj = md.load(args.trajectory, top=args.topology)
        atom_indices = trj.top.select(args.selection)
        sel_trj = trj.atom_slice(atom_indices)
        ave_xyz = sel_trj.xyz.mean(axis=0, keepdims=True)  # (1, n_sel_atoms, 3)
        avg_trj = md.Trajectory(ave_xyz, sel_trj.topology)
        avg_trj.save_pdb(args.output)
    LOGGER.info("Done")
