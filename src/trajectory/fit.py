import argparse
import mdtraj as md

from ..logger import generate_logger
from ..utils.gmx import gmx_index_args, gmx_tempfile
from ..utils.mdtraj import select_atom_indices
from ..utils.paths import ensure_output_parent
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx fit -f trj -p top -o trj -s selection
    """
    parser = subparsers.add_parser(
        "fit",
        help="Fit trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-f", "--file", required=True, type=str, help="Trajectory file")

    parser.add_argument(
        "-p",
        "--topology",
        required=True,
        type=str,
        help="Topology file and Reference structure",
    )

    parser.add_argument(
        "-o", "--output", default="out.xtc", type=str, help="Output trajectory file"
    )

    parser.add_argument(
        "-s",
        "--selection",
        default="protein",
        type=str,
        help="Selection (MDtraj atom selection language)",
    )
    parser.add_argument("--gmx", action="store_true", help="Use gmx instead of MDtraj")

    parser.add_argument(
        "--pbc",
        default="mol",
        type=str,
        help="PBC option for gmx trjconv",
        choices=["none", "mol", "res", "atom", "nojump", "cluster", "whole"],
    )

    parser.add_argument(
        "-idx", "--index", default="index.ndx", type=str, help="Index file"
    )

    parser.set_defaults(func=run)


def run(args):
    output_path = ensure_output_parent(args.output)

    if args.gmx:
        with gmx_tempfile(".xtc") as centered_path:
            cmd = [
                "gmx",
                "trjconv",
                "-f",
                args.file,
                "-s",
                args.topology,
                "-o",
                centered_path,
                "-pbc",
                args.pbc,
                *gmx_index_args(args.index),
                "-center",
            ]
            run_cmd(cmd, input=f"{args.selection}\nSystem\n")
            cmd = [
                "gmx",
                "trjconv",
                "-f",
                centered_path,
                "-s",
                args.topology,
                "-o",
                str(output_path),
                "-fit",
                "rot+trans",
                *gmx_index_args(args.index),
            ]
            run_cmd(cmd, input=f"{args.selection}\nSystem\n")
    else:
        trj = md.load(args.file, top=args.topology)
        ref = md.load(args.topology)[0]
        fit_trj_indices = select_atom_indices(
            trj.topology, args.selection, label="trajectory fit selection"
        )
        fit_ref_indices = select_atom_indices(
            ref.topology, args.selection, label="reference fit selection"
        )
        if len(fit_trj_indices) != len(fit_ref_indices):
            raise ValueError(
                "Fit selections contain different atom counts: "
                f"trajectory={len(fit_trj_indices)}, reference={len(fit_ref_indices)}"
            )
        fit_trj = trj.superpose(
            ref,
            atom_indices=fit_trj_indices,
            ref_atom_indices=fit_ref_indices,
        )
        fit_trj.save(str(output_path))
    LOGGER.info(f"{output_path} generated")
