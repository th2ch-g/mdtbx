import argparse
import mdtraj as md

from ..logger import generate_logger
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
    if args.gmx:
        INDEX_OPT = f"-n {args.index}"
        # cmd = f"echo {args.selection} System | gmx trjconv -f {args.file} -s {args.topology} -o tmp.xtc -pbc nojump -center"
        cmd = f"echo {args.selection} System | gmx trjconv -f {args.file} -s {args.topology} -o tmp.xtc -pbc {args.pbc} {INDEX_OPT} -center"
        run_cmd(cmd)
        cmd = f"echo {args.selection} System | gmx trjconv -f tmp.xtc -s {args.topology} -o {args.output} -fit rot+trans {INDEX_OPT}"
        run_cmd(cmd)
        cmd = "rm -f tmp.xtc"
        run_cmd(cmd, log="tmp.xtc removed")
    else:
        trj = md.load(args.file, top=args.topology)
        ref = md.load(args.topology)
        fit_trj = trj.superpose(
            ref,
            atom_indices=trj.top.select(args.selection),
            ref_atom_indices=ref.top.select(args.selection),
        )
        fit_trj.save(args.output)
    LOGGER.info(f"{args.output} generated")
