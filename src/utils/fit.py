import argparse
import mdtraj as md

from ..logger import generate_logger

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
        help="Centering selection (MDtraj atom selection language)",
    )


def run(args):
    trj = md.load(args.file, top=args.topology)
    ref = md.load(args.topology)
    fit_trj = trj.superpose(
        ref,
        atom_indices=trj.top.select(args.selection),
        ref_atom_indices=ref.top.select(args.selection),
    )
    fit_trj.save(args.output)
    LOGGER.info(f"{args.output} generated")
