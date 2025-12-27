import argparse
import mdtraj as md

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx show_mdtraj topology.pdb
    """
    parser = subparsers.add_parser(
        "show_mdtraj",
        help="show mdtraj's dataframe",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "topology", type=str, help="Topology file (.gro, .pdb)"
    )


def run(args):
    top = md.load(args.topology)
    atoms, bonds = top.to_dataframe()
    print(atoms)
