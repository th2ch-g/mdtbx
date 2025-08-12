import argparse

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx mutate
    """
    parser = subparsers.add_parser(
        "mutate",
        help="Mutate residues in PDB file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )


def run(args):
    pass
