import argparse
from pathlib import Path
from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx build_solution
    """
    parser = subparsers.add_parser(
        "build_solution",
        help="Build solution system",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    pass


def run(args):
    pass
