import argparse
from pathlib import Path
from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx build_membrane
    """
    parser = subparsers.add_parser(
        "build_membrane",
        help="Build membrane system",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    pass


def run(args):
    pass
