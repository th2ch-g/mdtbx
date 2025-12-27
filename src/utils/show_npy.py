import argparse
import numpy as np

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx show_npy --npy npy
    """
    parser = subparsers.add_parser(
        "show_npy",
        help="show npy file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("npy", type=str, help="npy file")


def run(args):
    npy = np.load(args.npy)
    print(npy)
    print(npy.shape)
