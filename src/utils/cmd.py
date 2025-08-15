import argparse
import subprocess
import sys

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx cmd
    """
    parser = subparsers.add_parser(
        "cmd",
        help="Run command",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("cmd", type=str, help="Command to run")


def run(args):
    cmd = " ".join(sys.argv[2:])
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{cmd} runned")
