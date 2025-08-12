import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen-sample-builder
    """
    parser = subparsers.add_parser(
        "gen-sample-builder",
        help="Generate sample builder scripts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default=".",
        help="Path to generate sample scripts",
    )

    parser.add_argument(
        "--type",
        type=str,
        default="solution",
        help="Type of sample scripts",
        choices=["solution", "membrane"],
    )


def run(args):
    for scrips in Path(__file__).parent.glob(f"{args.type}/sample_*"):
        cmd = f"cp {scrips} {args.path}"
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"chmod +x {args.path}/{scrips.name}"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info(f"{args.path}/{scrips.name} generated")
    for scrips in Path(__file__).parent.glob(f"{args.type}/*mdp"):
        cmd = f"cp {scrips} {args.path}"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info(f"{args.path}/{scrips.name} generated")
