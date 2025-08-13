import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_sample_mdp
    """
    parser = subparsers.add_parser(
        "gen_sample_mdp",
        help="Generate sample mdp files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default=".",
        help="Path to generate sample mdp",
    )

    parser.add_argument(
        "--type",
        type=str,
        default="solution",
        help="Type of sample mdp",
        choices=["solution", "membrane"],
    )


def run(args):
    for scrips in Path(__file__).parent.glob(f"{args.type}/*mdp"):
        cmd = f"cp {scrips} {args.path}"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info(f"{args.path}/{scrips.name} generated")
