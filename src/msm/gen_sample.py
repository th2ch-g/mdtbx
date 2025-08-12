import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_sample_msm
    """
    parser = subparsers.add_parser(
        "gen_sample_msm",
        help="Generate sample msm scripts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default=".",
        help="Path to generate sample scripts",
    )


def run(args):
    for scrips in Path(__file__).parent.glob("sample_*"):
        cmd = f"cp {scrips} {args.path}"
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"chmod +x {args.path}/{scrips.name}"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info(f"{args.path}/{scrips.name} generated")
