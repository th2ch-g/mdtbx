from pathlib import Path
import argparse

from ..logger import generate_logger
from .proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "rmfile",
        help="Remove unnecessary files related to MD simulation (e.g. .cpt, mdout.mdp, backup files)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--path", type=str, help="Path to the directory", default=".")

    parser.set_defaults(func=run)


def run(args):
    for suffix in ["#*#", "*.cpt", "mdout.mdp"]:
        # for p in Path(args.path).glob(suffix):
        for p in Path(args.path).rglob(suffix):
            cmd = f"rm -f '{p}'"
            run_cmd(cmd, log=f"{p} removed")
