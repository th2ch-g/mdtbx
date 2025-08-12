from pathlib import Path
import argparse

from ..logger import generate_logger

LOGGER = generate_logger(__name__)

def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "rmfile",
        help="Remove unnecessary files related to MD simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("path", type=str, help="Path to the directory", default=".")


def run(args):
    for suffix in ["\#*", "*cpt", "mdout.mdp"]:
        for p in Path(args.path).glob(suffix):
            LOGGER.info(f"rm -f {p}")
            p.unlink()
