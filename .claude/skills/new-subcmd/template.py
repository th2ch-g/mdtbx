import argparse

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "SUBCMD_NAME",
        help="HELP_TEXT",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Required arguments
    parser.add_argument("input", type=str, help="Input file")
    # Optional arguments
    parser.add_argument("--output", type=str, default="output.dat", help="Output file")

    parser.set_defaults(func=run)


def run(args):
    LOGGER.info(f"input: {args.input}")
    # TODO: implement
