import argparse

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "SUBCMD_NAME",
        help="HELP_TEXT",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # 必須引数
    parser.add_argument("input", type=str, help="Input file")
    # オプション引数
    parser.add_argument("--output", type=str, default="output.dat", help="Output file")

    parser.set_defaults(func=run)


def run(args):
    LOGGER.info(f"input: {args.input}")
    # TODO: 実装
