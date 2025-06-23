import argparse
import sys

from .logger import generate_logger

LOGGER = generate_logger(__name__)


def cli() -> None:
    # make parser
    parser = argparse.ArgumentParser(description=("ToolBox for MD simulation"))

    # subcommands
    subparsers = parser.add_subparsers()

    args = parser.parse_args()

    if len(sys.argv) == 1:
        LOGGER.error(f"use {sys.argv[0]} --help")
        sys.exit(1)

    LOGGER.info(f"{sys.argv[1]} called")

    LOGGER.info(f"{sys.argv[1]} finished")
