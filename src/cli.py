import argparse
import sys

from .logger import generate_logger

from .utils import rmfile
from .utils import addace

LOGGER = generate_logger(__name__)


def cli() -> None:
    # make parser
    parser = argparse.ArgumentParser(description=("ToolBox for MD simulation"))

    # subcommands
    subparsers = parser.add_subparsers()

    rmfile.add_subcmd(subparsers)
    addace.add_subcmd(subparsers)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        LOGGER.error(f"use {sys.argv[0]} --help")
        sys.exit(1)
    LOGGER.info(f"{sys.argv[1]} called")

    if sys.argv[1] == "rmfile":
        rmfile.run(args)

    if sys.argv[1] == "addace":
        addace.run(args)

    LOGGER.info(f"{sys.argv[1]} finished")
