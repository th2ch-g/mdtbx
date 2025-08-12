import argparse
import sys

from .logger import generate_logger

from .utils import rmfile
from .utils import addace
from .utils import addnme
from .utils import trjcat
from .utils import centering_gro
from .utils import convert
from .utils import find_bond
from .utils import gen_am1bcc
from .utils import gen_resp

LOGGER = generate_logger(__name__)


def cli() -> None:
    # make parser
    parser = argparse.ArgumentParser(description=("ToolBox for MD simulation"))

    # subcommands
    subparsers = parser.add_subparsers()

    rmfile.add_subcmd(subparsers)
    addace.add_subcmd(subparsers)
    addnme.add_subcmd(subparsers)
    trjcat.add_subcmd(subparsers)
    centering_gro.add_subcmd(subparsers)
    convert.add_subcmd(subparsers)
    find_bond.add_subcmd(subparsers)
    gen_am1bcc.add_subcmd(subparsers)
    gen_resp.add_subcmd(subparsers)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        LOGGER.error(f"use {sys.argv[0]} --help")
        sys.exit(1)
    LOGGER.info(f"{sys.argv[1]} called")

    if sys.argv[1] == "rmfile":
        rmfile.run(args)

    if sys.argv[1] == "addace":
        addace.run(args)

    if sys.argv[1] == "addnme":
        addnme.run(args)

    if sys.argv[1] == "trjcat":
        trjcat.run(args)

    if sys.argv[1] == "centering_gro":
        centering_gro.run(args)

    if sys.argv[1] == "convert":
        convert.run(args)

    if sys.argv[1] == "find_bond":
        find_bond.run(args)

    if sys.argv[1] == "gen_am1bcc":
        gen_am1bcc.run(args)

    if sys.argv[1] == "gen_resp":
        gen_resp.run(args)

    LOGGER.info(f"{sys.argv[1]} finished")
