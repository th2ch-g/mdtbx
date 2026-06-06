import argparse
import importlib
import pkgutil
import sys
from importlib import metadata

from . import analysis, build, cv, trajectory, utils
from .logger import generate_logger

LOGGER = generate_logger(__name__)

# Category packages scanned for subcommand modules. Each subcommand module
# exposes add_subcmd(subparsers); library modules (atom_selection_parser,
# parse_top, proc, tleap, gmx, common_args, pymol_session, ...) do not and are
# skipped automatically.
_SUBCOMMAND_PACKAGES = (utils, build, trajectory, analysis, cv)


def get_version():
    try:
        return metadata.version("mdtbx")
    except metadata.PackageNotFoundError:
        return "unknown"


def _register_subcommands(subparsers) -> None:
    """Import every module under the category packages exposing add_subcmd."""
    for pkg in _SUBCOMMAND_PACKAGES:
        for info in sorted(pkgutil.iter_modules(pkg.__path__), key=lambda m: m.name):
            module = importlib.import_module(f"{pkg.__name__}.{info.name}")
            add_subcmd = getattr(module, "add_subcmd", None)
            if add_subcmd is not None:
                add_subcmd(subparsers)


def cli() -> None:
    # make parser
    parser = argparse.ArgumentParser(description=("ToolBox for MD simulation"))

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {get_version()}",
        help="Print version",
    )

    # subcommands (auto-discovered from the category packages)
    subparsers = parser.add_subparsers()
    _register_subcommands(subparsers)

    args = parser.parse_args()

    if not hasattr(args, "func"):
        LOGGER.error(f"use {sys.argv[0]} --help")
        sys.exit(1)

    LOGGER.info(f"{sys.argv[1]} called")
    args.func(args)
    LOGGER.info(f"{sys.argv[1]} finished")
