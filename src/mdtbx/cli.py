"""Command-line entry point."""

import argparse
import importlib
import pkgutil
import sys
import warnings
from importlib import metadata
from pathlib import Path

from . import analysis, build, cv, trajectory, utils
from .logger import generate_logger

LOGGER = generate_logger(__name__)

# Category packages scanned for subcommand modules. Each subcommand module
# exposes add_subcmd(subparsers); library modules without it are skipped.
_SUBCOMMAND_PACKAGES = (utils, build, trajectory, analysis, cv)


def get_version():
    try:
        return metadata.version("mdtbx")
    except metadata.PackageNotFoundError:
        return "unknown"


def _module_for_command(name: str) -> str | None:
    """Resolve conventional command modules without importing unrelated code."""
    for package in _SUBCOMMAND_PACKAGES:
        if (Path(next(iter(package.__path__))) / f"{name}.py").is_file():
            return f"{package.__name__}.{name}"
    return None


def _register_module(subparsers, module_name: str) -> None:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"pkg_resources is deprecated as an API\.",
            category=UserWarning,
            module=r"mdtraj\.geometry\.order",
        )
        module = importlib.import_module(module_name)
    add_subcmd = getattr(module, "add_subcmd", None)
    if add_subcmd is None:
        return
    add_subcmd(subparsers)


def _register_subcommands(subparsers, commands: set[str] | None = None) -> None:
    """Import every module under the category packages exposing add_subcmd."""
    if commands is not None:
        modules = {_module_for_command(name) for name in commands}
        if None not in modules:
            for module_name in sorted(modules):
                _register_module(subparsers, module_name)
            return
    for package in _SUBCOMMAND_PACKAGES:
        for info in sorted(
            pkgutil.iter_modules(package.__path__), key=lambda item: item.name
        ):
            _register_module(subparsers, f"{package.__name__}.{info.name}")


def create_parser(commands: set[str] | None = None) -> argparse.ArgumentParser:
    """Create the top-level argument parser with every registered subcommand."""
    parser = argparse.ArgumentParser(description="ToolBox for MD simulation")
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {get_version()}",
        help="Print version",
    )
    subparsers = parser.add_subparsers(dest="_command")
    _register_subcommands(subparsers, commands)
    return parser


def cli() -> None:
    selected = None
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        selected = {sys.argv[1]}
    parser = create_parser(selected)
    args = parser.parse_args()
    if not hasattr(args, "func"):
        LOGGER.error(f"use {sys.argv[0]} --help")
        raise SystemExit(1)

    command = args._command or Path(sys.argv[0]).name
    LOGGER.info(f"{command} called")
    args.func(args)
    LOGGER.info(f"{command} finished")
