import argparse
import shlex
import subprocess
import sys
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx cmd
    """
    parser = subparsers.add_parser(
        "cmd",
        help="Run a command within the mdtbx project environment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help="The command and its arguments to execute.",
    )

    parser.set_defaults(func=run)


def run(args):
    if not args.command:
        LOGGER.error("No command provided to run.")
        return

    project_root = Path(__file__).parent.parent.parent.resolve()
    pixi_cmd = ["pixi", "run", "--manifest-path", str(project_root), *args.command]
    command_str = shlex.join(args.command)

    LOGGER.info(f"Executing command in pixi environment: {command_str}")
    try:
        subprocess.run(pixi_cmd, check=True)
        LOGGER.info("Command finished successfully.")
    except subprocess.CalledProcessError as e:
        LOGGER.error(f"Command failed with exit code {e.returncode}.")
        sys.exit(e.returncode)
