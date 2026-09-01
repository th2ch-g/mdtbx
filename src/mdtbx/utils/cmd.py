import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

from ..logger import generate_logger
from .proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """Register the project-environment command runner."""
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

    project_root = Path(__file__).resolve().parents[3]
    command = list(args.command)
    run_options = {}
    if command[0] == "packmol-memgen":
        compat_dir = Path(__file__).resolve().parents[1] / "_compat" / "packmol_memgen"
        env = os.environ.copy()
        current_pythonpath = env.get("PYTHONPATH")
        pythonpath_parts = [str(compat_dir)]
        if current_pythonpath:
            pythonpath_parts.append(current_pythonpath)
        env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)
        run_options["env"] = env
        amberhome = env.get("AMBERHOME")
        amber_packmol = Path(amberhome, "bin", "packmol") if amberhome else None
        if (
            "--packmol" not in command
            and amber_packmol is not None
            and amber_packmol.is_file()
            and os.access(amber_packmol, os.X_OK)
        ):
            command.extend(["--packmol", str(amber_packmol)])

    pixi_cmd = ["pixi", "run", "--manifest-path", str(project_root), *command]
    command_str = shlex.join(command)

    LOGGER.info(f"Executing command in pixi environment: {command_str}")
    try:
        run_cmd(pixi_cmd, **run_options)
        LOGGER.info("Command finished successfully.")
    except subprocess.CalledProcessError as error:
        LOGGER.error(f"Command failed with exit code {error.returncode}.")
        sys.exit(error.returncode)
