"""Command-line entry point for humans and autonomous agents."""

import argparse
import contextlib
import importlib
import io
import json
import pkgutil
import sys
import time
import traceback
import warnings
from importlib import metadata
from pathlib import Path

from . import agent, analysis, build, cv, trajectory, utils
from .agent.context import JSON_MODE
from .agent.model import json_value
from .agent.registry import AGENT_COMMANDS, MUTATING_PATH_COMMANDS, OUTPUT_DESTS
from .agent.runtime import direct_plan
from .logger import generate_logger

LOGGER = generate_logger(__name__)
_MAX_CAPTURE_BYTES = 1024 * 1024

# Category packages scanned for subcommand modules. Each subcommand module
# exposes add_subcmd(subparsers); library modules without it are skipped.
_SUBCOMMAND_PACKAGES = (utils, build, trajectory, analysis, cv, agent)
_AGENT_MODULES = {
    "agent_cancel": "mdtbx.agent.cancel",
    "agent_collect": "mdtbx.agent.collect",
    "agent_plan": "mdtbx.agent.plan",
    "agent_probe": "mdtbx.agent.probe",
    "agent_profile_save": "mdtbx.agent.profile_save",
    "agent_run": "mdtbx.agent.execute",
    "agent_schema": "mdtbx.agent.schema",
    "agent_status": "mdtbx.agent.status",
}


def get_version():
    try:
        return metadata.version("mdtbx")
    except metadata.PackageNotFoundError:
        return "unknown"


def _add_agent_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--json",
        dest="_agent_json",
        action="store_true",
        help="Emit one machine-readable result object",
    )
    parser.add_argument(
        "--dry-run",
        dest="_agent_dry_run",
        action="store_true",
        help="Validate arguments and emit an execution plan without running",
    )
    parser.add_argument(
        "--cluster-profile",
        dest="_agent_cluster_profile",
        help="External cluster profile JSON used by --dry-run",
    )


def _module_for_command(name: str) -> str | None:
    """Resolve conventional command modules without importing unrelated code."""
    if name in _AGENT_MODULES:
        return _AGENT_MODULES[name]
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
    before = set(subparsers.choices)
    add_subcmd(subparsers)
    for name in set(subparsers.choices) - before:
        if name not in AGENT_COMMANDS:
            _add_agent_flags(subparsers.choices[name])


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


def _capture_value(buffer: io.StringIO) -> str:
    value = buffer.getvalue()
    encoded = value.encode("utf-8")
    if len(encoded) <= _MAX_CAPTURE_BYTES:
        return value
    suffix = "\n[output truncated by mdtbx]\n"
    return encoded[:_MAX_CAPTURE_BYTES].decode("utf-8", errors="ignore") + suffix


def _artifact_paths(args) -> list[str]:
    paths = []
    names = set(OUTPUT_DESTS)
    # Mirror registry.descriptor(): these commands mutate their --path
    # directory in place, so it belongs in the artifacts list too.
    if getattr(args, "_command", None) in MUTATING_PATH_COMMANDS:
        names.add("path")
    for name in names:
        value = getattr(args, name, None)
        values = value if isinstance(value, list) else [value]
        for item in values:
            if isinstance(item, str) and item:
                paths.append(str(Path(item).expanduser()))
    return sorted(set(paths))


def _envelope(
    *,
    command: str,
    started: float,
    data=None,
    error: BaseException | None = None,
    captured: str = "",
    artifacts: list[str] | None = None,
) -> dict:
    payload = {
        "schema_version": 1,
        "command": command,
        "ok": error is None,
        "exit_code": 0,
        "duration_seconds": round(time.monotonic() - started, 6),
        "data": json_value(data),
        "artifacts": artifacts or [],
    }
    if captured:
        payload["stdout"] = captured
    if error is not None:
        if isinstance(error, SystemExit):
            exit_code = error.code
        else:
            exit_code = getattr(error, "returncode", None)
        if not isinstance(exit_code, int):
            exit_code = 130 if isinstance(error, KeyboardInterrupt) else 1
        payload["exit_code"] = exit_code
        payload["error"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
    return payload


def _run_agent_mode(args) -> int:
    started = time.monotonic()
    command = args._command
    capture = io.StringIO()
    token = JSON_MODE.set(True)
    result = None
    error = None
    try:
        with contextlib.redirect_stdout(capture):
            if getattr(args, "_agent_dry_run", False):
                result = direct_plan(command, args)
            else:
                result = args.func(args)
    except (Exception, KeyboardInterrupt, SystemExit) as caught:
        error = caught
        LOGGER.debug(traceback.format_exc())
    finally:
        JSON_MODE.reset(token)
    payload = _envelope(
        command=command,
        started=started,
        data=result,
        error=error,
        captured=_capture_value(capture),
        artifacts=_artifact_paths(args),
    )
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    return int(payload["exit_code"])


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
    agent_mode = command in AGENT_COMMANDS or getattr(args, "_agent_json", False)
    agent_mode = agent_mode or getattr(args, "_agent_dry_run", False)
    if agent_mode:
        raise SystemExit(_run_agent_mode(args))

    LOGGER.info(f"{command} called")
    args.func(args)
    LOGGER.info(f"{command} finished")
