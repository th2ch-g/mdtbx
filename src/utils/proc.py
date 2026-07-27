"""Thin subprocess wrapper shared across subcommands."""

import subprocess
import sys
from collections.abc import Sequence

from ..agent.context import JSON_MODE
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def run_cmd(
    cmd: str | Sequence[str],
    *,
    log: str | None = None,
    check: bool = True,
    **kwargs,
):
    """Run an external command from a string or preferred argv sequence."""
    shell = isinstance(cmd, str)
    if isinstance(kwargs.get("input"), str) and not any(
        option in kwargs for option in ("text", "universal_newlines", "encoding")
    ):
        kwargs["text"] = True
    if JSON_MODE.get() and not kwargs.get("capture_output"):
        kwargs.setdefault("stdout", sys.stderr)
        kwargs.setdefault("stderr", sys.stderr)
    result = subprocess.run(cmd, shell=shell, check=check, **kwargs)
    if log:
        LOGGER.info(log)
    return result
