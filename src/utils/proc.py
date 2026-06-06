"""Thin subprocess wrapper shared across subcommands.

Replaces the repeated `subprocess.run(cmd, shell=True, check=True)` + LOGGER
boilerplate with a single helper that accepts either an argv list (preferred,
no shell) or a command string (shell), and optionally logs a success message.
"""

import subprocess

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def run_cmd(cmd, *, log=None, check=True, **kwargs):
    """Run an external command.

    Parameters
    ----------
    cmd : list[str] | str
        argv list (run without a shell) or a command string (run via the shell,
        for cases that need redirection or pipes).
    log : str | None
        Message logged at INFO level on success.
    check : bool
        Raise CalledProcessError on a non-zero exit (default True).
    **kwargs
        Forwarded to subprocess.run (e.g. cwd, env, stdout).
    """
    shell = isinstance(cmd, str)
    result = subprocess.run(cmd, shell=shell, check=check, **kwargs)
    if log:
        LOGGER.info(log)
    return result
