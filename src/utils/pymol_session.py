"""Shared PyMOL session helper for the build subcommands.

Standardizes the `cmd.reinitialize()` + `cmd.load(structure, name)` idiom so
every PyMOL-driven subcommand starts from a clean state (avoiding leftover
object state when several run in the same process, e.g. under pytest).
"""

from contextlib import contextmanager


@contextmanager
def pymol_session(cmd, structure=None, name="target"):
    """Reinitialize PyMOL, optionally load ``structure`` into ``name``, yield cmd.

    ``cmd`` is the caller's ``pymol.cmd`` handle, passed in so each subcommand
    keeps its own (mockable) reference instead of this module binding pymol.
    """
    cmd.reinitialize()
    if structure is not None:
        cmd.load(structure, name)
    yield cmd
