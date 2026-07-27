"""Helpers shared by the GROMACS-backed subcommands.

Centralizes three patterns that were copy-pasted across cv/ and trajectory/:
the optional `-n <index>` flag, the 0-based -> 1-based atom-id conversion, and
the "write a temp file, read it, delete it" dance around gmx output.
"""

import os
import tempfile
from contextlib import contextmanager
from pathlib import Path


def gmx_index_flag(index):
    """Return the gmx ``-n <index>`` option, or ``""`` when no index file is set."""
    return f"-n {index}" if index is not None else ""


def gmx_index_args(index):
    """Return the gmx index option as argv elements."""
    return ["-n", str(index)] if index is not None else []


def to_gmx_index(indices):
    """Convert 0-based atom indices to a sorted list of 1-based GROMACS atom ids."""
    return sorted(i + 1 for i in indices)


def remove_gmx_backups(directory="."):
    """Remove GROMACS backup files from one directory."""
    for path in Path(directory).glob("#*#"):
        if path.is_file() or path.is_symlink():
            path.unlink()


@contextmanager
def gmx_tempfile(suffix=".xvg"):
    """Yield a unique temp filename in the current directory, removed on exit.

    The placeholder is removed before yielding so gmx creates the output
    without generating a backup. The file is always cleaned up, even if
    reading it raises. Unique names also avoid collisions between concurrent
    runs in the same directory.
    """
    fd, path = tempfile.mkstemp(suffix=suffix, dir=".")
    os.close(fd)
    Path(path).unlink()
    try:
        yield path
    finally:
        Path(path).unlink(missing_ok=True)
