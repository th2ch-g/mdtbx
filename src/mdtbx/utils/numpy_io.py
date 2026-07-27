"""NumPy output helpers shared by analysis commands."""

from pathlib import Path
from typing import Any

import numpy as np

from .paths import ensure_output_parent


def _prepare_numpy_path(path: str | Path, suffix: str) -> Path:
    output_path = Path(path).expanduser()
    if output_path.suffix != suffix:
        output_path = Path(f"{output_path}{suffix}")
    return ensure_output_parent(output_path)


def save_npy(path: str | Path, array: Any) -> Path:
    """Save an array and create missing output directories."""
    output_path = _prepare_numpy_path(path, ".npy")
    np.save(output_path, array)
    return output_path


def save_npz(path: str | Path, **arrays: Any) -> Path:
    """Save named arrays in an uncompressed archive and create output directories."""
    output_path = _prepare_numpy_path(path, ".npz")
    np.savez(output_path, **arrays)
    return output_path
