"""Shared validation helpers for MDTraj-backed commands."""

from typing import Any

import numpy as np


def select_atom_indices(
    topology: Any,
    selection: str,
    *,
    label: str = "selection",
) -> np.ndarray:
    """Return selected atom indices and reject empty selections."""
    atom_indices = np.asarray(topology.select(selection), dtype=int)
    if atom_indices.size == 0:
        raise ValueError(f"No atoms selected for {label}: {selection}")
    return atom_indices
