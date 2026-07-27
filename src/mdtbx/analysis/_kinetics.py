"""Validated NumPy trajectory I/O for kinetic analysis commands."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def _split_concatenated(
    values: np.ndarray, lengths: np.ndarray, label: str
) -> list[np.ndarray]:
    if lengths.ndim != 1 or lengths.size == 0:
        raise ValueError(f"{label} lengths must be a non-empty one-dimensional array")
    if not np.issubdtype(lengths.dtype, np.integer) or np.any(lengths <= 0):
        raise ValueError(f"{label} lengths must contain positive integers")
    if int(lengths.sum()) != len(values):
        raise ValueError(f"{label} lengths do not match the concatenated data")
    boundaries = np.cumsum(lengths)[:-1]
    return list(np.split(values, boundaries))


def _feature_array(value: np.ndarray, label: str) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim != 2 or array.shape[0] == 0 or array.shape[1] == 0:
        raise ValueError(f"{label} must be a non-empty two-dimensional array")
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"{label} must contain numeric features")
    result = np.asarray(array, dtype=np.float64)
    if not np.isfinite(result).all():
        raise ValueError(f"{label} contains NaN or infinite values")
    return result


def load_feature_trajectories(
    paths: Iterable[str],
) -> tuple[list[np.ndarray], list[str]]:
    trajectories: list[np.ndarray] = []
    sources: list[str] = []
    feature_count: int | None = None
    for value in paths:
        path = Path(value).expanduser()
        sources.append(str(path.resolve()))
        if path.suffix == ".npy":
            loaded = [_feature_array(np.load(path, allow_pickle=False), str(path))]
        elif path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as archive:
                if "scores" not in archive or "lengths" not in archive:
                    raise ValueError(f"{path} is not a supported tICA feature archive")
                scores = _feature_array(archive["scores"], f"{path}:scores")
                loaded = _split_concatenated(scores, archive["lengths"], str(path))
        else:
            raise ValueError(f"Unsupported feature input suffix: {path.suffix}")
        for index, trajectory in enumerate(loaded):
            if feature_count is None:
                feature_count = trajectory.shape[1]
            elif trajectory.shape[1] != feature_count:
                raise ValueError(
                    f"Feature dimension mismatch in {path} trajectory {index}"
                )
            trajectories.append(trajectory)
    if not trajectories:
        raise ValueError("At least one feature trajectory is required")
    return trajectories, sources


def _discrete_array(value: np.ndarray, label: str) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{label} must be a non-empty one-dimensional array")
    if not np.issubdtype(array.dtype, np.integer):
        raise ValueError(f"{label} must contain integer state labels")
    result = np.asarray(array, dtype=np.int64)
    if np.any(result < 0):
        raise ValueError(f"{label} contains negative state labels")
    return result


def load_discrete_trajectories(
    paths: Iterable[str],
) -> tuple[list[np.ndarray], list[str]]:
    trajectories: list[np.ndarray] = []
    sources: list[str] = []
    for value in paths:
        path = Path(value).expanduser()
        sources.append(str(path.resolve()))
        if path.suffix == ".npy":
            loaded = [_discrete_array(np.load(path, allow_pickle=False), str(path))]
        elif path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as archive:
                if "dtrajs" not in archive or "lengths" not in archive:
                    raise ValueError(f"{path} is not a supported cluster archive")
                dtrajs = _discrete_array(archive["dtrajs"], f"{path}:dtrajs")
                loaded = _split_concatenated(dtrajs, archive["lengths"], str(path))
        else:
            raise ValueError(f"Unsupported discrete input suffix: {path.suffix}")
        trajectories.extend(loaded)
    if not trajectories:
        raise ValueError("At least one discrete trajectory is required")
    return trajectories, sources


def concatenate(trajectories: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    lengths = np.asarray(
        [len(trajectory) for trajectory in trajectories], dtype=np.int64
    )
    return np.concatenate(trajectories, axis=0), lengths
