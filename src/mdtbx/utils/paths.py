"""Shared filesystem path helpers."""

from pathlib import Path


def ensure_output_parent(path: str | Path) -> Path:
    """Expand a user path and create its parent directory."""
    output_path = Path(path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path
