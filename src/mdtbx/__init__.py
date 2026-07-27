"""Molecular-dynamics toolbox."""

import os
from pathlib import Path


def prepare_runtime() -> None:
    """Expose project-local executables when running from a source checkout."""
    project_root = Path(__file__).resolve().parents[2]
    pixi_bin = project_root / ".pixi/envs/default/bin"
    current_path = os.environ.get("PATH", "")
    path_entries = current_path.split(os.pathsep) if current_path else []
    if pixi_bin.is_dir() and str(pixi_bin) not in path_entries:
        os.environ["PATH"] = os.pathsep.join([str(pixi_bin), *path_entries])
