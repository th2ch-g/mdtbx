"""External cluster profile loading and validation."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from .model import fingerprint, read_json

SCHEDULERS = {"local", "slurm", "age", "pjm"}


def resolve_profile_path(value: str | None) -> Path | None:
    selected = value or os.environ.get("MDTBX_CLUSTER_PROFILE")
    if not selected:
        return None
    return Path(selected).expanduser().resolve()


def local_profile() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "name": "local",
        "scheduler": "local",
        "work_root": ".mdtbx/runs",
        "environment": {"prologue": [], "env": {}},
        "policy": {
            "batch_classes": ["data", "quantum", "gpu", "md"],
            "memory_safety_factor": 1.25,
            "walltime_safety_factor": 1.5,
            "minimum_confidence": 0.7,
            "pilot_walltime_seconds": 600,
            "poll_min_seconds": 60,
            "poll_max_seconds": 300,
        },
        "resources": [
            {
                "name": "local",
                "capabilities": {
                    "cpus_per_node": 1,
                    "gpus_per_node": 0,
                    "memory_mb_per_node": 4096,
                },
                "limits": {"max_nodes": 1, "max_walltime_seconds": 86400},
                "scheduler_options": {},
            }
        ],
    }


def validate_profile(data: Any) -> dict[str, Any]:
    if not isinstance(data, dict):
        raise ValueError("Cluster profile must be a JSON object")
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported cluster profile schema")
    if not isinstance(data.get("name"), str) or not data["name"]:
        raise ValueError("Cluster profile requires a name")
    if data.get("scheduler") not in SCHEDULERS:
        raise ValueError("scheduler must be local, slurm, age, or pjm")
    resources = data.get("resources")
    if not isinstance(resources, list) or not resources:
        raise ValueError("Cluster profile requires at least one resource")
    names = set()
    for resource in resources:
        if not isinstance(resource, dict) or not isinstance(resource.get("name"), str):
            raise ValueError("Each resource requires a name")
        if resource["name"] in names:
            raise ValueError(f"Duplicate resource: {resource['name']}")
        names.add(resource["name"])
        if resource.get("incomplete"):
            raise ValueError(f"Cluster resource is incomplete: {resource['name']}")
        capabilities = resource.get("capabilities", {})
        for key in ("cpus_per_node", "gpus_per_node", "memory_mb_per_node"):
            value = capabilities.get(key, 0)
            if not isinstance(value, int) or value < 0:
                raise ValueError(f"{resource['name']}.{key} must be non-negative")
        for key in ("cpus_per_node", "memory_mb_per_node"):
            if capabilities.get(key, 0) == 0:
                raise ValueError(f"{resource['name']}.{key} must be positive")
        extra = resource.get("scheduler_options", {}).get("extra_directives", [])
        if not isinstance(extra, list) or not all(
            isinstance(item, str) for item in extra
        ):
            raise ValueError("extra_directives must be a list of strings")
    environment = data.setdefault("environment", {})
    prologue = environment.setdefault("prologue", [])
    variables = environment.setdefault("env", {})
    if not isinstance(prologue, list) or not all(
        isinstance(item, str) for item in prologue
    ):
        raise ValueError("environment.prologue must be a list of strings")
    if not isinstance(variables, dict) or not all(
        isinstance(key, str) and isinstance(value, str)
        for key, value in variables.items()
    ):
        raise ValueError("environment.env must be a string mapping")
    data.setdefault("work_root", ".mdtbx/runs")
    data.setdefault("policy", {})
    return data


def load_profile(value: str | None) -> tuple[Path | None, dict[str, Any]]:
    path = resolve_profile_path(value)
    if path is None:
        return None, local_profile()
    if not path.is_file():
        raise FileNotFoundError(f"Cluster profile not found: {path}")
    return path, validate_profile(read_json(path))


def profile_fingerprint(profile: dict[str, Any]) -> str:
    return fingerprint(profile)


def default_profile_dir() -> Path:
    root = os.environ.get("XDG_CONFIG_HOME")
    base = Path(root).expanduser() if root else Path.home() / ".config"
    return base / "mdtbx" / "clusters"
