"""External cluster profile loading and validation."""

from __future__ import annotations

import math
import os
import re
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

from .model import (
    SCHEMA_VERSION,
    fingerprint,
    positive_int as _positive_int,
    read_json,
    reject_unknown as _reject_unknown,
)

SCHEDULERS = {"local", "slurm", "age", "pjm"}
RESOURCE_CLASSES = {"light", "data", "quantum", "gpu", "md"}
_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
_PROFILE_KEYS = {
    "schema_version",
    "name",
    "scheduler",
    "work_root",
    "runner_argv",
    "environment",
    "policy",
    "resources",
}
_RESOURCE_KEYS = {
    "name",
    "classes",
    "capabilities",
    "limits",
    "scheduler_options",
    "incomplete",
}


def _finite_number(value: Any, label: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0):
        qualifier = "positive finite" if positive else "finite"
        raise ValueError(f"{label} must be a {qualifier} number")
    return result


def _safe_string(value: Any, label: str, *, single_line: bool = False) -> str:
    if not isinstance(value, str) or not value or "\x00" in value:
        raise ValueError(f"{label} must be a non-empty string")
    if single_line and any(character in value for character in "\r\n"):
        raise ValueError(f"{label} must be a single-line string")
    return value


def resolve_profile_path(value: str | None) -> Path | None:
    selected = value or os.environ.get("MDTBX_CLUSTER_PROFILE")
    if not selected:
        return None
    return Path(selected).expanduser().resolve()


def local_profile() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "name": "local",
        "scheduler": "local",
        "work_root": ".mdtbx/runs",
        "runner_argv": [sys.executable, "-m", "mdtbx"],
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
    result = deepcopy(data)
    _reject_unknown(result, _PROFILE_KEYS, "cluster profile")
    if result.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported cluster profile schema")
    _safe_string(result.get("name"), "Cluster profile name", single_line=True)
    if result.get("scheduler") not in SCHEDULERS:
        raise ValueError("scheduler must be local, slurm, age, or pjm")

    work_root = result.setdefault("work_root", ".mdtbx/runs")
    _safe_string(work_root, "work_root")
    runner = result.setdefault("runner_argv", [sys.executable, "-m", "mdtbx"])
    if not isinstance(runner, list) or not runner:
        raise ValueError("runner_argv must be a non-empty string list")
    for index, item in enumerate(runner):
        _safe_string(item, f"runner_argv[{index}]", single_line=True)

    resources = result.get("resources")
    if not isinstance(resources, list) or not resources:
        raise ValueError("Cluster profile requires at least one resource")
    names: set[str] = set()
    for resource in resources:
        if not isinstance(resource, dict):
            raise ValueError("Each resource must be an object")
        _reject_unknown(resource, _RESOURCE_KEYS, "resource")
        name = _safe_string(resource.get("name"), "Resource name", single_line=True)
        if not _NAME_PATTERN.fullmatch(name):
            raise ValueError(f"Invalid resource name: {name}")
        if name in names:
            raise ValueError(f"Duplicate resource: {name}")
        names.add(name)
        if resource.get("incomplete"):
            raise ValueError(f"Cluster resource is incomplete: {name}")

        classes = resource.get("classes")
        if classes is not None:
            if (
                not isinstance(classes, list)
                or not classes
                or not all(item in RESOURCE_CLASSES for item in classes)
            ):
                raise ValueError(f"{name}.classes contains an invalid resource class")

        capabilities = resource.get("capabilities")
        if not isinstance(capabilities, dict):
            raise ValueError(f"{name}.capabilities must be an object")
        _reject_unknown(
            capabilities,
            {
                "cpus_per_node",
                "gpus_per_node",
                "memory_mb_per_node",
                "gpu_type",
            },
            f"{name}.capabilities",
        )
        _positive_int(capabilities.get("cpus_per_node"), f"{name}.cpus_per_node")
        _positive_int(
            capabilities.get("gpus_per_node", 0),
            f"{name}.gpus_per_node",
            zero=True,
        )
        _positive_int(
            capabilities.get("memory_mb_per_node"), f"{name}.memory_mb_per_node"
        )
        gpu_type = capabilities.get("gpu_type")
        if gpu_type is not None:
            _safe_string(gpu_type, f"{name}.gpu_type", single_line=True)

        limits = resource.setdefault("limits", {})
        if not isinstance(limits, dict):
            raise ValueError(f"{name}.limits must be an object")
        _reject_unknown(limits, {"max_nodes", "max_walltime_seconds"}, f"{name}.limits")
        for key in ("max_nodes", "max_walltime_seconds"):
            if key in limits:
                _positive_int(limits[key], f"{name}.{key}")

        options = resource.setdefault("scheduler_options", {})
        if not isinstance(options, dict):
            raise ValueError(f"{name}.scheduler_options must be an object")
        extra = options.get("extra_directives", [])
        if not isinstance(extra, list):
            raise ValueError("extra_directives must be a list of strings")
        for index, item in enumerate(extra):
            _safe_string(item, f"extra_directives[{index}]", single_line=True)
        dependency_args = options.get("dependency_args")
        if dependency_args is not None:
            if not isinstance(dependency_args, list):
                raise ValueError("dependency_args must be a list of strings")
            for index, item in enumerate(dependency_args):
                _safe_string(item, f"dependency_args[{index}]", single_line=True)

    environment = result.setdefault("environment", {})
    if not isinstance(environment, dict):
        raise ValueError("environment must be an object")
    _reject_unknown(environment, {"prologue", "env"}, "environment")
    prologue = environment.setdefault("prologue", [])
    variables = environment.setdefault("env", {})
    if not isinstance(prologue, list):
        raise ValueError("environment.prologue must be a list of strings")
    for index, item in enumerate(prologue):
        _safe_string(item, f"environment.prologue[{index}]")
    if not isinstance(variables, dict):
        raise ValueError("environment.env must be a string mapping")
    for key, value in variables.items():
        if not isinstance(key, str) or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
            raise ValueError(f"Invalid environment variable name: {key}")
        if not isinstance(value, str) or "\x00" in value:
            raise ValueError("environment.env must be a string mapping")

    policy = result.setdefault("policy", {})
    if not isinstance(policy, dict):
        raise ValueError("policy must be an object")
    allowed_policy = {
        "batch_classes",
        "memory_safety_factor",
        "walltime_safety_factor",
        "minimum_confidence",
        "pilot_walltime_seconds",
        "poll_min_seconds",
        "poll_max_seconds",
    }
    _reject_unknown(policy, allowed_policy, "policy")
    batch_classes = policy.setdefault("batch_classes", ["md", "gpu", "quantum"])
    if not isinstance(batch_classes, list) or not all(
        item in RESOURCE_CLASSES for item in batch_classes
    ):
        raise ValueError("policy.batch_classes contains an invalid resource class")
    for key, default in (
        ("memory_safety_factor", 1.25),
        ("walltime_safety_factor", 1.5),
    ):
        policy[key] = _finite_number(
            policy.get(key, default), f"policy.{key}", positive=True
        )
    confidence = _finite_number(
        policy.get("minimum_confidence", 0.7), "policy.minimum_confidence"
    )
    if not 0 <= confidence <= 1:
        raise ValueError("policy.minimum_confidence must be between zero and one")
    policy["minimum_confidence"] = confidence
    for key, default in (
        ("pilot_walltime_seconds", 600),
        ("poll_min_seconds", 60),
        ("poll_max_seconds", 300),
    ):
        policy[key] = _positive_int(policy.get(key, default), f"policy.{key}")
    if policy["poll_max_seconds"] < policy["poll_min_seconds"]:
        raise ValueError("policy.poll_max_seconds must not be below poll_min_seconds")
    return result


def load_profile(value: str | None) -> tuple[Path | None, dict[str, Any]]:
    path = resolve_profile_path(value)
    if path is None:
        return None, validate_profile(local_profile())
    if not path.is_file():
        raise FileNotFoundError(f"Cluster profile not found: {path}")
    return path, validate_profile(read_json(path))


def profile_fingerprint(profile: dict[str, Any]) -> str:
    return fingerprint(profile)


def default_profile_dir() -> Path:
    root = os.environ.get("XDG_CONFIG_HOME")
    base = Path(root).expanduser() if root else Path.home() / ".config"
    return base / "mdtbx" / "clusters"
