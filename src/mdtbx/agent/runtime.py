"""Planning, resource selection, execution, and run-state collection."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shlex
import stat
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any

from .json_schemas import schemas as protocol_schemas
from .model import (
    SCHEMA_VERSION,
    finalize_plan,
    fingerprint,
    json_value,
    positive_int as _positive_int,
    read_json,
    reject_unknown as _reject_unknown,
    utc_now,
    verify_plan,
    write_json,
    write_text as _write_text,
)
from .profile import load_profile, profile_fingerprint
from .registry import (
    AGENT_COMMANDS,
    all_descriptors,
    arguments_to_argv,
    artifact_paths,
    command_parser,
    input_paths,
    normalized_arguments,
)
from .schedulers import detect_scheduler, scheduler

DEFAULT_MEMORY_MB = {
    "light": 1024,
    "data": 4096,
    "quantum": 8192,
    "gpu": 16384,
    "md": 8192,
}
DEFAULT_WALLTIME_SECONDS = {
    "light": 600,
    "data": 3600,
    "quantum": 14400,
    "gpu": 14400,
    "md": 21600,
}
_STEP_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
_REQUEST_KEYS = {
    "schema_version",
    "name",
    "cwd",
    "cluster_profile",
    "execution",
    "steps",
}
_STEP_KEYS = {
    "id",
    "command",
    "arguments",
    "depends_on",
    "resources",
    "evidence",
    "confidence",
}
_RESOURCE_KEYS = {
    "resource",
    "nodes",
    "cpus_per_node",
    "tasks_per_node",
    "gpus_per_node",
    "memory_mb",
    "walltime_seconds",
}


def _root_parser(commands: set[str] | None = None):
    from ..cli import create_parser

    return create_parser(commands)


def schema(command: str | None = None) -> dict[str, Any]:
    commands = all_descriptors(_root_parser({command} if command is not None else None))
    if command is not None:
        try:
            commands = {command: commands[command]}
        except KeyError as error:
            raise ValueError(f"Unknown mdtbx command: {command}") from error
    return {
        "schema_version": SCHEMA_VERSION,
        "schemas": protocol_schemas(),
        "commands": commands,
    }


def _validate_resources(identifier: str, value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError(f"{identifier}.resources must be an object")
    _reject_unknown(value, _RESOURCE_KEYS, f"{identifier}.resources")
    result = dict(value)
    if "resource" in result:
        resource = result["resource"]
        if not isinstance(resource, str) or not resource or "\x00" in resource:
            raise ValueError(f"{identifier}.resources.resource must be a string")
    for key in (
        "nodes",
        "cpus_per_node",
        "tasks_per_node",
        "memory_mb",
        "walltime_seconds",
    ):
        if key in result:
            _positive_int(result[key], f"{identifier}.resources.{key}")
    if "gpus_per_node" in result:
        _positive_int(
            result["gpus_per_node"],
            f"{identifier}.resources.gpus_per_node",
            zero=True,
        )
    if (
        "tasks_per_node" in result
        and "cpus_per_node" in result
        and result["tasks_per_node"] > result["cpus_per_node"]
    ):
        raise ValueError(
            f"{identifier}.resources.tasks_per_node must not exceed cpus_per_node"
        )
    return result


def _validated_steps(request: dict[str, Any]) -> list[dict[str, Any]]:
    raw_steps = request.get("steps")
    if not isinstance(raw_steps, list) or not raw_steps:
        raise ValueError("Request requires a non-empty steps list")
    command_names = {
        raw.get("command")
        for raw in raw_steps
        if isinstance(raw, dict) and isinstance(raw.get("command"), str)
    }
    root = _root_parser(command_names)
    descriptors = all_descriptors(root)
    steps = []
    identifiers: set[str] = set()
    for raw in raw_steps:
        if not isinstance(raw, dict):
            raise ValueError("Each request step must be an object")
        _reject_unknown(raw, _STEP_KEYS, "step")
        identifier = raw.get("id")
        command = raw.get("command")
        arguments = raw.get("arguments", {})
        if not isinstance(identifier, str) or not _STEP_ID.fullmatch(identifier):
            raise ValueError(f"Invalid step id: {identifier!r}")
        if identifier in {".", ".."}:
            raise ValueError(f"Invalid step id: {identifier!r}")
        if identifier in identifiers:
            raise ValueError(f"Duplicate step id: {identifier}")
        if not isinstance(command, str) or command not in descriptors:
            raise ValueError(f"Unknown mdtbx command: {command}")
        if command in AGENT_COMMANDS:
            raise ValueError("Agent commands cannot be nested in an agent plan")
        if not isinstance(arguments, dict):
            raise ValueError(f"{identifier}.arguments must be an object")
        resource_evidence = raw.get("evidence", [])
        if not isinstance(resource_evidence, list) or not all(
            isinstance(item, dict) for item in resource_evidence
        ):
            raise ValueError(f"{identifier}.evidence must be an object list")
        requested_confidence = raw.get("confidence")
        if requested_confidence is not None:
            if (
                isinstance(requested_confidence, bool)
                or not isinstance(requested_confidence, (int, float))
                or not math.isfinite(float(requested_confidence))
                or not 0 <= requested_confidence <= 1
            ):
                raise ValueError(
                    f"{identifier}.confidence must be between zero and one"
                )
        parser = command_parser(root, command)
        argv = arguments_to_argv(parser, arguments)
        normalized = normalized_arguments(parser, argv)
        depends_on = raw.get("depends_on", [])
        if not isinstance(depends_on, list) or not all(
            isinstance(item, str) and _STEP_ID.fullmatch(item) for item in depends_on
        ):
            raise ValueError(f"{identifier}.depends_on must be a step-id list")
        if len(depends_on) != len(set(depends_on)):
            raise ValueError(f"{identifier}.depends_on contains duplicates")
        steps.append(
            {
                "id": identifier,
                "command": command,
                "arguments": normalized,
                "argv": argv,
                "depends_on": depends_on,
                "descriptor": descriptors[command],
                "requested_resources": _validate_resources(
                    identifier, raw.get("resources")
                ),
                "resource_evidence": resource_evidence,
                "requested_confidence": requested_confidence,
            }
        )
        identifiers.add(identifier)
    for step in steps:
        missing = sorted(set(step["depends_on"]) - identifiers)
        if missing:
            raise ValueError(
                f"{step['id']} depends on unknown steps: {', '.join(missing)}"
            )
    _topological(steps)
    return steps


def _topological(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    remaining = {step["id"]: step for step in steps}
    completed: set[str] = set()
    ordered = []
    while remaining:
        ready = sorted(
            (
                step
                for step in remaining.values()
                if set(step["depends_on"]) <= completed
            ),
            key=lambda item: item["id"],
        )
        if not ready:
            raise ValueError("Workflow dependencies contain a cycle")
        for step in ready:
            ordered.append(step)
            completed.add(step["id"])
            del remaining[step["id"]]
    return ordered


def _absolute_path(cwd: Path, value: str) -> Path:
    expanded = Path(value).expanduser()
    if not expanded.is_absolute():
        expanded = cwd / expanded
    return Path(os.path.abspath(expanded))


def _snapshot(path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {"path": str(path), "exists": os.path.lexists(path)}
    if not result["exists"]:
        return result
    details = path.lstat()
    if stat.S_ISLNK(details.st_mode):
        kind = "symlink"
        result["symlink_target"] = os.readlink(path)
    elif stat.S_ISREG(details.st_mode):
        kind = "file"
    elif stat.S_ISDIR(details.st_mode):
        kind = "directory"
    else:
        kind = "other"
    result.update(
        {
            "kind": kind,
            "size": details.st_size,
            "mtime_ns": details.st_mtime_ns,
        }
    )
    return result


def _snapshot_identity(snapshot: dict[str, Any]) -> dict[str, Any]:
    return {
        key: snapshot.get(key)
        for key in ("path", "exists", "kind", "size", "mtime_ns", "symlink_target")
        if key in snapshot
    }


def _verify_snapshots(plan: dict[str, Any]) -> None:
    for step in plan["steps"]:
        for group in ("artifact_snapshots", "destructive_snapshots"):
            for approved in step.get(group, []):
                current = _snapshot(Path(approved["path"]))
                if _snapshot_identity(current) != _snapshot_identity(approved):
                    raise ValueError(
                        "Filesystem changed after plan approval target was captured: "
                        f"{approved['path']}; create and approve a new plan"
                    )


def _input_bytes(step: dict[str, Any], cwd: Path) -> int:
    total = 0
    for name in step["descriptor"]["inputs"]:
        value = step["arguments"].get(name)
        values = value if isinstance(value, list) else [value]
        for item in values:
            if not isinstance(item, str):
                continue
            path = _absolute_path(cwd, item)
            if path.is_file():
                total += path.stat().st_size
    return total


def _policy(profile: dict[str, Any], key: str, default: Any) -> Any:
    return profile.get("policy", {}).get(key, default)


def _resource_candidates(
    profile: dict[str, Any], resource_class: str
) -> list[dict[str, Any]]:
    # Resources without a classes list are wildcards that match every class.
    return [
        resource
        for resource in profile["resources"]
        if not resource.get("classes") or resource_class in resource["classes"]
    ]


def _select_resources(
    step: dict[str, Any], profile: dict[str, Any], cwd: Path
) -> tuple[dict[str, Any], dict[str, Any], float, list[dict[str, Any]]]:
    requested = step.get("requested_resources")
    resource_class = step["descriptor"]["resource_class"]
    input_bytes = _input_bytes(step, cwd)
    evidence = [
        {"source": "command_registry", "resource_class": resource_class},
        {"source": "input_files", "total_bytes": input_bytes},
    ]
    evidence.extend(json_value(step["resource_evidence"]))

    explicit = requested is not None
    memory_estimate = DEFAULT_MEMORY_MB[resource_class]
    if resource_class == "data":
        memory_estimate = max(memory_estimate, int(input_bytes * 3 / 1_000_000) + 1024)
    memory_factor = float(_policy(profile, "memory_safety_factor", 1.25))
    walltime_factor = float(_policy(profile, "walltime_safety_factor", 1.5))
    memory_mb = max(1, int(memory_estimate * memory_factor))
    walltime = max(1, int(DEFAULT_WALLTIME_SECONDS[resource_class] * walltime_factor))
    nodes = 1
    minimum_gpus = 1 if resource_class == "gpu" else 0
    minimum_cpus = 1
    tasks_per_node = 1
    resource_name = None
    if requested is not None:
        nodes = requested.get("nodes", nodes)
        memory_mb = requested.get("memory_mb", memory_mb)
        walltime = requested.get("walltime_seconds", walltime)
        minimum_gpus = requested.get("gpus_per_node", minimum_gpus)
        minimum_cpus = requested.get("cpus_per_node", minimum_cpus)
        tasks_per_node = requested.get("tasks_per_node", tasks_per_node)
        minimum_cpus = max(minimum_cpus, tasks_per_node)
        resource_name = requested.get("resource")
        evidence.append({"source": "request", "resources": json_value(requested)})

    candidates = _resource_candidates(profile, resource_class)
    suitable = []
    for candidate in candidates:
        if resource_name and candidate["name"] != resource_name:
            continue
        capabilities = candidate["capabilities"]
        if capabilities["memory_mb_per_node"] < memory_mb:
            continue
        if capabilities.get("gpus_per_node", 0) < minimum_gpus:
            continue
        if capabilities["cpus_per_node"] < minimum_cpus:
            continue
        limits = candidate.get("limits", {})
        if nodes > limits.get("max_nodes", nodes):
            continue
        if walltime > limits.get("max_walltime_seconds", walltime):
            continue
        score = (
            capabilities.get("gpus_per_node", 0) * 10**12
            + capabilities["memory_mb_per_node"] * 10**4
            + capabilities["cpus_per_node"]
        )
        suitable.append((score, candidate))
    if not suitable:
        raise ValueError(f"No cluster resource can satisfy step {step['id']}")
    resource = min(suitable, key=lambda item: item[0])[1]
    if minimum_cpus % tasks_per_node:
        raise ValueError(
            f"Step {step['id']} cpus_per_node must be divisible by tasks_per_node"
        )
    allocation = {
        "resource": resource["name"],
        "nodes": nodes,
        "cpus_per_node": minimum_cpus,
        "tasks_per_node": tasks_per_node,
        "gpus_per_node": minimum_gpus,
        "memory_mb": memory_mb,
        "walltime_seconds": walltime,
    }
    requested_confidence = step["requested_confidence"]
    trusted_sources = {"approved_pilot", "prior_run", "benchmark"}
    has_measured_evidence = any(
        item.get("source") in trusted_sources for item in step["resource_evidence"]
    )
    if requested_confidence is not None:
        confidence = float(requested_confidence)
    elif has_measured_evidence:
        confidence = 0.85
    else:
        confidence = (
            0.65 if explicit else (0.6 if profile["scheduler"] != "local" else 0.5)
        )
    evidence.append(
        {
            "source": "cluster_profile",
            "selected": resource["name"],
            "capabilities": resource["capabilities"],
            "safety_factors": {"memory": memory_factor, "walltime": walltime_factor},
        }
    )
    return allocation, resource, confidence, evidence


def _execution_mode(
    request: dict[str, Any], steps: list[dict[str, Any]], profile: dict[str, Any]
) -> str:
    requested = request.get("execution", "auto")
    if requested not in {"auto", "local", "batch"}:
        raise ValueError("execution must be auto, local, or batch")
    batch_classes = set(_policy(profile, "batch_classes", ["md", "gpu", "quantum"]))
    requires_batch = any(
        step["descriptor"]["resource_class"] in batch_classes for step in steps
    )
    if requested == "local" and requires_batch:
        classes = sorted(
            {
                step["descriptor"]["resource_class"]
                for step in steps
                if step["descriptor"]["resource_class"] in batch_classes
            }
        )
        raise ValueError(
            "Local execution is prohibited for configured batch classes: "
            + ", ".join(classes)
        )
    selected = (
        "batch"
        if requested == "batch" or (requested == "auto" and requires_batch)
        else "local"
    )
    if selected == "batch" and profile["scheduler"] == "local":
        raise ValueError("Batch execution requires a non-local cluster profile")
    return selected


def build_plan(request: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(request, dict):
        raise ValueError("Request must be a JSON object")
    _reject_unknown(request, _REQUEST_KEYS, "request")
    if request.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported request schema")
    cwd_value = request.get("cwd", ".")
    if not isinstance(cwd_value, str) or not cwd_value or "\x00" in cwd_value:
        raise ValueError("cwd must be a non-empty path string")
    cwd = Path(cwd_value).expanduser().resolve()
    if not cwd.is_dir():
        raise FileNotFoundError(f"Working directory not found: {cwd}")
    name = request.get("name", "mdtbx-agent-plan")
    if not isinstance(name, str) or not name or "\x00" in name:
        raise ValueError("name must be a non-empty string")
    profile_value = request.get("cluster_profile")
    if profile_value is not None and not isinstance(profile_value, str):
        raise ValueError("cluster_profile must be a path string or null")
    profile_path, profile = load_profile(profile_value)
    steps = _validated_steps(request)
    execution = _execution_mode(request, steps, profile)

    planned_steps = []
    plan_kind = "production"
    minimum_confidence = float(_policy(profile, "minimum_confidence", 0.7))
    for step in steps:
        allocation, resource, confidence, evidence = _select_resources(
            step, profile, cwd
        )
        arguments = dict(step["arguments"])
        argv = list(step["argv"])
        if confidence < minimum_confidence and step["descriptor"]["pilot_capable"]:
            plan_kind = "pilot"
            requested_resources = step.get("requested_resources") or {}
            if "walltime_seconds" not in requested_resources:
                allocation["walltime_seconds"] = min(
                    allocation["walltime_seconds"],
                    int(_policy(profile, "pilot_walltime_seconds", 600)),
                )
            if (
                step["command"] in {"run_fep", "run_abfe"}
                and arguments.get("nsteps") is None
            ):
                arguments["nsteps"] = 1000
                parser = command_parser(
                    _root_parser({step["command"]}), step["command"]
                )
                argv = arguments_to_argv(parser, arguments)
                arguments = normalized_arguments(parser, argv)
                evidence.append(
                    {"source": "pilot_policy", "override": {"nsteps": 1000}}
                )

        inputs = [
            str(_absolute_path(cwd, value))
            for value in input_paths(arguments, step["descriptor"])
        ]
        outputs = [
            str(_absolute_path(cwd, value))
            for value in artifact_paths(arguments, step["descriptor"])
        ]
        artifact_snapshots = [_snapshot(Path(value)) for value in outputs]
        existing_artifacts = [
            item["path"] for item in artifact_snapshots if item["exists"]
        ]
        destructive_targets = (
            inputs if step["descriptor"]["risk"] == "destructive" else []
        )
        destructive_snapshots = [
            _snapshot(Path(value)) for value in destructive_targets
        ]
        planned_steps.append(
            {
                "id": step["id"],
                "command": step["command"],
                "arguments": arguments,
                "argv": argv,
                "depends_on": step["depends_on"],
                "execution": execution,
                "resource_class": step["descriptor"]["resource_class"],
                "risk": step["descriptor"]["risk"],
                "resources": allocation,
                "resource_profile": resource,
                "confidence": confidence,
                "evidence": evidence,
                "inputs": inputs,
                "artifacts": outputs,
                "artifact_snapshots": artifact_snapshots,
                "existing_artifacts": existing_artifacts,
                "destructive_targets": destructive_targets,
                "destructive_snapshots": destructive_snapshots,
            }
        )

    result = {
        "schema_version": SCHEMA_VERSION,
        "created_at": utc_now(),
        "name": name,
        "plan_kind": plan_kind,
        "cwd": str(cwd),
        "cluster_profile": str(profile_path) if profile_path else None,
        "profile_fingerprint": profile_fingerprint(profile),
        "scheduler": profile["scheduler"],
        "execution": execution,
        "runner_argv": profile["runner_argv"],
        "work_root": profile["work_root"],
        "steps": planned_steps,
        "approval": {
            "required": True,
            "unsafe": any(step["risk"] == "unsafe" for step in planned_steps),
            "destructive": any(step["risk"] == "destructive" for step in planned_steps),
            "overwrite": any(step["existing_artifacts"] for step in planned_steps),
        },
    }
    return finalize_plan(result)


def direct_plan(command: str, namespace: Any) -> dict[str, Any]:
    values = {
        key: json_value(value)
        for key, value in vars(namespace).items()
        if key
        not in {
            "_command",
            "func",
            "_agent_json",
            "_agent_dry_run",
            "_agent_cluster_profile",
        }
    }
    request = {
        "schema_version": SCHEMA_VERSION,
        "name": f"dry-run-{command}",
        "cwd": str(Path.cwd()),
        "cluster_profile": getattr(namespace, "_agent_cluster_profile", None),
        "execution": "auto",
        "steps": [
            {
                "id": command,
                "command": command,
                "arguments": values,
                "depends_on": [],
            }
        ],
    }
    return build_plan(request)


def _run_root(plan: dict[str, Any]) -> Path:
    root = Path(plan["work_root"]).expanduser()
    if not root.is_absolute():
        root = Path(plan["cwd"]) / root
    timestamp = utc_now().replace(":", "").replace("-", "").replace(".", "")
    run_id = f"{timestamp}-{plan['plan_id'][:10]}-{uuid.uuid4().hex[:8]}"
    return root / run_id


def _append_event(run_dir: Path, event: dict[str, Any]) -> None:
    path = run_dir / "events.jsonl"
    payload = json.dumps(json_value(event), ensure_ascii=False) + "\n"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def _persist_state(run_dir: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = utc_now()
    write_json(run_dir / "state.json", state)


def _render_script(
    plan: dict[str, Any],
    step: dict[str, Any],
    run_dir: Path,
    profile: dict[str, Any],
) -> Path:
    adapter = scheduler(plan["scheduler"])
    step_dir = run_dir / "steps" / step["id"]
    step_dir.mkdir(parents=True, exist_ok=True)
    resource = step["resource_profile"]
    directives = adapter.directives(step, resource, profile)
    command = [*plan["runner_argv"], step["command"], *step["argv"], "--json"]
    environment = profile.get("environment", {})
    lines = ["#!/bin/bash", *directives, "set -euo pipefail"]
    for key, value in environment.get("env", {}).items():
        lines.append(f"export {key}={shlex.quote(value)}")
    lines.extend(environment.get("prologue", []))
    lines.append(f"cd {shlex.quote(plan['cwd'])}")
    lines.append(
        f"{shlex.join(command)} > {shlex.quote(str(step_dir / 'result.json'))} "
        f"2> {shlex.quote(str(step_dir / 'stderr.log'))}"
    )
    script = step_dir / "job.sh"
    _write_text(script, "\n".join(lines) + "\n")
    script.chmod(0o700)
    return script


def _validate_executable_plan(plan: dict[str, Any]) -> str:
    if not isinstance(plan, dict) or plan.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported plan schema")
    plan_id = verify_plan(plan)
    if plan.get("execution") not in {"local", "batch"}:
        raise ValueError("Plan execution must be local or batch")
    if not isinstance(plan.get("steps"), list) or not plan["steps"]:
        raise ValueError("Plan requires steps")
    identifiers = set()
    for step in plan["steps"]:
        identifier = step.get("id")
        if not isinstance(identifier, str) or not _STEP_ID.fullmatch(identifier):
            raise ValueError(f"Invalid planned step id: {identifier!r}")
        if identifier in identifiers:
            raise ValueError(f"Duplicate planned step id: {identifier}")
        identifiers.add(identifier)
    _topological(plan["steps"])
    return plan_id


def execute_plan(
    plan: dict[str, Any],
    *,
    approval: str,
    approve_unsafe: bool = False,
    approve_destructive: bool = False,
    approve_overwrite: bool = False,
) -> dict[str, Any]:
    plan_id = _validate_executable_plan(plan)
    if approval != plan_id:
        raise ValueError("--approve must exactly match plan_id")
    if plan["approval"]["unsafe"] and not approve_unsafe:
        raise ValueError("Unsafe plan requires --approve-unsafe")
    if plan["approval"]["destructive"] and not approve_destructive:
        raise ValueError("Destructive plan requires --approve-destructive")
    if plan["approval"].get("overwrite") and not approve_overwrite:
        raise ValueError("Existing artifacts require --approve-overwrite")
    profile_path, profile = load_profile(plan.get("cluster_profile"))
    if profile_fingerprint(profile) != plan["profile_fingerprint"]:
        raise ValueError("Cluster profile changed after plan creation")
    _verify_snapshots(plan)

    run_dir = _run_root(plan)
    run_dir.mkdir(parents=True, exist_ok=False)
    write_json(run_dir / "plan.json", plan)
    state = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_dir.name,
        "run_dir": str(run_dir),
        "plan_id": plan_id,
        "scheduler": plan["scheduler"] if plan["execution"] == "batch" else "local",
        "cluster_profile": str(profile_path) if profile_path else None,
        "status": "running",
        "created_at": utc_now(),
        "steps": {step["id"]: {"state": "pending"} for step in plan["steps"]},
    }
    _persist_state(run_dir, state)
    _append_event(run_dir, {"at": utc_now(), "event": "approved", "plan_id": plan_id})

    ordered = _topological(plan["steps"])
    if plan["execution"] == "local":
        for step in ordered:
            step_dir = run_dir / "steps" / step["id"]
            step_dir.mkdir(parents=True, exist_ok=True)
            command = [*plan["runner_argv"], step["command"], *step["argv"], "--json"]
            state["steps"][step["id"]] = {
                "state": "running",
                "started_at": utc_now(),
                "argv": command,
            }
            _persist_state(run_dir, state)
            _append_event(
                run_dir,
                {"at": utc_now(), "event": "step_started", "step": step["id"]},
            )
            try:
                result = subprocess.run(
                    command,
                    cwd=plan["cwd"],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                _write_text(step_dir / "result.json", result.stdout)
                _write_text(step_dir / "stderr.log", result.stderr)
                status = "succeeded" if result.returncode == 0 else "failed"
                state["steps"][step["id"]].update(
                    {
                        "state": status,
                        "returncode": result.returncode,
                        "finished_at": utc_now(),
                    }
                )
                _persist_state(run_dir, state)
                _append_event(
                    run_dir,
                    {
                        "at": utc_now(),
                        "event": "step_finished",
                        "step": step["id"],
                        "state": status,
                    },
                )
                if result.returncode:
                    state["status"] = "failed"
                    _persist_state(run_dir, state)
                    break
            except BaseException as error:
                state["steps"][step["id"]].update(
                    {
                        "state": "failed",
                        "finished_at": utc_now(),
                        "error": {"type": type(error).__name__, "message": str(error)},
                    }
                )
                state["status"] = "failed"
                _persist_state(run_dir, state)
                _append_event(
                    run_dir,
                    {"at": utc_now(), "event": "step_failed", "step": step["id"]},
                )
                raise
        else:
            state["status"] = "succeeded"
            _persist_state(run_dir, state)
    else:
        adapter = scheduler(plan["scheduler"])
        submitted: dict[str, str] = {}
        for step in ordered:
            try:
                dependencies = [submitted[item] for item in step["depends_on"]]
                script = _render_script(plan, step, run_dir, profile)
                state["steps"][step["id"]] = {
                    "state": "submitting",
                    "script": str(script),
                }
                _persist_state(run_dir, state)
                job_id = adapter.submit(
                    script,
                    dependencies=dependencies,
                    resource=step["resource_profile"],
                )
                submitted[step["id"]] = job_id
                state["steps"][step["id"]] = {
                    "state": "queued",
                    "job_id": job_id,
                    "script": str(script),
                    "submitted_at": utc_now(),
                }
                _persist_state(run_dir, state)
                _append_event(
                    run_dir,
                    {
                        "at": utc_now(),
                        "event": "job_submitted",
                        "step": step["id"],
                        "job_id": job_id,
                    },
                )
            except BaseException as error:
                state["steps"][step["id"]] = {
                    **state["steps"].get(step["id"], {}),
                    "state": "failed",
                    "failed_at": utc_now(),
                    "error": {"type": type(error).__name__, "message": str(error)},
                }
                state["status"] = "submission_failed"
                _persist_state(run_dir, state)
                _append_event(
                    run_dir,
                    {
                        "at": utc_now(),
                        "event": "submission_failed",
                        "step": step["id"],
                    },
                )
                raise
        state["status"] = "submitted"
        _persist_state(run_dir, state)
    return state


def load_run(value: str) -> tuple[Path, dict[str, Any]]:
    path = Path(value).expanduser()
    if path.is_dir():
        run_dir = path.resolve()
    else:
        roots = (Path.cwd(), *Path.cwd().parents)
        matches = [
            item
            for root in roots
            for item in (root / ".mdtbx" / "runs").glob("*/state.json")
        ]
        selected = [
            item.parent
            for item in matches
            if item.parent.name == value or read_json(item).get("run_id") == value
        ]
        if len(selected) != 1:
            raise FileNotFoundError(f"Run not found: {value}")
        run_dir = selected[0].resolve()
    state_path = run_dir / "state.json"
    if not state_path.is_file():
        raise FileNotFoundError(f"Run state not found: {state_path}")
    state = read_json(state_path)
    if state.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported run state schema")
    return run_dir, state


def _overall_status(states: set[str]) -> str:
    if states and states <= {"succeeded"}:
        return "succeeded"
    if states & {"failed", "cancelled"}:
        return "failed"
    if "running" in states:
        return "running"
    if states & {"queued", "submitting", "pending"}:
        return "queued"
    return "unknown"


def run_status(value: str) -> dict[str, Any]:
    run_dir, state = load_run(value)
    if state["scheduler"] == "local":
        return state
    adapter = scheduler(state["scheduler"])
    statuses = {}
    for step, item in state["steps"].items():
        if item.get("job_id"):
            statuses[step] = adapter.status(item["job_id"])
        else:
            statuses[step] = {
                "schema_version": SCHEMA_VERSION,
                "job_id": None,
                "state": item.get("state", "unknown"),
                "raw_state": "",
                "reason": "not submitted",
                "elapsed": "",
                "checked_at": utc_now(),
            }
    overall = _overall_status({item["state"] for item in statuses.values()})
    return {
        **state,
        "run_dir": str(run_dir),
        "status": overall,
        "checked_at": utc_now(),
        "scheduler_status": statuses,
    }


def cancel_run(value: str, approval: str) -> dict[str, Any]:
    """Cancel queued or running jobs from one immutable run."""
    run_dir, state = load_run(value)
    if approval != state["plan_id"]:
        raise ValueError("Approval must exactly match the run plan_id")
    if state["scheduler"] == "local":
        raise ValueError("Local runs cannot be cancelled through a batch scheduler")

    adapter = scheduler(state["scheduler"])
    cancelled = {}
    for step, item in state["steps"].items():
        job_id = item.get("job_id")
        if not job_id:
            continue
        current = adapter.status(job_id)
        item["state"] = current["state"]
        if current["state"] not in {"queued", "running"}:
            continue
        adapter.cancel(job_id)
        item["state"] = "cancelled"
        cancelled[step] = str(job_id)

    if cancelled:
        state["status"] = "cancelled"
        state["updated_at"] = utc_now()
        _persist_state(run_dir, state)
        _append_event(
            run_dir,
            {"at": utc_now(), "event": "cancelled", "steps": cancelled},
        )
    return {
        **state,
        "run_dir": str(run_dir),
        "cancelled": cancelled,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_result(path: str) -> dict[str, Any]:
    snapshot = _snapshot(Path(path))
    if snapshot.get("kind") == "file":
        snapshot["sha256"] = _sha256(Path(path))
    return snapshot


def collect_run(value: str) -> dict[str, Any]:
    run_dir, state = load_run(value)
    status = run_status(str(run_dir))
    plan = read_json(run_dir / "plan.json")
    plan_steps = {step["id"]: step for step in plan["steps"]}
    results = {}
    for step in state["steps"]:
        result_path = run_dir / "steps" / step / "result.json"
        entry: dict[str, Any] = {
            "result_path": str(result_path),
            "valid": False,
            "result": None,
            "artifacts": [
                _artifact_result(path) for path in plan_steps[step].get("artifacts", [])
            ],
        }
        if result_path.is_file() and result_path.read_text().strip():
            try:
                payload = read_json(result_path)
                entry["result"] = payload
                entry["valid"] = (
                    isinstance(payload, dict)
                    and payload.get("schema_version") == 1
                    and isinstance(payload.get("ok"), bool)
                    and isinstance(payload.get("exit_code"), int)
                )
            except json.JSONDecodeError as error:
                entry["error"] = {"type": type(error).__name__, "message": str(error)}
        results[step] = entry
    result = {
        "schema_version": SCHEMA_VERSION,
        "run_id": state["run_id"],
        "plan_id": state["plan_id"],
        "status": status["status"],
        "collected_at": utc_now(),
        "scheduler_status": status.get("scheduler_status", {}),
        "steps": results,
    }
    write_json(run_dir / "result.json", result)
    for step, scheduler_state in status.get("scheduler_status", {}).items():
        if step in state["steps"]:
            state["steps"][step]["state"] = scheduler_state["state"]
    state["status"] = status["status"]
    state["collected_at"] = result["collected_at"]
    _persist_state(run_dir, state)
    _append_event(
        run_dir, {"at": utc_now(), "event": "collected", "status": result["status"]}
    )
    return result


def probe_cluster(
    scheduler_name: str | None, profile_value: str | None
) -> dict[str, Any]:
    profile_path = None
    profile = None
    if profile_value:
        profile_path, profile = load_profile(profile_value)
        scheduler_name = scheduler_name or profile["scheduler"]
    detected, matches = detect_scheduler()
    selected = scheduler_name or detected
    if selected is None:
        if matches:
            raise ValueError(
                "Multiple schedulers detected; select one explicitly: "
                + ", ".join(matches)
            )
        raise ValueError("No supported scheduler detected")
    if selected == "local":
        raise ValueError("agent_probe requires slurm, age, or pjm")
    result = scheduler(selected).probe(profile)
    result["detected_schedulers"] = matches
    result["cluster_profile"] = str(profile_path) if profile_path else None
    result["draft_profile"] = draft_profile(result)
    result["draft_id"] = fingerprint(result["draft_profile"])
    return result


def draft_profile(probe: dict[str, Any]) -> dict[str, Any]:
    scheduler_name = probe["scheduler"]
    resources = probe.get("resources") or []
    drafted = []
    for index, item in enumerate(resources):
        capabilities = item.get("capabilities", {})
        gpus = int(capabilities.get("gpus_per_node", 0))
        gpu_type = capabilities.get("gpu_type")
        drafted_capabilities = {
            "cpus_per_node": int(capabilities.get("cpus_per_node", 1)),
            "gpus_per_node": gpus,
            "memory_mb_per_node": int(capabilities.get("memory_mb_per_node", 4096)),
        }
        if gpu_type:
            drafted_capabilities["gpu_type"] = str(gpu_type)
        drafted.append(
            {
                "name": item.get("name", f"resource-{index + 1}"),
                "classes": ["light", "data", "gpu", "md"]
                if gpus
                else ["light", "data"],
                "capabilities": drafted_capabilities,
                "limits": {"max_nodes": 1, "max_walltime_seconds": 86400},
                "scheduler_options": dict(item.get("scheduler_options", {})),
            }
        )
    if not drafted:
        drafted = [
            {
                "name": "configure-me",
                "capabilities": {
                    "cpus_per_node": 1,
                    "gpus_per_node": 0,
                    "memory_mb_per_node": 4096,
                },
                "limits": {"max_nodes": 1, "max_walltime_seconds": 3600},
                "scheduler_options": {},
                "incomplete": True,
            }
        ]
    return {
        "schema_version": SCHEMA_VERSION,
        "name": f"{scheduler_name}-cluster",
        "scheduler": scheduler_name,
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
        "resources": drafted,
    }


def save_profile(
    draft_path: Path, output: Path, approval: str, *, replace: bool = False
) -> dict[str, Any]:
    draft = read_json(draft_path)
    draft_id = fingerprint(draft)
    if approval != draft_id:
        raise ValueError("--approve must exactly match the profile draft hash")
    from .profile import validate_profile

    profile = validate_profile(draft)
    target = output.expanduser()
    if os.path.lexists(target) and not replace:
        raise FileExistsError(f"Cluster profile already exists: {target}")
    write_json(target, profile)
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "saved",
        "path": str(output.expanduser().resolve()),
        "profile_fingerprint": profile_fingerprint(profile),
    }
