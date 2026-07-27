"""Planning, resource selection, execution, and run-state collection."""

from __future__ import annotations

import json
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

from .model import (
    SCHEMA_VERSION,
    finalize_plan,
    fingerprint,
    json_value,
    read_json,
    utc_now,
    verify_plan,
    write_json,
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


def _root_parser():
    from ..cli import create_parser

    return create_parser()


def schema(command: str | None = None) -> dict[str, Any]:
    commands = all_descriptors(_root_parser())
    if command is not None:
        try:
            commands = {command: commands[command]}
        except KeyError as error:
            raise ValueError(f"Unknown mdtbx command: {command}") from error
    return {
        "schema_version": SCHEMA_VERSION,
        "commands": commands,
    }


def _validated_steps(request: dict[str, Any]) -> list[dict[str, Any]]:
    raw_steps = request.get("steps")
    if not isinstance(raw_steps, list) or not raw_steps:
        raise ValueError("Request requires a non-empty steps list")
    root = _root_parser()
    descriptors = all_descriptors(root)
    steps = []
    identifiers = set()
    for raw in raw_steps:
        if not isinstance(raw, dict):
            raise ValueError("Each request step must be an object")
        identifier = raw.get("id")
        command = raw.get("command")
        arguments = raw.get("arguments", {})
        if not isinstance(identifier, str) or not identifier:
            raise ValueError("Each request step requires an id")
        if identifier in identifiers:
            raise ValueError(f"Duplicate step id: {identifier}")
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
                not isinstance(requested_confidence, (int, float))
                or isinstance(requested_confidence, bool)
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
            isinstance(item, str) for item in depends_on
        ):
            raise ValueError(f"{identifier}.depends_on must be a string list")
        descriptor = descriptors[command]
        steps.append(
            {
                "id": identifier,
                "command": command,
                "arguments": normalized,
                "argv": argv,
                "depends_on": depends_on,
                "descriptor": descriptor,
                "requested_resources": raw.get("resources"),
                "requested_execution": raw.get("execution", "auto"),
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
    completed = set()
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


def _input_bytes(step: dict[str, Any]) -> int:
    total = 0
    for name in step["descriptor"]["inputs"]:
        value = step["arguments"].get(name)
        values = value if isinstance(value, list) else [value]
        for item in values:
            if not isinstance(item, str):
                continue
            path = Path(item).expanduser()
            if path.is_file():
                total += path.stat().st_size
    return total


def _policy(profile: dict[str, Any], key: str, default: Any) -> Any:
    return profile.get("policy", {}).get(key, default)


def _resource_candidates(
    profile: dict[str, Any], resource_class: str
) -> list[dict[str, Any]]:
    candidates = []
    unclassified = []
    for resource in profile["resources"]:
        classes = resource.get("classes")
        if not classes:
            unclassified.append(resource)
        if classes and resource_class not in classes:
            continue
        candidates.append(resource)
    return candidates or unclassified


def _select_resources(
    step: dict[str, Any], profile: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any], float, list[dict[str, Any]]]:
    requested = step.get("requested_resources")
    resource_class = step["descriptor"]["resource_class"]
    input_bytes = _input_bytes(step)
    evidence = [
        {
            "source": "command_registry",
            "resource_class": resource_class,
        },
        {
            "source": "input_files",
            "total_bytes": input_bytes,
        },
    ]
    evidence.extend(json_value(step["resource_evidence"]))

    explicit = isinstance(requested, dict)
    memory_estimate = DEFAULT_MEMORY_MB[resource_class]
    if resource_class == "data":
        memory_estimate = max(memory_estimate, int(input_bytes * 3 / 1_000_000) + 1024)
    memory_factor = float(_policy(profile, "memory_safety_factor", 1.25))
    walltime_factor = float(_policy(profile, "walltime_safety_factor", 1.5))
    memory_mb = int(memory_estimate * memory_factor)
    walltime = int(DEFAULT_WALLTIME_SECONDS[resource_class] * walltime_factor)
    nodes = 1
    minimum_gpus = 1 if resource_class == "gpu" else 0
    minimum_cpus = 1
    resource_name = None
    if explicit:
        nodes = int(requested.get("nodes", nodes))
        memory_mb = int(requested.get("memory_mb", memory_mb))
        walltime = int(requested.get("walltime_seconds", walltime))
        minimum_gpus = int(requested.get("gpus_per_node", minimum_gpus))
        minimum_cpus = int(requested.get("cpus_per_node", minimum_cpus))
        resource_name = requested.get("resource")
        evidence.append({"source": "request", "resources": json_value(requested)})

    candidates = _resource_candidates(profile, resource_class)
    suitable = []
    for candidate in candidates:
        if resource_name and candidate["name"] != resource_name:
            continue
        capabilities = candidate.get("capabilities", {})
        if capabilities.get("memory_mb_per_node", 0) < memory_mb:
            continue
        if capabilities.get("gpus_per_node", 0) < minimum_gpus:
            continue
        if capabilities.get("cpus_per_node", 0) < minimum_cpus:
            continue
        limits = candidate.get("limits", {})
        if nodes > limits.get("max_nodes", nodes):
            continue
        if walltime > limits.get("max_walltime_seconds", walltime):
            continue
        score = (
            capabilities.get("gpus_per_node", 0) * 10**12
            + capabilities.get("memory_mb_per_node", 0) * 10**4
            + capabilities.get("cpus_per_node", 0)
        )
        suitable.append((score, candidate))
    if not suitable:
        raise ValueError(f"No cluster resource can satisfy step {step['id']}")
    resource = min(suitable, key=lambda item: item[0])[1]
    capabilities = resource.get("capabilities", {})
    allocation = {
        "resource": resource["name"],
        "nodes": nodes,
        "cpus_per_node": min(
            int(capabilities.get("cpus_per_node", minimum_cpus)),
            (
                max(minimum_cpus, int(requested.get("cpus_per_node", minimum_cpus)))
                if explicit
                else minimum_cpus
            ),
        ),
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
            "capabilities": capabilities,
            "safety_factors": {
                "memory": memory_factor,
                "walltime": walltime_factor,
            },
        }
    )
    return allocation, resource, confidence, evidence


def build_plan(request: dict[str, Any]) -> dict[str, Any]:
    if request.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported request schema")
    cwd = Path(request.get("cwd", ".")).expanduser().resolve()
    if not cwd.is_dir():
        raise FileNotFoundError(f"Working directory not found: {cwd}")
    profile_value = request.get("cluster_profile")
    profile_path, profile = load_profile(profile_value)
    steps = _validated_steps(request)
    batch_classes = set(_policy(profile, "batch_classes", ["md", "gpu", "quantum"]))
    any_batch = any(
        step["requested_execution"] == "batch"
        or (
            step["requested_execution"] == "auto"
            and step["descriptor"]["resource_class"] in batch_classes
        )
        for step in steps
    )
    if any_batch and profile["scheduler"] == "local":
        raise ValueError("Batch execution requires a non-local cluster profile")

    planned_steps = []
    plan_kind = "production"
    minimum_confidence = float(_policy(profile, "minimum_confidence", 0.7))
    for step in steps:
        allocation, resource, confidence, evidence = _select_resources(step, profile)
        execution = "batch" if any_batch else "local"
        arguments = dict(step["arguments"])
        argv = list(step["argv"])
        if confidence < minimum_confidence and step["descriptor"]["pilot_capable"]:
            plan_kind = "pilot"
            allocation["walltime_seconds"] = min(
                allocation["walltime_seconds"],
                int(_policy(profile, "pilot_walltime_seconds", 600)),
            )
            if (
                step["command"] in {"run_fep", "run_abfe"}
                and arguments.get("nsteps") is None
            ):
                arguments["nsteps"] = 1000
                parser = command_parser(_root_parser(), step["command"])
                argv = arguments_to_argv(parser, arguments)
                arguments = normalized_arguments(parser, argv)
                evidence.append(
                    {
                        "source": "pilot_policy",
                        "override": {"nsteps": 1000},
                    }
                )
        inputs = input_paths(arguments, step["descriptor"])
        outputs = artifact_paths(arguments, step["descriptor"])
        resolved_outputs = []
        existing_artifacts = []
        for value in outputs:
            target = Path(value)
            if not target.is_absolute():
                target = cwd / target
            target = target.resolve()
            resolved_outputs.append(str(target))
            if target.is_file() or target.is_symlink():
                existing_artifacts.append(str(target))
        destructive_targets = []
        if step["descriptor"]["risk"] == "destructive":
            destructive_targets = inputs
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
                "artifacts": resolved_outputs,
                "existing_artifacts": existing_artifacts,
                "destructive_targets": destructive_targets,
            }
        )

    profile_hash = profile_fingerprint(profile)
    result = {
        "schema_version": SCHEMA_VERSION,
        "created_at": utc_now(),
        "name": request.get("name", "mdtbx-agent-plan"),
        "plan_kind": plan_kind,
        "cwd": str(cwd),
        "cluster_profile": str(profile_path) if profile_path else None,
        "profile_fingerprint": profile_hash,
        "scheduler": profile["scheduler"],
        "execution": "batch" if any_batch else "local",
        "runner_argv": profile.get("runner_argv", [sys.executable, "-m", "src"]),
        "work_root": profile.get("work_root", ".mdtbx/runs"),
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
    run_id = (
        f"{utc_now().replace(':', '').replace('-', '')[:15]}-{plan['plan_id'][:12]}"
    )
    return root / run_id


def _append_event(run_dir: Path, event: dict[str, Any]) -> None:
    path = run_dir / "events.jsonl"
    with path.open("a") as handle:
        handle.write(json.dumps(json_value(event), ensure_ascii=False) + "\n")


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
    command = [
        *plan["runner_argv"],
        step["command"],
        *step["argv"],
        "--json",
    ]
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
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o700)
    return script


def execute_plan(
    plan: dict[str, Any],
    *,
    approval: str,
    approve_unsafe: bool = False,
    approve_destructive: bool = False,
    approve_overwrite: bool = False,
) -> dict[str, Any]:
    plan_id = verify_plan(plan)
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

    run_dir = _run_root(plan)
    run_dir.mkdir(parents=True)
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
        "steps": {},
    }
    write_json(run_dir / "state.json", state)
    _append_event(
        run_dir,
        {"at": utc_now(), "event": "approved", "plan_id": plan_id},
    )

    ordered = _topological(plan["steps"])
    if plan["execution"] == "local":
        for step in ordered:
            step_dir = run_dir / "steps" / step["id"]
            step_dir.mkdir(parents=True)
            command = [
                *plan["runner_argv"],
                step["command"],
                *step["argv"],
                "--json",
            ]
            result = subprocess.run(
                command,
                cwd=plan["cwd"],
                capture_output=True,
                text=True,
            )
            (step_dir / "result.json").write_text(result.stdout)
            (step_dir / "stderr.log").write_text(result.stderr)
            status = "succeeded" if result.returncode == 0 else "failed"
            state["steps"][step["id"]] = {
                "state": status,
                "returncode": result.returncode,
            }
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
                break
        else:
            state["status"] = "succeeded"
    else:
        adapter = scheduler(plan["scheduler"])
        submitted: dict[str, str] = {}
        for step in ordered:
            dependencies = [submitted[item] for item in step["depends_on"]]
            script = _render_script(plan, step, run_dir, profile)
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
            }
            _append_event(
                run_dir,
                {
                    "at": utc_now(),
                    "event": "job_submitted",
                    "step": step["id"],
                    "job_id": job_id,
                },
            )
        state["status"] = "submitted"
    write_json(run_dir / "state.json", state)
    return state


def load_run(value: str) -> tuple[Path, dict[str, Any]]:
    path = Path(value).expanduser()
    if path.is_dir():
        run_dir = path.resolve()
    else:
        matches = list(Path(".mdtbx/runs").glob("*/state.json"))
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
    return run_dir, read_json(state_path)


def run_status(value: str) -> dict[str, Any]:
    run_dir, state = load_run(value)
    if state["scheduler"] == "local":
        return state
    adapter = scheduler(state["scheduler"])
    statuses = {}
    for step, item in state["steps"].items():
        statuses[step] = adapter.status(item["job_id"])
    values = {item["state"] for item in statuses.values()}
    if values and values <= {"succeeded"}:
        overall = "succeeded"
    elif "failed" in values or "cancelled" in values:
        overall = "failed"
    elif "running" in values:
        overall = "running"
    elif "queued" in values:
        overall = "queued"
    else:
        overall = "unknown"
    return {
        **state,
        "run_dir": str(run_dir),
        "status": overall,
        "checked_at": utc_now(),
        "scheduler_status": statuses,
    }


def collect_run(value: str) -> dict[str, Any]:
    run_dir, state = load_run(value)
    status = run_status(str(run_dir))
    results = {}
    for step in state["steps"]:
        result_path = run_dir / "steps" / step / "result.json"
        if not result_path.is_file() or not result_path.read_text().strip():
            results[step] = None
            continue
        try:
            results[step] = read_json(result_path)
        except json.JSONDecodeError:
            results[step] = {
                "status": "invalid",
                "path": str(result_path),
            }
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
    state["status"] = status["status"]
    state["collected_at"] = result["collected_at"]
    write_json(run_dir / "state.json", state)
    _append_event(
        run_dir,
        {"at": utc_now(), "event": "collected", "status": result["status"]},
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
        drafted.append(
            {
                "name": item.get("name", f"resource-{index + 1}"),
                "capabilities": {
                    "cpus_per_node": int(capabilities.get("cpus_per_node", 1)),
                    "gpus_per_node": int(capabilities.get("gpus_per_node", 0)),
                    "memory_mb_per_node": int(
                        capabilities.get("memory_mb_per_node", 4096)
                    ),
                },
                "limits": {
                    "max_nodes": 1,
                    "max_walltime_seconds": 86400,
                },
                "scheduler_options": {},
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
                "limits": {
                    "max_nodes": 1,
                    "max_walltime_seconds": 3600,
                },
                "scheduler_options": {},
                "incomplete": True,
            }
        ]
    return {
        "schema_version": SCHEMA_VERSION,
        "name": f"{scheduler_name}-cluster",
        "scheduler": scheduler_name,
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
    if target.exists() and not replace:
        raise FileExistsError(f"Cluster profile already exists: {target}")
    write_json(target, profile)
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "saved",
        "path": str(output.expanduser().resolve()),
        "profile_fingerprint": profile_fingerprint(profile),
    }
