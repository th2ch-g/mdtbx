"""Formal JSON Schema documents for the agent protocol."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from .model import SCHEMA_VERSION

_STEP_ID = r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$"
_POSITIVE_INTEGER = {"type": "integer", "minimum": 1}

REQUEST_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://mdtbx.dev/schemas/agent-request-v2.json",
    "title": "mdtbx agent request",
    "type": "object",
    "additionalProperties": False,
    "required": ["schema_version", "steps"],
    "properties": {
        "schema_version": {"const": SCHEMA_VERSION},
        "name": {"type": "string", "minLength": 1},
        "cwd": {"type": "string", "minLength": 1},
        "cluster_profile": {"type": ["string", "null"]},
        "execution": {"enum": ["auto", "local", "batch"]},
        "steps": {
            "type": "array",
            "minItems": 1,
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": ["id", "command"],
                "properties": {
                    "id": {"type": "string", "pattern": _STEP_ID},
                    "command": {"type": "string", "minLength": 1},
                    "arguments": {"type": "object"},
                    "depends_on": {
                        "type": "array",
                        "uniqueItems": True,
                        "items": {"type": "string", "pattern": _STEP_ID},
                    },
                    "resources": {
                        "type": "object",
                        "additionalProperties": False,
                        "properties": {
                            "resource": {"type": "string", "minLength": 1},
                            "nodes": _POSITIVE_INTEGER,
                            "cpus_per_node": _POSITIVE_INTEGER,
                            "tasks_per_node": _POSITIVE_INTEGER,
                            "gpus_per_node": {"type": "integer", "minimum": 0},
                            "memory_mb": _POSITIVE_INTEGER,
                            "walltime_seconds": _POSITIVE_INTEGER,
                        },
                    },
                    "evidence": {
                        "type": "array",
                        "items": {"type": "object"},
                    },
                    "confidence": {
                        "type": "number",
                        "minimum": 0,
                        "maximum": 1,
                    },
                },
            },
        },
    },
}

PROFILE_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://mdtbx.dev/schemas/cluster-profile-v2.json",
    "title": "mdtbx cluster profile",
    "type": "object",
    "additionalProperties": False,
    "required": ["schema_version", "name", "scheduler", "resources"],
    "properties": {
        "schema_version": {"const": SCHEMA_VERSION},
        "name": {"type": "string", "minLength": 1},
        "scheduler": {"enum": ["local", "slurm", "age", "pjm"]},
        "work_root": {"type": "string", "minLength": 1},
        "runner_argv": {
            "type": "array",
            "minItems": 1,
            "items": {"type": "string", "minLength": 1},
        },
        "environment": {"type": "object"},
        "policy": {"type": "object"},
        "resources": {"type": "array", "minItems": 1},
    },
}

PLAN_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://mdtbx.dev/schemas/agent-plan-v2.json",
    "title": "mdtbx immutable agent plan",
    "type": "object",
    "required": [
        "schema_version",
        "plan_id",
        "cwd",
        "execution",
        "steps",
        "approval",
    ],
    "properties": {
        "schema_version": {"const": SCHEMA_VERSION},
        "plan_id": {"type": "string", "pattern": "^[0-9a-f]{64}$"},
        "execution": {"enum": ["local", "batch"]},
        "steps": {"type": "array", "minItems": 1},
    },
}

STATE_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://mdtbx.dev/schemas/agent-state-v2.json",
    "title": "mdtbx durable run state",
    "type": "object",
    "required": ["schema_version", "run_id", "plan_id", "status", "steps"],
    "properties": {"schema_version": {"const": SCHEMA_VERSION}},
}

RESULT_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://mdtbx.dev/schemas/agent-result-v2.json",
    "title": "mdtbx collected run result",
    "type": "object",
    "required": ["schema_version", "run_id", "plan_id", "status", "steps"],
    "properties": {"schema_version": {"const": SCHEMA_VERSION}},
}


def schemas() -> dict[str, dict[str, Any]]:
    """Return isolated schema documents safe for JSON serialization."""
    return deepcopy(
        {
            "request": REQUEST_SCHEMA,
            "profile": PROFILE_SCHEMA,
            "plan": PLAN_SCHEMA,
            "state": STATE_SCHEMA,
            "result": RESULT_SCHEMA,
        }
    )
