import json
import subprocess
import sys
from argparse import _SubParsersAction
from pathlib import Path

import pytest

from src.agent.model import finalize_plan, fingerprint, verify_plan
from src.agent.profile import validate_profile
from src.agent.runtime import build_plan
from src.agent.schedulers import AgeScheduler, PjmScheduler, SlurmScheduler
from src.cli import create_parser


def _profile(scheduler="slurm"):
    return {
        "schema_version": 1,
        "name": "test-cluster",
        "scheduler": scheduler,
        "work_root": ".mdtbx/runs",
        "environment": {"prologue": [], "env": {}},
        "policy": {
            "batch_classes": ["data"],
            "memory_safety_factor": 1.0,
            "walltime_safety_factor": 1.0,
            "minimum_confidence": 0.7,
        },
        "resources": [
            {
                "name": "small",
                "classes": ["light", "data"],
                "capabilities": {
                    "cpus_per_node": 16,
                    "gpus_per_node": 0,
                    "memory_mb_per_node": 8192,
                },
                "limits": {
                    "max_nodes": 1,
                    "max_walltime_seconds": 86400,
                },
                "scheduler_options": {},
            }
        ],
    }


def test_agent_commands_and_common_flags_are_registered():
    parser = create_parser()
    subparsers = next(
        action for action in parser._actions if isinstance(action, _SubParsersAction)
    )
    assert {
        "agent_collect",
        "agent_plan",
        "agent_probe",
        "agent_profile_save",
        "agent_run",
        "agent_schema",
        "agent_status",
    } <= set(subparsers.choices)
    option_strings = {
        option
        for action in subparsers.choices["show_npy"]._actions
        for option in action.option_strings
    }
    assert {"--json", "--dry-run", "--cluster-profile"} <= option_strings


def test_plan_hash_rejects_changes():
    plan = finalize_plan({"schema_version": 1, "steps": []})
    assert verify_plan(plan) == plan["plan_id"]
    plan["steps"].append({"id": "changed"})
    with pytest.raises(ValueError, match="does not match"):
        verify_plan(plan)


def test_resource_plan_uses_external_profile(tmp_path):
    profile_path = tmp_path / "cluster.json"
    profile_path.write_text(json.dumps(_profile()))
    request = {
        "schema_version": 1,
        "cwd": str(tmp_path),
        "cluster_profile": str(profile_path),
        "steps": [
            {
                "id": "inspect",
                "command": "show_npy",
                "arguments": {"npy": "missing.npy"},
                "depends_on": [],
            }
        ],
    }
    plan = build_plan(request)
    assert plan["scheduler"] == "slurm"
    assert plan["steps"][0]["resource_profile"]["name"] == "small"
    assert plan["steps"][0]["execution"] == "local"
    assert verify_plan(plan) == plan["plan_id"]


def test_incomplete_profile_is_rejected():
    profile = _profile()
    profile["resources"][0]["incomplete"] = True
    with pytest.raises(ValueError, match="incomplete"):
        validate_profile(profile)


@pytest.mark.parametrize(
    ("adapter", "marker"),
    [
        (SlurmScheduler(), "#SBATCH"),
        (AgeScheduler(), "#$"),
        (PjmScheduler(), "#PJM"),
    ],
)
def test_scheduler_directives(adapter, marker):
    profile = _profile(adapter.name)
    resource = profile["resources"][0]
    step = {
        "id": "test",
        "resources": {
            "nodes": 1,
            "cpus_per_node": 1,
            "gpus_per_node": 0,
            "memory_mb": 4096,
            "walltime_seconds": 600,
        },
    }
    assert any(
        line.startswith(marker) for line in adapter.directives(step, resource, profile)
    )


def test_approval_hook_asks_and_denies():
    hook = Path("agent-skills/md-research/scripts/approval_hook.py")
    ask = subprocess.run(
        [sys.executable, str(hook)],
        input=json.dumps(
            {"tool_input": {"command": "pixi run mdtbx agent_run --plan p.json"}}
        ),
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(ask.stdout)["hookSpecificOutput"]["permissionDecision"] == "ask"
    deny = subprocess.run(
        [sys.executable, str(hook)],
        input=json.dumps({"tool_input": {"command": "sbatch job.sh"}}),
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(deny.stdout)["hookSpecificOutput"]["permissionDecision"] == "deny"


def test_profile_fingerprint_is_stable():
    profile = _profile()
    assert fingerprint(profile) == fingerprint(json.loads(json.dumps(profile)))
