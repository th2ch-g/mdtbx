import json
import subprocess
import sys
from argparse import _SubParsersAction
from pathlib import Path

import pytest

from mdtbx.agent.model import finalize_plan, fingerprint, verify_plan
from mdtbx.agent.profile import validate_profile
from mdtbx.agent.runtime import build_plan, collect_run, load_run
from mdtbx.agent.schedulers import AgeScheduler, PjmScheduler, SlurmScheduler
from mdtbx.cli import create_parser


def _profile(scheduler="slurm"):
    return {
        "schema_version": 2,
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
        "agent_cancel",
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
    plan = finalize_plan({"schema_version": 2, "steps": []})
    assert verify_plan(plan) == plan["plan_id"]
    plan["steps"].append({"id": "changed"})
    with pytest.raises(ValueError, match="does not match"):
        verify_plan(plan)


def test_resource_plan_uses_external_profile(tmp_path):
    profile_path = tmp_path / "cluster.json"
    profile_path.write_text(json.dumps(_profile()))
    request = {
        "schema_version": 2,
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


def test_slurm_scheduler_maps_node_cpus_to_tasks():
    profile = _profile("slurm")
    resource = profile["resources"][0]
    step = {
        "id": "hrex",
        "resources": {
            "nodes": 1,
            "cpus_per_node": 16,
            "tasks_per_node": 16,
            "gpus_per_node": 1,
            "memory_mb": 8192,
            "walltime_seconds": 600,
        },
    }

    lines = SlurmScheduler().directives(step, resource, profile)

    assert "#SBATCH --ntasks-per-node=16" in lines
    assert "#SBATCH --cpus-per-task=1" in lines


def test_slurm_scheduler_cancels_exact_job(monkeypatch):
    calls = []

    def fake_run(argv, *, check=False):
        calls.append((argv, check))

    monkeypatch.setattr("mdtbx.agent.schedulers._run", fake_run)

    SlurmScheduler().cancel("12345")

    assert calls == [(["scancel", "12345"], True)]


def test_slurm_scheduler_writes_logs_next_to_job_script(monkeypatch, tmp_path):
    calls = []

    def fake_run(argv, *, check=False):
        calls.append((argv, check))
        return subprocess.CompletedProcess(argv, 0, stdout="12345\n", stderr="")

    monkeypatch.setattr("mdtbx.agent.schedulers._run", fake_run)
    script = tmp_path / "step" / "job.sh"
    script.parent.mkdir()
    script.touch()

    job_id = SlurmScheduler().submit(script, dependencies=["100"], resource={})

    assert job_id == "12345"
    assert calls == [
        (
            [
                "sbatch",
                "--parsable",
                "--dependency=afterok:100",
                f"--output={script.parent / 'scheduler.stdout.log'}",
                f"--error={script.parent / 'scheduler.stderr.log'}",
                str(script),
            ],
            True,
        )
    ]


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


def test_agent_schema_exposes_formal_v2_documents():
    from mdtbx.agent.runtime import schema

    document = schema("show_npy")
    assert document["schema_version"] == 2
    assert set(document["schemas"]) == {"request", "profile", "plan", "state", "result"}
    assert document["schemas"]["request"]["additionalProperties"] is False
    assert document["schemas"]["request"]["properties"]["schema_version"] == {
        "const": 2
    }


def test_run_fep_path_is_tracked_as_input_and_output():
    from mdtbx.cli import create_parser
    from mdtbx.agent.registry import all_descriptors

    item = all_descriptors(create_parser({"run_fep"}))["run_fep"]

    assert "path" in item["inputs"]
    assert "path" in item["outputs"]


def test_fep_schedule_paths_are_tracked_by_agent_schema():
    from mdtbx.agent.registry import all_descriptors

    parser = create_parser({"setup_fep", "setup_abfe", "optimize_fep_schedule"})
    items = all_descriptors(parser)

    assert "schedule" in items["setup_fep"]["inputs"]
    assert "anchor_trajectory" in items["setup_abfe"]["inputs"]
    assert "solvent_charge_schedule" in items["setup_abfe"]["inputs"]
    assert "log" in items["optimize_fep_schedule"]["inputs"]
    assert "output" in items["optimize_fep_schedule"]["outputs"]
    assert "path" in items["optimize_fep_schedule"]["outputs"]


def test_load_run_finds_default_root_from_nested_directory(tmp_path, monkeypatch):
    run_id = "20260101T000000000000Z-example"
    run_dir = tmp_path / ".mdtbx" / "runs" / run_id
    run_dir.mkdir(parents=True)
    (run_dir / "state.json").write_text(
        json.dumps({"schema_version": 2, "run_id": run_id})
    )
    nested = tmp_path / "project" / "workflow"
    nested.mkdir(parents=True)
    monkeypatch.chdir(nested)

    resolved, state = load_run(run_id)

    assert resolved == run_dir.resolve()
    assert state["run_id"] == run_id


def test_collect_run_persists_scheduler_step_state(tmp_path, monkeypatch):
    run_dir = tmp_path / "run-id"
    run_dir.mkdir()
    state = {
        "schema_version": 2,
        "run_id": "run-id",
        "plan_id": "plan-id",
        "scheduler": "slurm",
        "steps": {"simulation": {"state": "queued", "job_id": "123"}},
    }
    (run_dir / "state.json").write_text(json.dumps(state))
    (run_dir / "plan.json").write_text(
        json.dumps({"steps": [{"id": "simulation", "artifacts": []}]})
    )

    class FakeScheduler:
        def status(self, job_id):
            return {"job_id": job_id, "state": "succeeded"}

    monkeypatch.setattr("mdtbx.agent.runtime.scheduler", lambda _name: FakeScheduler())

    result = collect_run(str(run_dir))
    persisted = json.loads((run_dir / "state.json").read_text())

    assert result["status"] == "succeeded"
    assert persisted["steps"]["simulation"]["state"] == "succeeded"
