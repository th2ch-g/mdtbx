import json

import pytest

from mdtbx.agent.model import fingerprint
from mdtbx.agent.runtime import build_plan, execute_plan, save_profile


def _batch_profile():
    return {
        "schema_version": 2,
        "name": "batch",
        "scheduler": "slurm",
        "work_root": ".mdtbx/runs",
        "environment": {"prologue": [], "env": {}},
        "policy": {
            "batch_classes": ["md"],
            "memory_safety_factor": 1.0,
            "walltime_safety_factor": 1.0,
            "minimum_confidence": 0.7,
            "pilot_walltime_seconds": 600,
        },
        "resources": [
            {
                "name": "md-node",
                "classes": ["md"],
                "capabilities": {
                    "cpus_per_node": 8,
                    "gpus_per_node": 1,
                    "memory_mb_per_node": 16384,
                },
                "limits": {
                    "max_nodes": 1,
                    "max_walltime_seconds": 86400,
                },
                "scheduler_options": {},
            }
        ],
    }


def test_unbenchmarked_md_requires_pilot(tmp_path):
    profile = tmp_path / "cluster.json"
    profile.write_text(json.dumps(_batch_profile()))
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "cluster_profile": str(profile),
        "steps": [
            {
                "id": "fep",
                "command": "run_fep",
                "arguments": {"path": "fep"},
                "depends_on": [],
            }
        ],
    }
    plan = build_plan(request)
    assert plan["plan_kind"] == "pilot"
    assert plan["steps"][0]["arguments"]["nsteps"] == 1000
    assert plan["steps"][0]["confidence"] < 0.7


def test_explicit_pilot_walltime_is_preserved(tmp_path):
    profile = tmp_path / "cluster.json"
    profile.write_text(json.dumps(_batch_profile()))
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "cluster_profile": str(profile),
        "steps": [
            {
                "id": "fep",
                "command": "run_fep",
                "arguments": {"path": "fep", "nsteps": 50000},
                "depends_on": [],
                "resources": {
                    "resource": "md-node",
                    "nodes": 1,
                    "cpus_per_node": 8,
                    "tasks_per_node": 8,
                    "gpus_per_node": 1,
                    "memory_mb": 8192,
                    "walltime_seconds": 3600,
                },
            }
        ],
    }

    plan = build_plan(request)

    assert plan["plan_kind"] == "pilot"
    assert plan["steps"][0]["resources"]["walltime_seconds"] == 3600
    assert plan["steps"][0]["arguments"]["nsteps"] == 50000


def test_approved_pilot_evidence_enables_production_plan(tmp_path):
    profile = tmp_path / "cluster.json"
    profile.write_text(json.dumps(_batch_profile()))
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "cluster_profile": str(profile),
        "steps": [
            {
                "id": "fep",
                "command": "run_fep",
                "arguments": {"path": "fep", "nsteps": 5000},
                "depends_on": [],
                "evidence": [
                    {
                        "source": "approved_pilot",
                        "run_id": "pilot-1",
                    }
                ],
            }
        ],
    }
    plan = build_plan(request)
    assert plan["plan_kind"] == "production"
    assert plan["steps"][0]["confidence"] >= 0.7


def test_existing_output_requires_overwrite_approval(tmp_path):
    output = tmp_path / "converted.pdb"
    output.write_text("existing\n")
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "steps": [
            {
                "id": "convert",
                "command": "convert",
                "arguments": {
                    "structure": "input.pdb",
                    "output": str(output),
                },
                "depends_on": [],
            }
        ],
    }
    plan = build_plan(request)
    assert plan["approval"]["overwrite"] is True
    assert plan["steps"][0]["existing_artifacts"] == [str(output)]
    with pytest.raises(ValueError, match="approve-overwrite"):
        execute_plan(plan, approval=plan["plan_id"])


def test_profile_save_refuses_unapproved_replacement(tmp_path):
    draft = tmp_path / "draft.json"
    output = tmp_path / "cluster.json"
    profile = _batch_profile()
    draft.write_text(json.dumps(profile))
    output.write_text("{}\n")
    approval = fingerprint(profile)

    with pytest.raises(FileExistsError):
        save_profile(draft, output, approval)
    result = save_profile(draft, output, approval, replace=True)
    assert result["status"] == "saved"
    assert json.loads(output.read_text())["name"] == "batch"


def test_profile_classes_do_not_fall_back_to_wrong_node(tmp_path):
    profile_data = _batch_profile()
    profile_data["resources"][0]["classes"] = ["data"]
    profile = tmp_path / "cluster.json"
    profile.write_text(json.dumps(profile_data))
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "cluster_profile": str(profile),
        "steps": [
            {
                "id": "fep",
                "command": "run_fep",
                "arguments": {"path": "fep"},
                "depends_on": [],
            }
        ],
    }
    with pytest.raises(ValueError, match="No cluster resource"):
        build_plan(request)


def test_step_id_cannot_escape_run_directory(tmp_path):
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "steps": [
            {
                "id": "../escape",
                "command": "show_npy",
                "arguments": {"npy": "values.npy"},
            }
        ],
    }
    with pytest.raises(ValueError, match="Invalid step id"):
        build_plan(request)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("nodes", 0),
        ("cpus_per_node", -1),
        ("gpus_per_node", -1),
        ("memory_mb", 0),
        ("walltime_seconds", -10),
    ],
)
def test_requested_resources_must_be_nonnegative(field, value, tmp_path):
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "steps": [
            {
                "id": "inspect",
                "command": "show_npy",
                "arguments": {"npy": "values.npy"},
                "resources": {field: value},
            }
        ],
    }
    with pytest.raises(ValueError, match=field):
        build_plan(request)


def test_per_step_execution_mode_is_rejected(tmp_path):
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "steps": [
            {
                "id": "inspect",
                "command": "show_npy",
                "arguments": {"npy": "values.npy"},
                "execution": "local",
            }
        ],
    }
    with pytest.raises(ValueError, match="Unknown step fields"):
        build_plan(request)


def test_tasks_per_node_cannot_exceed_requested_cpus(tmp_path):
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "steps": [
            {
                "id": "inspect",
                "command": "show_npy",
                "arguments": {"npy": "values.npy"},
                "resources": {"cpus_per_node": 4, "tasks_per_node": 8},
            }
        ],
    }
    with pytest.raises(ValueError, match="tasks_per_node"):
        build_plan(request)


def test_artifact_appearance_after_planning_requires_new_plan(tmp_path):
    output = tmp_path / "converted.pdb"
    request = {
        "schema_version": 2,
        "cwd": str(tmp_path),
        "steps": [
            {
                "id": "convert",
                "command": "convert",
                "arguments": {"structure": "input.pdb", "output": str(output)},
            }
        ],
    }
    plan = build_plan(request)
    output.write_text("appeared after planning\n")
    with pytest.raises(ValueError, match="Filesystem changed"):
        execute_plan(plan, approval=plan["plan_id"])


def test_directories_and_broken_symlinks_require_overwrite_approval(tmp_path):
    directory = tmp_path / "directory-output"
    directory.mkdir()
    directory_plan = build_plan(
        {
            "schema_version": 2,
            "cwd": str(tmp_path),
            "steps": [
                {
                    "id": "convert-dir",
                    "command": "convert",
                    "arguments": {
                        "structure": "input.pdb",
                        "output": str(directory),
                    },
                }
            ],
        }
    )
    assert directory_plan["steps"][0]["artifact_snapshots"][0]["kind"] == "directory"
    assert directory_plan["approval"]["overwrite"] is True

    symlink = tmp_path / "broken-output"
    symlink.symlink_to(tmp_path / "missing-target")
    symlink_plan = build_plan(
        {
            "schema_version": 2,
            "cwd": str(tmp_path),
            "steps": [
                {
                    "id": "convert-link",
                    "command": "convert",
                    "arguments": {
                        "structure": "input.pdb",
                        "output": str(symlink),
                    },
                }
            ],
        }
    )
    assert symlink_plan["steps"][0]["artifact_snapshots"][0]["kind"] == "symlink"
    assert symlink_plan["approval"]["overwrite"] is True


def test_partial_batch_submission_persists_issued_job_ids(monkeypatch, tmp_path):
    from mdtbx.agent import runtime
    from mdtbx.agent.model import read_json

    profile_data = _batch_profile()
    profile_data["resources"][0]["classes"] = ["light"]
    profile = tmp_path / "cluster.json"
    profile.write_text(json.dumps(profile_data))
    plan = build_plan(
        {
            "schema_version": 2,
            "cwd": str(tmp_path),
            "cluster_profile": str(profile),
            "execution": "batch",
            "steps": [
                {
                    "id": "first",
                    "command": "show_npy",
                    "arguments": {"npy": "first.npy"},
                },
                {
                    "id": "second",
                    "command": "show_npy",
                    "arguments": {"npy": "second.npy"},
                    "depends_on": ["first"],
                },
            ],
        }
    )

    class FakeScheduler:
        def __init__(self):
            self.calls = 0

        def directives(self, step, resource, profile):
            return ["# scheduler"]

        def submit(self, script, *, dependencies, resource):
            self.calls += 1
            if self.calls == 1:
                return "job-101"
            raise RuntimeError("submission failed")

    adapter = FakeScheduler()
    monkeypatch.setattr(runtime, "scheduler", lambda _name: adapter)
    with pytest.raises(RuntimeError, match="submission failed"):
        execute_plan(plan, approval=plan["plan_id"])

    state_paths = list((tmp_path / ".mdtbx" / "runs").glob("*/state.json"))
    assert len(state_paths) == 1
    state = read_json(state_paths[0])
    assert state["status"] == "submission_failed"
    assert state["steps"]["first"]["job_id"] == "job-101"
    assert state["steps"]["second"]["state"] == "failed"
