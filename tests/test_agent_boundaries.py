import json

import pytest

from src.agent.model import fingerprint
from src.agent.runtime import build_plan, execute_plan, save_profile


def _batch_profile():
    return {
        "schema_version": 1,
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
        "schema_version": 1,
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


def test_approved_pilot_evidence_enables_production_plan(tmp_path):
    profile = tmp_path / "cluster.json"
    profile.write_text(json.dumps(_batch_profile()))
    request = {
        "schema_version": 1,
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
        "schema_version": 1,
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
        "schema_version": 1,
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
