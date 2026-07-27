import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from mdtbx.agent.profile import validate_profile
from mdtbx.agent.schedulers import AgeScheduler, PjmScheduler, SlurmScheduler
from mdtbx.cli import _run_agent_mode


def _completed(stdout):
    return subprocess.CompletedProcess([], 0, stdout=stdout, stderr="")


def test_json_mode_wraps_system_exit(capsys):
    def stop(_args):
        raise SystemExit(7)

    args = SimpleNamespace(
        _command="test",
        _agent_dry_run=False,
        func=stop,
    )
    assert _run_agent_mode(args) == 7
    payload = json.loads(capsys.readouterr().out)
    assert payload["ok"] is False
    assert payload["exit_code"] == 7
    assert payload["error"]["type"] == "SystemExit"


@pytest.mark.parametrize("key", ["cpus_per_node", "memory_mb_per_node"])
def test_profile_requires_positive_capacity(key):
    profile = {
        "schema_version": 2,
        "name": "invalid",
        "scheduler": "slurm",
        "resources": [
            {
                "name": "node",
                "capabilities": {
                    "cpus_per_node": 1,
                    "gpus_per_node": 0,
                    "memory_mb_per_node": 1,
                },
                "scheduler_options": {},
            }
        ],
    }
    profile["resources"][0]["capabilities"][key] = 0
    with pytest.raises(ValueError, match="must be positive"):
        validate_profile(profile)


def test_slurm_submit_uses_afterok(monkeypatch, tmp_path):
    calls = []

    def fake_run(argv, check=False):
        calls.append((argv, check))
        return _completed("123;cluster\n")

    monkeypatch.setattr("mdtbx.agent.schedulers._run", fake_run)
    job_id = SlurmScheduler().submit(
        tmp_path / "job.sh",
        dependencies=["100", "101"],
        resource={},
    )
    assert job_id == "123"
    assert calls[0][0][2] == "--dependency=afterok:100:101"


def test_age_submit_uses_hold_jid(monkeypatch, tmp_path):
    calls = []

    def fake_run(argv, check=False):
        calls.append((argv, check))
        return _completed("456.cluster\n")

    monkeypatch.setattr("mdtbx.agent.schedulers._run", fake_run)
    job_id = AgeScheduler().submit(
        tmp_path / "job.sh",
        dependencies=["200", "201"],
        resource={"scheduler_options": {}},
    )
    assert job_id == "456"
    assert calls[0][0][2:4] == ["-hold_jid", "200,201"]


def test_pjm_submit_requires_profile_dependency_syntax(tmp_path):
    with pytest.raises(ValueError, match="dependency_args"):
        PjmScheduler().submit(
            Path(tmp_path) / "job.sh",
            dependencies=["300"],
            resource={"scheduler_options": {}},
        )


def test_pjm_submit_parses_job_id(monkeypatch, tmp_path):
    calls = []

    def fake_run(argv, check=False):
        calls.append((argv, check))
        return _completed("[INFO] PJM 0000 pjsub Job 789 submitted.\n")

    monkeypatch.setattr("mdtbx.agent.schedulers._run", fake_run)
    job_id = PjmScheduler().submit(
        tmp_path / "job.sh",
        dependencies=["300"],
        resource={
            "scheduler_options": {
                "dependency_args": ["--sparam", "jid={job_ids}"],
            }
        },
    )
    assert job_id == "789"
    assert calls[0][0][1:3] == ["--sparam", "jid=300"]
