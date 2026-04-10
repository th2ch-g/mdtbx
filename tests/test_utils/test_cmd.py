import subprocess
from argparse import Namespace
from pathlib import Path

import pytest

from src.utils import cmd


def test_cmd_run_uses_argv_list(monkeypatch):
    captured = {}

    def fake_run(command, check, cwd):
        captured["command"] = command
        captured["check"] = check
        captured["cwd"] = cwd

    monkeypatch.setattr(cmd.subprocess, "run", fake_run)

    cmd.run(Namespace(command=["python", "-c", "print('hello world')"]))

    assert isinstance(captured["command"], list)
    assert captured["command"][:3] == ["pixi", "run", "--manifest-path"]
    assert captured["command"][4:] == ["python", "-c", "print('hello world')"]
    assert captured["check"] is True
    assert captured["cwd"] == Path(__file__).resolve().parents[2]


def test_cmd_run_exits_with_subprocess_return_code(monkeypatch):
    def fake_run(command, check, cwd):
        raise subprocess.CalledProcessError(returncode=7, cmd=command)

    monkeypatch.setattr(cmd.subprocess, "run", fake_run)

    with pytest.raises(SystemExit) as exc_info:
        cmd.run(Namespace(command=["python", "-V"]))

    assert exc_info.value.code == 7
