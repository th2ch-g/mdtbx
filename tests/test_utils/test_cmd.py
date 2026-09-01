import subprocess
from argparse import Namespace
from pathlib import Path

import pytest

from mdtbx.utils import cmd


def test_cmd_run_uses_argv_list(monkeypatch):
    captured = {}

    def fake_run(command, check, **kwargs):
        captured["command"] = command
        captured["check"] = check

    monkeypatch.setattr(cmd.subprocess, "run", fake_run)

    cmd.run(Namespace(command=["python", "-c", "print('hello world')"]))

    assert isinstance(captured["command"], list)
    assert captured["command"][:3] == ["pixi", "run", "--manifest-path"]
    assert captured["command"][4:] == ["python", "-c", "print('hello world')"]
    assert captured["check"] is True


def test_cmd_run_exits_with_subprocess_return_code(monkeypatch):
    def fake_run(command, check, **kwargs):
        raise subprocess.CalledProcessError(returncode=7, cmd=command)

    monkeypatch.setattr(cmd.subprocess, "run", fake_run)

    with pytest.raises(SystemExit) as exc_info:
        cmd.run(Namespace(command=["python", "-V"]))

    assert exc_info.value.code == 7


def test_cmd_packmol_memgen_adds_python_compatibility_path(monkeypatch):
    captured = {}

    def fake_run(command, check, **kwargs):
        captured["env"] = kwargs["env"]

    monkeypatch.setattr(cmd.subprocess, "run", fake_run)

    cmd.run(Namespace(command=["packmol-memgen", "--help"]))

    compatibility_path = captured["env"]["PYTHONPATH"].split(":", maxsplit=1)[0]
    assert Path(compatibility_path).name == "packmol_memgen"
    assert (Path(compatibility_path) / "sitecustomize.py").is_file()


def test_cmd_packmol_memgen_prefers_amber_packmol(monkeypatch, tmp_path):
    captured = {}
    packmol = tmp_path / "bin" / "packmol"
    packmol.parent.mkdir()
    packmol.write_text("")
    packmol.chmod(0o755)

    def fake_run(command, check, **kwargs):
        captured["command"] = command

    monkeypatch.setenv("AMBERHOME", str(tmp_path))
    monkeypatch.setattr(cmd.subprocess, "run", fake_run)

    cmd.run(Namespace(command=["packmol-memgen", "--help"]))

    assert captured["command"][-2:] == ["--packmol", str(packmol)]
