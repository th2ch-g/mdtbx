import sys

from src.agent.context import JSON_MODE
from src.utils import proc


def test_json_mode_redirects_uncaptured_subprocess(monkeypatch):
    captured = {}

    def fake_run(_command, **kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(proc.subprocess, "run", fake_run)
    token = JSON_MODE.set(True)
    try:
        proc.run_cmd(["tool"])
    finally:
        JSON_MODE.reset(token)

    assert captured["stdout"] is sys.stderr
    assert captured["stderr"] is sys.stderr


def test_json_mode_preserves_capture_output(monkeypatch):
    captured = {}

    def fake_run(_command, **kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(proc.subprocess, "run", fake_run)
    token = JSON_MODE.set(True)
    try:
        proc.run_cmd(["tool"], capture_output=True)
    finally:
        JSON_MODE.reset(token)

    assert captured["capture_output"] is True
    assert "stdout" not in captured
    assert "stderr" not in captured
