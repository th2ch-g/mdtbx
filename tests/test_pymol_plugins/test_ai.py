import importlib.util
from pathlib import Path

import pytest


PLUGIN_PATH = (
    Path(__file__).resolve().parents[2] / "pymol-plugins" / "pymol_plugins" / "ai.py"
)


def _load_ai_module():
    spec = importlib.util.spec_from_file_location("ai_under_test", PLUGIN_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    "code",
    [
        "import os",
        "open('output.txt', 'w')",
        "cmd.do('system whoami')",
        "cmd.__class__",
        "(lambda: 1)()",
    ],
)
def test_validate_python_code_rejects_escape_paths(code):
    ai = _load_ai_module()

    with pytest.raises(RuntimeError):
        ai._validate_python_code(code)


def test_execute_python_block_allows_public_cmd_calls():
    ai = _load_ai_module()
    calls = []
    ai.cmd.color = lambda color, selection: calls.append((color, selection))

    ai._execute_blocks([("python", "cmd.color('red', 'all')")])

    assert calls == [("red", "all")]


@pytest.mark.parametrize(
    "command",
    [
        "run malicious.py",
        "@malicious.pml",
        "/print('unsafe')",
        "python",
        "system whoami",
        "color red, all; system whoami",
    ],
)
def test_validate_pymol_command_rejects_unsafe_commands(command):
    ai = _load_ai_module()

    with pytest.raises(RuntimeError):
        ai._validate_pymol_command(command)


def test_validate_pymol_command_allows_registered_command():
    ai = _load_ai_module()
    ai.cmd.keyword = {"color": object()}

    ai._validate_pymol_command("color red, all")


def test_prompt_includes_shared_cross_backend_history(tmp_path):
    ai = _load_ai_module()
    ai._record_ai_history(
        {
            "id": 1,
            "type": "claude",
            "instruction": "color the ligand red",
            "response": "color red, ligand",
            "status": "done",
            "attempts": 1,
            "executable_blocks": [{"language": "pymol", "code": "color red, ligand"}],
        }
    )
    ai._record_ai_history(
        {
            "id": 2,
            "type": "codex",
            "instruction": "show it as sticks",
            "response": "show sticks, ligand",
            "status": "done",
            "attempts": 1,
        }
    )

    prompt = ai._build_prompt(
        "make it transparent",
        "Loaded objects: ['protein']",
        tmp_path / "scene.png",
        ai._snapshot_ai_history(),
    )

    assert "color the ligand red" in prompt
    assert "show it as sticks" in prompt
    assert "claude" in prompt
    assert "codex" in prompt
    assert "make it transparent" in prompt


def test_ai_history_is_session_bounded_and_clearable(capsys):
    ai = _load_ai_module()
    for index in range(ai.AI_HISTORY_MAX_ENTRIES + 2):
        ai._record_ai_history(
            {
                "id": index,
                "type": "claude",
                "instruction": f"request {index}",
                "response": "color red, all",
                "status": "done",
                "attempts": 1,
            }
        )
    history = ai._snapshot_ai_history()
    assert len(history) == ai.AI_HISTORY_MAX_ENTRIES
    assert history[0]["instruction"] == "request 2"

    ai.ai_clear()
    assert ai._snapshot_ai_history() == []
    assert "Cleared" in capsys.readouterr().out


def test_ai_requests_are_executed_by_one_shared_worker(monkeypatch):
    import threading
    import time

    ai = _load_ai_module()
    guard = threading.Lock()
    active = 0
    maximum_active = 0
    order = []

    def fake_run(job_id, ai_type, instruction):
        nonlocal active, maximum_active
        with guard:
            active += 1
            maximum_active = max(maximum_active, active)
            order.append((job_id, ai_type, instruction))
        time.sleep(0.02)
        ai._update_ai_job(
            job_id,
            status="done",
            response="color red, all",
            attempts=1,
            executable_blocks=[{"language": "pymol", "code": "color red, all"}],
        )
        with guard:
            active -= 1

    monkeypatch.setattr(ai, "_run_ai_job", fake_run)
    ai._submit_ai_request("first", type="claude", async_="1")
    ai._submit_ai_request("second", type="codex", async_="1")
    with ai.AI_JOBS_LOCK:
        events = [job["_event"] for job in ai.AI_JOBS.values()]
    assert all(event.wait(2) for event in events)
    assert maximum_active == 1
    assert [item[2] for item in order] == ["first", "second"]
    assert [item["backend"] for item in ai._snapshot_ai_history()] == [
        "claude",
        "codex",
    ]


def test_claude_cli_has_noninteractive_tool_free_flags(monkeypatch, tmp_path):
    import subprocess

    ai = _load_ai_module()
    image = tmp_path / "scene.png"
    image.write_bytes(b"png")
    captured = {}

    def fake_run(args, stdin=None, cwd=None):
        captured.update(args=args, stdin=stdin, cwd=cwd)
        return subprocess.CompletedProcess(
            args,
            0,
            stdout='{"type":"result","result":"color red, all"}\n',
            stderr="",
        )

    monkeypatch.setattr(ai, "_run_subprocess", fake_run)
    assert ai._run_claude("prompt", image) == "color red, all"
    assert "--bare" in captured["args"]
    assert captured["args"][captured["args"].index("--tools") + 1] == ""
    assert (
        captured["args"][captured["args"].index("--permission-mode") + 1] == "dontAsk"
    )


def test_codex_cli_is_ephemeral_read_only_and_noninteractive(monkeypatch, tmp_path):
    import subprocess

    ai = _load_ai_module()
    image = tmp_path / "scene.png"
    image.write_bytes(b"png")
    captured = {}

    def fake_run(args, stdin=None, cwd=None):
        captured.update(args=args, stdin=stdin, cwd=cwd)
        output = Path(args[args.index("--output-last-message") + 1])
        output.write_text("color red, all")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    monkeypatch.setattr(ai, "_run_subprocess", fake_run)
    assert ai._run_codex("prompt", image) == "color red, all"
    assert "--ephemeral" in captured["args"]
    assert captured["args"][captured["args"].index("--sandbox") + 1] == "read-only"
    assert captured["args"][captured["args"].index("--ask-for-approval") + 1] == "never"
    assert captured["cwd"] != str(tmp_path)
