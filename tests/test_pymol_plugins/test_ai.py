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
