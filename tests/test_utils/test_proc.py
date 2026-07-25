from src.utils import proc


def test_run_cmd_enables_text_mode_for_string_input(monkeypatch):
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return object()

    monkeypatch.setattr(proc.subprocess, "run", fake_run)

    proc.run_cmd(["gmx", "trjconv"], input="Protein\n")

    assert captured["command"] == ["gmx", "trjconv"]
    assert captured["kwargs"]["input"] == "Protein\n"
    assert captured["kwargs"]["text"] is True
    assert captured["kwargs"]["shell"] is False


def test_run_cmd_leaves_bytes_input_in_binary_mode(monkeypatch):
    captured = {}

    def fake_run(_command, **kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(proc.subprocess, "run", fake_run)

    proc.run_cmd(["tool"], input=b"data")

    assert "text" not in captured


def test_run_cmd_respects_explicit_text_mode(monkeypatch):
    captured = {}

    def fake_run(_command, **kwargs):
        captured.update(kwargs)
        return object()

    monkeypatch.setattr(proc.subprocess, "run", fake_run)

    proc.run_cmd(["tool"], input="data", text=False)

    assert captured["text"] is False
