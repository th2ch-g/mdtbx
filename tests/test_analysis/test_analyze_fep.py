import json
from types import SimpleNamespace

import pytest

from mdtbx.analysis import analyze_fep


def _fep_setup(tmp_path):
    base = tmp_path / "fep"
    base.mkdir()
    windows = []
    for index in range(2):
        name = f"lambda_{index:03d}"
        directory = base / name
        directory.mkdir()
        (directory / "fep.xvg").write_text("Delta H data\n")
        windows.append({"index": index, "directory": name})
    manifest = {
        "schema_version": 1,
        "mode": "transform",
        "deffnm": "fep",
        "windows": windows,
    }
    (base / "fep_manifest.json").write_text(json.dumps(manifest))
    return base


def _args(base, **overrides):
    values = {
        "path": str(base),
        "begin": 0.0,
        "end": None,
        "temperature": None,
        "output_prefix": None,
        "gmx": "gmx",
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_parse_total_accepts_gromacs_bar_table():
    stdout = "total  0 - 2,  DG 1.250 +/- 0.200\n"

    assert analyze_fep._parse_total(stdout) == (1.25, 0.2)


def test_run_executes_bar_and_writes_machine_readable_result(tmp_path, monkeypatch):
    base = _fep_setup(tmp_path)
    captured = {}

    def fake_run_cmd(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return SimpleNamespace(
            stdout="BAR output\n",
            stderr="total  0 - 2,  DG -4.184 +/- 0.4184\n",
        )

    monkeypatch.setattr(analyze_fep, "run_cmd", fake_run_cmd)

    analyze_fep.run(
        _args(
            base,
            begin=100.0,
            end=500.0,
            temperature=310.0,
            output_prefix=str(tmp_path / "result.v1"),
        )
    )

    command = captured["command"]
    assert command[:2] == ["gmx", "bar"]
    assert command[command.index("-b") + 1] == "100.0"
    assert command[command.index("-e") + 1] == "500.0"
    assert command[command.index("-temp") + 1] == "310.0"
    assert captured["kwargs"] == {"capture_output": True, "text": True}

    output = json.loads((tmp_path / "result.v1.json").read_text())
    assert output["delta_g_kj_mol"] == -4.184
    assert output["delta_g_kcal_mol"] == pytest.approx(-1.0)
    assert output["windows"] == 2
    assert (tmp_path / "result.v1.log").is_file()


def test_run_reports_missing_delta_h_file(tmp_path):
    base = _fep_setup(tmp_path)
    (base / "lambda_001" / "fep.xvg").unlink()

    with pytest.raises(FileNotFoundError, match="Missing Delta H"):
        analyze_fep.run(_args(base))


def test_run_rejects_invalid_time_range(tmp_path):
    base = _fep_setup(tmp_path)

    with pytest.raises(ValueError, match="greater than"):
        analyze_fep.run(_args(base, begin=10.0, end=5.0))
