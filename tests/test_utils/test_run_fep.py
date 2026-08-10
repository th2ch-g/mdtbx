import json
from types import SimpleNamespace

import pytest

from mdtbx.utils import run_fep


def _fep_setup(tmp_path, count=3):
    base = tmp_path / "fep"
    base.mkdir()
    windows = []
    for index in range(count):
        name = f"lambda_{index:03d}"
        directory = base / name
        directory.mkdir()
        (directory / "fep.tpr").write_text("tpr")
        windows.append({"index": index, "directory": name})
    manifest = {
        "schema_version": 1,
        "mode": "decouple",
        "deffnm": "fep",
        "windows": windows,
    }
    (base / "fep_manifest.json").write_text(json.dumps(manifest))
    return base


def _args(base, **overrides):
    values = {
        "path": str(base),
        "start_window": 0,
        "end_window": None,
        "multidir": False,
        "replex": 0,
        "ntmpi": None,
        "ntomp": None,
        "maxh": None,
        "nsteps": None,
        "no_continue": False,
        "gmx": "gmx",
        "mpi_launcher": None,
        "gpu_offload": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_run_selected_windows_sequentially_and_continues(tmp_path, monkeypatch):
    base = _fep_setup(tmp_path)
    (base / "lambda_001" / "fep.cpt").write_text("checkpoint")
    calls = []
    monkeypatch.setattr(
        run_fep,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )

    run_fep.run(
        _args(
            base,
            start_window=1,
            end_window=2,
            ntomp=2,
            nsteps=10,
        )
    )

    assert len(calls) == 2
    assert calls[0][1]["cwd"] == (base / "lambda_001").resolve()
    assert calls[0][0][-2:] == ["-cpi", "fep.cpt"]
    assert "-cpi" not in calls[1][0]
    assert calls[0][0][calls[0][0].index("-ntomp") + 1] == "2"
    assert calls[0][0][calls[0][0].index("-nsteps") + 1] == "10"


def test_multidir_uses_one_mdrun_invocation(tmp_path, monkeypatch):
    base = _fep_setup(tmp_path)
    calls = []
    monkeypatch.setattr(
        run_fep,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )

    run_fep.run(_args(base, multidir=True, replex=100))

    assert len(calls) == 1
    command = calls[0][0]
    assert command[:3] == ["gmx", "mdrun", "-multidir"]
    assert "-replex" in command
    assert command[command.index("-replex") + 1] == "100"


def test_gpu_offload_adds_validated_gromacs_options(tmp_path, monkeypatch):
    base = _fep_setup(tmp_path)
    calls = []
    monkeypatch.setattr(
        run_fep,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )

    run_fep.run(_args(base, gpu_offload=True))

    command = calls[0][0]
    assert command[command.index("-nb") + 1] == "gpu"
    assert command[command.index("-pme") + 1] == "gpu"
    assert command[command.index("-update") + 1] == "cpu"


def test_multidir_rejects_partial_checkpoint_set(tmp_path):
    base = _fep_setup(tmp_path)
    (base / "lambda_000" / "fep.cpt").write_text("checkpoint")

    with pytest.raises(ValueError, match="All selected windows"):
        run_fep.run(_args(base, multidir=True))


def test_run_rejects_invalid_window_range(tmp_path):
    base = _fep_setup(tmp_path)

    with pytest.raises(ValueError, match="Window range"):
        run_fep.run(_args(base, start_window=2, end_window=1))


def test_fep_rest_adds_plumed_hrex_options(tmp_path, monkeypatch):
    base = _fep_setup(tmp_path)
    manifest_path = base / "fep_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    plumed = base / "plumed.dat"
    plumed.write_text("# HREX\n")
    manifest.update(
        {
            "workflow": "fep-rest",
            "plumed_file": str(plumed),
            "gmx_executable": "gmx_mpi",
        }
    )
    manifest_path.write_text(json.dumps(manifest))
    calls = []
    monkeypatch.setattr(
        run_fep,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )

    run_fep.run(_args(base, multidir=True, replex=100, gmx=None))

    command = calls[0][0]
    assert command[:5] == ["mpirun", "-np", "3", "gmx_mpi", "mdrun"]
    assert "-plumed" in command
    assert command[command.index("-plumed") + 1] == str(plumed.resolve())
    assert "-hrex" in command
    assert command[command.index("-dlb") + 1] == "no"


def test_fep_rest_accepts_scheduler_launcher_template(tmp_path, monkeypatch):
    base = _fep_setup(tmp_path)
    manifest_path = base / "fep_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    plumed = base / "plumed.dat"
    plumed.write_text("# HREX\n")
    manifest.update({"workflow": "fep-rest", "plumed_file": str(plumed)})
    manifest_path.write_text(json.dumps(manifest))
    calls = []
    monkeypatch.setattr(
        run_fep,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )

    run_fep.run(
        _args(
            base,
            multidir=True,
            replex=100,
            ntmpi=6,
            mpi_launcher="srun --ntasks {ntmpi}",
        )
    )

    command = calls[0][0]
    assert command[:4] == ["srun", "--ntasks", "6", "gmx"]
    assert "-ntmpi" not in command


def test_fep_rest_rejects_incompatible_external_rank_count(tmp_path):
    base = _fep_setup(tmp_path)
    manifest_path = base / "fep_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["workflow"] = "fep-rest"
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="multiple"):
        run_fep.run(_args(base, multidir=True, replex=100, ntmpi=4))


def test_fep_rest_requires_exchange(tmp_path):
    base = _fep_setup(tmp_path)
    manifest_path = base / "fep_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["workflow"] = "fep-rest"
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="--multidir"):
        run_fep.run(_args(base))
