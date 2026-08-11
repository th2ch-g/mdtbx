import os
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mdtbx.analysis import analyze_fep_rest


def _energy(path, values):
    path.write_text(
        '@ title "Potential"\n'
        + "".join(f"{index}.0 {value}\n" for index, value in enumerate(values))
    )
    return path


def test_calculate_bar_uses_cross_hamiltonian_energies(tmp_path):
    paths = {
        (0, 0): _energy(tmp_path / "00.xvg", [0.0, 1.0, 2.0]),
        (0, 1): _energy(tmp_path / "01.xvg", [2.0, 3.0, 4.0]),
        (1, 1): _energy(tmp_path / "11.xvg", [5.0, 6.0, 7.0]),
        (1, 0): _energy(tmp_path / "10.xvg", [3.0, 4.0, 5.0]),
    }

    pairs, total, uncertainty = analyze_fep_rest._calculate_bar(
        paths,
        2,
        300.0,
        0.0,
        None,
        1,
    )

    assert total == pytest.approx(2.0)
    assert uncertainty == pytest.approx(0.0, abs=1e-7)
    assert pairs[0]["forward_samples"] == 3


def test_rerun_energy_uses_cache(tmp_path):
    output_dir = tmp_path / "rerun"
    output_dir.mkdir()
    trajectory = tmp_path / "trajectory.trr"
    tpr = tmp_path / "state.tpr"
    trajectory.write_text("trajectory")
    tpr.write_text("tpr")
    cached = output_dir / "potential.xvg"
    cached.write_text("0 1\n")

    result = analyze_fep_rest._rerun_energy(
        type("Args", (), {"force": False, "gmx": "gmx"})(),
        trajectory,
        tpr,
        output_dir,
    )

    assert result == cached


def test_rerun_energy_invalidates_cache_after_trajectory_update(tmp_path, monkeypatch):
    output_dir = tmp_path / "rerun"
    output_dir.mkdir()
    trajectory = tmp_path / "trajectory.trr"
    tpr = tmp_path / "state.tpr"
    cached = output_dir / "potential.xvg"
    trajectory.write_text("trajectory")
    tpr.write_text("tpr")
    cached.write_text("0 1\n")
    os.utime(cached, ns=(1_000_000_000, 1_000_000_000))
    os.utime(trajectory, ns=(2_000_000_000, 2_000_000_000))
    calls = []

    def fake_run(command, cwd=None, **_kwargs):
        calls.append(command)
        if command[1] == "energy":
            Path(cwd, "potential.xvg").write_text("0 2\n")

    monkeypatch.setattr(analyze_fep_rest, "run_cmd", fake_run)
    args = type(
        "Args",
        (),
        {"force": False, "gmx": "gmx_mpi", "gpu_offload": False, "ntomp": 1},
    )()

    analyze_fep_rest._rerun_energy(args, trajectory, tpr, output_dir)

    assert [command[1] for command in calls] == ["mdrun", "energy"]


def test_rerun_energy_sets_openmp_threads(tmp_path, monkeypatch):
    output_dir = tmp_path / "rerun"
    calls = []

    def fake_run(command, cwd=None, **_kwargs):
        calls.append((command, cwd))
        if command[1] == "mdrun":
            Path(cwd, "rerun.trr").write_text("temporary trajectory")
        elif command[1] == "energy":
            Path(cwd, "potential.xvg").write_text("0 1\n")

    monkeypatch.setattr(analyze_fep_rest, "run_cmd", fake_run)
    args = type(
        "Args",
        (),
        {"force": False, "gmx": "gmx_mpi", "gpu_offload": False, "ntomp": 8},
    )()

    analyze_fep_rest._rerun_energy(
        args,
        Path("trajectory.trr"),
        Path("state.tpr"),
        output_dir,
    )

    assert calls[0][0][calls[0][0].index("-ntomp") + 1] == "8"
    assert not (output_dir / "rerun.trr").exists()


def test_run_reuses_rerun_energies_for_convergence(tmp_path, monkeypatch):
    windows = []
    for index in range(2):
        directory = tmp_path / f"lambda_{index:03d}"
        directory.mkdir()
        (directory / "rest.trr").write_text("trajectory")
        (directory / "rest.tpr").write_text("tpr")
        windows.append({"index": index, "directory": directory.name})
    (tmp_path / "fep_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "workflow": "fep-rest",
                "mode": "transform",
                "deffnm": "rest",
                "physical_temperature": 300.0,
                "windows": windows,
            }
        )
    )
    energy = _energy(tmp_path / "energy.xvg", [0, 0, 0, 0, 0])
    rerun_calls = []
    bar_calls = []

    def fake_rerun(*rerun_args):
        rerun_calls.append(rerun_args)
        return energy

    def fake_bar(_paths, _count, _temperature, begin, end, _subsample):
        bar_calls.append((begin, end))
        return [], float(begin), 0.1

    monkeypatch.setattr(analyze_fep_rest, "_rerun_energy", fake_rerun)
    monkeypatch.setattr(analyze_fep_rest, "_calculate_bar", fake_bar)
    analyze_fep_rest.run(
        SimpleNamespace(
            path=str(tmp_path),
            begin=0.0,
            end=None,
            subsample=1,
            output=None,
            gmx="gmx_mpi",
            force=False,
            gpu_offload=False,
            ntomp=1,
            convergence_blocks=2,
        )
    )

    output = json.loads((tmp_path / "fep_rest_result.json").read_text())
    assert len(rerun_calls) == 4
    assert len(bar_calls) == 4
    assert len(output["convergence"]["block_estimates"]) == 2
