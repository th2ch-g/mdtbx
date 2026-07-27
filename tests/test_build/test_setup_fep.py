import json
from types import SimpleNamespace

import pytest

from mdtbx.build import setup_fep


def _args(tmp_path, **overrides):
    mdp = tmp_path / "base.mdp"
    topology = tmp_path / "system.top"
    structure = tmp_path / "system.gro"
    mdp.write_text("integrator = md\nnsteps = 10\nfree-energy = no\n")
    topology.write_text("topology\n")
    structure.write_text("structure\n")
    values = {
        "mdp": str(mdp),
        "topology": str(topology),
        "structure": str(structure),
        "outdir": str(tmp_path / "fep"),
        "mode": "decouple",
        "moltype": "LIG",
        "windows": 12,
        "coul_windows": 2,
        "vdw_windows": 2,
        "fep_lambdas": None,
        "coul_lambdas": None,
        "vdw_lambdas": None,
        "calc_lambda_neighbors": 1,
        "nstdhdl": 10,
        "checkpoint": None,
        "reference": None,
        "b_reference": None,
        "index": None,
        "deffnm": "fep",
        "gmx": "gmx",
        "maxwarn": 0,
        "no_grompp": True,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_run_generates_staged_windows_and_manifest(tmp_path):
    args = _args(tmp_path)

    setup_fep.run(args)

    manifest = json.loads((tmp_path / "fep" / "fep_manifest.json").read_text())
    assert manifest["mode"] == "decouple"
    assert manifest["molecule_type"] == "LIG"
    assert manifest["prepared"] is False
    assert manifest["gmx_executable"] == "gmx"
    assert len(manifest["windows"]) == 3
    assert manifest["lambda_components"] == {
        "coul": [0.0, 1.0, 1.0],
        "vdw": [0.0, 0.0, 1.0],
    }

    first_mdp = (tmp_path / "fep" / "lambda_000" / "fep.mdp").read_text()
    last_mdp = (tmp_path / "fep" / "lambda_002" / "fep.mdp").read_text()
    assert "couple-moltype" in first_mdp
    assert "couple-lambda0" in first_mdp
    assert "init-lambda-state" in last_mdp
    assert "= 2" in last_mdp
    assert "free-energy" in first_mdp
    assert "free-energy = no" not in first_mdp


def test_transform_mode_uses_dual_state_topology_without_moltype(tmp_path):
    args = _args(
        tmp_path,
        mode="transform",
        moltype=None,
        windows=3,
    )

    setup_fep.run(args)

    mdp = (tmp_path / "fep" / "lambda_001" / "fep.mdp").read_text()
    assert "fep-lambdas" in mdp
    assert "couple-moltype" not in mdp


def test_run_invokes_grompp_from_topology_directory(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        setup_fep,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )
    args = _args(tmp_path, no_grompp=False)

    setup_fep.run(args)

    assert len(calls) == 3
    command, kwargs = calls[0]
    assert command[:2] == ["gmx", "grompp"]
    assert command[command.index("-o") + 1].endswith("lambda_000/fep.tpr")
    assert kwargs["cwd"] == tmp_path
    manifest = json.loads((tmp_path / "fep" / "fep_manifest.json").read_text())
    assert manifest["prepared"] is True


def test_run_rejects_nonempty_output_directory(tmp_path):
    outdir = tmp_path / "fep"
    outdir.mkdir()
    (outdir / "existing.txt").write_text("data")

    with pytest.raises(FileExistsError, match="not empty"):
        setup_fep.run(_args(tmp_path))


def test_run_requires_moltype_for_decoupling(tmp_path):
    with pytest.raises(ValueError, match="--moltype is required"):
        setup_fep.run(_args(tmp_path, moltype=None))


def test_run_rejects_zero_lambda_neighbors(tmp_path):
    with pytest.raises(ValueError, match="must be -1 or positive"):
        setup_fep.run(_args(tmp_path, calc_lambda_neighbors=0))
