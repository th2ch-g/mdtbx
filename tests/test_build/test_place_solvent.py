import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mdtbx.build import place_solvent


def _rism_args(**overrides):
    values = {
        "closure": "kh",
        "grdspc": 0.5,
        "tolerance": 1e-5,
        "buffer": 14.0,
        "solvcut": 14.0,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def _run_args(tmp_path, **overrides):
    prmtop = tmp_path / "solute.parm7"
    coord = tmp_path / "solute.rst7"
    xvv = tmp_path / "water.xvv"
    prmtop.write_text("topology")
    coord.write_text("coordinates")
    xvv.write_text("susceptibility")
    values = {
        "prmtop": str(prmtop),
        "coord": str(coord),
        "output": str(tmp_path / "placed.pdb"),
        "xvv": str(xvv),
        "xvv_output": None,
        "solvent": "water",
        "solvent_model": "cSPCE",
        "solvent_density": 55.5,
        "dielectric": 78.44,
        "temperature": 300.0,
        "closure": "kh",
        "grdspc": 0.5,
        "tolerance": 1e-5,
        "buffer": 14.0,
        "solvcut": 14.0,
        "threshold": 1.5,
        "exclusion_radius": 2.6,
        "max_sites": None,
        "use_sander": False,
        "keepfiles": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_run_rism1d_writes_current_namelist_and_uses_argv(tmp_path, monkeypatch):
    amberhome = tmp_path / "amber"
    model = amberhome / "dat" / "rism1d" / "mdl" / "cSPCE.mdl"
    model.parent.mkdir(parents=True)
    model.write_text("model")
    workdir = tmp_path / "work"
    workdir.mkdir()
    monkeypatch.setenv("AMBERHOME", str(amberhome))

    calls = []

    def fake_run_cmd(command, **kwargs):
        calls.append((command, kwargs))
        (Path(kwargs["cwd"]) / f"{command[1]}.xvv").write_text("xvv")

    monkeypatch.setattr(place_solvent, "run_cmd", fake_run_cmd)

    result = place_solvent._run_rism1d(
        "cSPCE",
        300.0,
        str(workdir),
        solvent_density=54.0,
        dielectric=77.0,
    )

    assert result == str(workdir / "cSPCE_300.00.xvv")
    assert calls[0][0] == ["rism1d", "cSPCE_300.00"]
    assert calls[0][1]["stderr"] is subprocess.STDOUT
    input_text = (workdir / "cSPCE_300.00.inp").read_text()
    assert "OUTLIST='x'" in input_text
    assert "TEMPERATURE=300.0, DIEPS=77.0, NSP=1" in input_text
    assert "DENSITY=54.0, UNITS='M'" in input_text
    assert f"MODEL='{model.resolve()}'" in input_text
    assert "OUTLST" not in input_text


def test_run_rism3d_generates_pdb_and_requests_dx(tmp_path, monkeypatch):
    prmtop = tmp_path / "solute.parm7"
    coord = tmp_path / "solute.rst7"
    xvv = tmp_path / "water.xvv"
    for path in (prmtop, coord, xvv):
        path.write_text("data")
    workdir = tmp_path / "work"
    workdir.mkdir()
    calls = []

    def fake_run_cmd(command, **kwargs):
        calls.append((command, kwargs))

    monkeypatch.setattr(place_solvent, "run_cmd", fake_run_cmd)

    place_solvent._run_rism3d_snglpnt(
        str(prmtop),
        str(coord),
        str(xvv),
        _rism_args(),
        str(workdir),
    )

    assert calls[0][0] == [
        "ambpdb",
        "-p",
        str(prmtop),
        "-c",
        str(coord),
    ]
    command = calls[1][0]
    assert command[command.index("--pdb") + 1] == str(workdir / "solute.pdb")
    assert command[command.index("--volfmt") + 1] == "dx"
    assert command[command.index("--xvv") + 1] == str(xvv)
    assert command[command.index("--guv") + 1] == str(workdir / "solute")
    assert "--ntwrism" not in command


def test_run_sander_supplies_xvv_and_guv(tmp_path, monkeypatch):
    prmtop = tmp_path / "solute.parm7"
    coord = tmp_path / "solute.rst7"
    xvv = tmp_path / "water.xvv"
    for path in (prmtop, coord, xvv):
        path.write_text("data")
    workdir = tmp_path / "work"
    workdir.mkdir()
    calls = []
    monkeypatch.setattr(
        place_solvent,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )

    place_solvent._run_sander_rism(
        str(prmtop),
        str(coord),
        str(xvv),
        _rism_args(),
        str(workdir),
    )

    command = calls[0][0]
    assert command[command.index("-xvv") + 1] == str(xvv)
    assert command[command.index("-guv") + 1] == str(workdir / "solute.kh")


def test_parse_dx_validates_value_count(tmp_path):
    dx_path = tmp_path / "grid.dx"
    dx_path.write_text(
        "object 1 class gridpositions counts 2 1 1\n"
        "origin 0 0 0\n"
        "delta 1 0 0\n"
        "delta 0 1 0\n"
        "delta 0 0 1\n"
        "object 2 class gridconnections counts 2 1 1\n"
        "object 3 class array type double rank 0 items 2 data follows\n"
        "1.0\n"
    )

    with pytest.raises(ValueError, match="Expected 2 grid values"):
        place_solvent._parse_dx(dx_path)


def test_find_oxygen_dx_ignores_other_sites_and_correlations(tmp_path):
    (tmp_path / "solute.H1.1.dx").write_text("hydrogen")
    (tmp_path / "solute.cuv.O.1.dx").write_text("correlation")
    oxygen = tmp_path / "solute.O.1.dx"
    oxygen.write_text("oxygen")

    result = place_solvent._find_oxygen_dx(tmp_path, "solute", "kh")

    assert result == oxygen


def test_find_oxygen_dx_does_not_fall_back_to_wrong_site(tmp_path):
    (tmp_path / "solute.H1.1.dx").write_text("hydrogen")

    with pytest.raises(FileNotFoundError, match="oxygen pair-distribution"):
        place_solvent._find_oxygen_dx(tmp_path, "solute", "kh")


def test_extract_peaks_excludes_points_on_radius_boundary():
    data = np.array([[[2.0]], [[3.0]]])

    coords, values = place_solvent._extract_peaks_greedy(
        data,
        np.zeros(3),
        np.eye(3),
        threshold=1.0,
        exclusion_radius=1.0,
    )

    np.testing.assert_allclose(coords, [[1.0, 0.0, 0.0]])
    np.testing.assert_allclose(values, [3.0])


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"exclusion_radius": 0.0}, "exclusion_radius"),
        ({"max_sites": 0}, "max_sites"),
        ({"threshold": float("nan")}, "threshold"),
    ],
)
def test_extract_peaks_rejects_invalid_options(kwargs, message):
    options = {
        "threshold": 1.0,
        "exclusion_radius": 1.0,
        "max_sites": None,
    }
    options.update(kwargs)

    with pytest.raises(ValueError, match=message):
        place_solvent._extract_peaks_greedy(
            np.ones((1, 1, 1)),
            np.zeros(3),
            np.eye(3),
            **options,
        )


def test_run_copies_xvv_creates_output_parent_and_cleans_workdir(tmp_path, monkeypatch):
    workdir = tmp_path / "work"
    workdir.mkdir()
    dx_path = tmp_path / "density.dx"
    dx_path.write_text(
        "object 1 class gridpositions counts 1 1 1\n"
        "origin 0 0 0\n"
        "delta 1 0 0\n"
        "delta 0 1 0\n"
        "delta 0 0 1\n"
        "object 3 class array type double rank 0 items 1 data follows\n"
        "2.0\n"
    )
    args = _run_args(
        tmp_path,
        output=str(tmp_path / "output" / "placed.pdb"),
        xvv_output=str(tmp_path / "cache" / "water.xvv"),
    )
    monkeypatch.setattr(
        place_solvent.tempfile, "mkdtemp", lambda **kwargs: str(workdir)
    )
    monkeypatch.setattr(place_solvent, "_run_rism3d_snglpnt", lambda *args: None)
    monkeypatch.setattr(place_solvent, "_find_oxygen_dx", lambda *args: dx_path)

    place_solvent.run(args)

    assert (tmp_path / "output" / "placed.pdb").is_file()
    assert (tmp_path / "cache" / "water.xvv").read_text() == "susceptibility"
    assert not workdir.exists()


def test_run_rejects_non_positive_numeric_option_before_creating_workdir(
    tmp_path, monkeypatch
):
    args = _run_args(tmp_path, temperature=0.0)
    monkeypatch.setattr(
        place_solvent.tempfile,
        "mkdtemp",
        lambda **kwargs: pytest.fail("workdir should not be created"),
    )

    with pytest.raises(ValueError, match="--temperature"):
        place_solvent.run(args)
