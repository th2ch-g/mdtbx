import argparse

import pytest

from mdtbx.build import build_solution, build_vacuum


def _parse_args(add_subcmd, argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    add_subcmd(subparsers)
    return parser.parse_args(argv)


def _fake_tleap_factory():
    """Capture the tleap input passed to run_tleap (which is mocked out)."""
    calls = {}

    def fake_run_tleap(input_text, *, keepfiles=False, extra_cleanup=(), cwd=None):
        calls["input_text"] = input_text
        calls["keepfiles"] = keepfiles
        calls["cwd"] = cwd

    return calls, fake_run_tleap


def test_build_solution_default_template_exists_and_outdir_is_created(
    tmp_path, monkeypatch
):
    calls, fake = _fake_tleap_factory()

    outdir = tmp_path / "build_solution" / "output"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(build_solution, "run_tleap", fake)
    monkeypatch.setattr(build_solution, "repack_water", lambda *_args, **_kwargs: None)

    args = _parse_args(build_solution.add_subcmd, ["build_solution", "-o", str(outdir)])

    assert args.template_tleap.endswith("src/mdtbx/utils/template_tleap.in")

    build_solution.run(args)

    assert outdir.exists()
    assert "input_text" in calls
    assert calls["cwd"] == outdir
    assert "source leaprc.gaff2" in calls["input_text"]
    exact_box_command = f"set {build_solution.SYSTEM_NAME} box {{ 100 100 100 }}"
    assert calls["input_text"].count(exact_box_command) == 2
    assert f"setbox {build_solution.SYSTEM_NAME} vdw" not in calls["input_text"]


def test_build_solution_repacks_water_with_packmol(tmp_path, monkeypatch):
    calls, fake = _fake_tleap_factory()
    repack_calls = []
    outdir = tmp_path / "output"
    monkeypatch.setattr(build_solution, "run_tleap", fake)
    monkeypatch.setattr(
        build_solution,
        "repack_water",
        lambda *args, **kwargs: (
            repack_calls.append((args, kwargs))
            or {"vectorized_transfer": True, "saved_max_abs_error_A": 0.0}
        ),
    )
    args = _parse_args(
        build_solution.add_subcmd,
        [
            "build_solution",
            "-o",
            str(outdir),
            "--water-seed",
            "42",
            "--packmol-tolerance",
            "1.8",
        ],
    )

    build_solution.run(args)

    assert "input_text" in calls
    assert calls["cwd"] == outdir
    assert "source leaprc.gaff2" in calls["input_text"]
    assert repack_calls == [
        (
            (
                outdir / "leap.parm7",
                outdir / "leap.rst7",
                outdir / "leap.pdb",
            ),
            {"seed": 42, "tolerance": 1.8},
        )
    ]
    assert (outdir / "packmol_transfer.json").exists()


def test_build_solution_rejects_non_positive_packmol_tolerance(tmp_path):
    args = _parse_args(
        build_solution.add_subcmd,
        [
            "build_solution",
            "-o",
            str(tmp_path / "output"),
            "--packmol-tolerance",
            "0",
        ],
    )

    with pytest.raises(ValueError, match="--packmol-tolerance must be positive"):
        build_solution.run(args)


@pytest.mark.parametrize(
    ("extra_args", "message"),
    [
        (["--ion_conc", "-0.1"], "--ion_conc must be non-negative"),
        (["--water-seed", "-1"], "--water-seed must be non-negative"),
        (
            ["--boxsize", "2", "10", "10", "--packmol-tolerance", "2"],
            "Each --boxsize edge",
        ),
    ],
)
def test_build_solution_rejects_invalid_numeric_inputs(tmp_path, extra_args, message):
    args = _parse_args(
        build_solution.add_subcmd,
        ["build_solution", "-o", str(tmp_path / "output"), *extra_args],
    )

    with pytest.raises(ValueError, match=message):
        build_solution.run(args)


def test_build_vacuum_creates_outdir(tmp_path, sample_pdb_path, monkeypatch):
    calls, fake = _fake_tleap_factory()

    outdir = tmp_path / "build_vacuum" / "output"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(build_vacuum, "run_tleap", fake)

    args = _parse_args(
        build_vacuum.add_subcmd,
        ["build_vacuum", "-i", str(sample_pdb_path), "-o", str(outdir)],
    )

    build_vacuum.run(args)

    assert outdir.exists()
    assert "input_text" in calls
    assert calls["cwd"] == outdir
    assert "source leaprc.gaff2" in calls["input_text"]


def test_center_pdb_in_box_moves_bounding_box_to_center(tmp_path):
    source = tmp_path / "source.pdb"
    output = tmp_path / "centered.pdb"
    source.write_text(
        "ATOM      1  C1  LIG A   1      10.000  20.000  30.000  1.00  0.00           C\n"
        "ATOM      2  C2  LIG A   1      14.000  26.000  38.000  1.00  0.00           C\n"
    )

    build_solution._center_pdb_in_box(source, output, [20.0, 30.0, 40.0])

    lines = output.read_text().splitlines()
    coordinates = [
        tuple(float(line[start : start + 8]) for start in (30, 38, 46))
        for line in lines
    ]
    assert coordinates == [(8.0, 12.0, 16.0), (12.0, 18.0, 24.0)]


def test_center_pdb_in_box_rejects_oversized_solute(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(
        "ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  C2  LIG A   1      20.000   1.000   1.000  1.00  0.00           C\n"
    )

    with pytest.raises(ValueError, match="smaller than --boxsize"):
        build_solution._center_pdb_in_box(
            source, tmp_path / "centered.pdb", [20.0, 30.0, 40.0]
        )
