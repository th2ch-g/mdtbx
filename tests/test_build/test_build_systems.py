import argparse

from src.build import build_solution, build_vacuum


def _parse_args(add_subcmd, argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    add_subcmd(subparsers)
    return parser.parse_args(argv)


def test_build_solution_default_template_exists_and_outdir_is_created(
    tmp_path, monkeypatch
):
    commands = []

    def fake_run(command, shell, check):
        commands.append(command)

    outdir = tmp_path / "build_solution" / "output"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(build_solution.subprocess, "run", fake_run)

    args = _parse_args(build_solution.add_subcmd, ["build_solution", "-o", str(outdir)])

    assert args.template_tleap.endswith("src/utils/template_tleap.in")

    build_solution.run(args)

    assert outdir.exists()
    assert (tmp_path / "tleap.in").exists()
    assert commands[0] == "tleap -f tleap.in"


def test_build_vacuum_creates_outdir(tmp_path, sample_pdb_path, monkeypatch):
    commands = []

    def fake_run(command, shell, check):
        commands.append(command)

    outdir = tmp_path / "build_vacuum" / "output"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(build_vacuum.subprocess, "run", fake_run)

    args = _parse_args(
        build_vacuum.add_subcmd,
        ["build_vacuum", "-i", str(sample_pdb_path), "-o", str(outdir)],
    )

    build_vacuum.run(args)

    assert outdir.exists()
    assert (tmp_path / "tleap.in").exists()
    assert commands[0] == "tleap -f tleap.in"
