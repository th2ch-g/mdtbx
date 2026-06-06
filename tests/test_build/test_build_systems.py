import argparse

from src.build import build_solution, build_vacuum


def _parse_args(add_subcmd, argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    add_subcmd(subparsers)
    return parser.parse_args(argv)


def _fake_tleap_factory():
    """Capture the tleap input passed to run_tleap (which is mocked out)."""
    calls = {}

    def fake_run_tleap(input_text, *, keepfiles=False, extra_cleanup=()):
        calls["input_text"] = input_text
        calls["keepfiles"] = keepfiles

    return calls, fake_run_tleap


def test_build_solution_default_template_exists_and_outdir_is_created(
    tmp_path, monkeypatch
):
    calls, fake = _fake_tleap_factory()

    outdir = tmp_path / "build_solution" / "output"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(build_solution, "run_tleap", fake)

    args = _parse_args(build_solution.add_subcmd, ["build_solution", "-o", str(outdir)])

    assert args.template_tleap.endswith("src/utils/template_tleap.in")

    build_solution.run(args)

    assert outdir.exists()
    assert "input_text" in calls


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
