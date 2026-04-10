import argparse

import pytest

from src.build import gen_distres


def _parse_args(argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    gen_distres.add_subcmd(subparsers)
    return parser.parse_args(argv)


def test_gen_distres_uses_default_bounds(tmp_path, trajectory_files, monkeypatch):
    topology = tmp_path / "system.top"
    topology.write_text("; test topology\n")

    monkeypatch.chdir(tmp_path)

    args = _parse_args(
        [
            "gen_distres",
            "-g",
            trajectory_files["pdb"],
            "-p",
            str(topology),
            "-s1",
            "index 0",
            "-s2",
            "index 1",
        ]
    )

    gen_distres.run(args)

    output = tmp_path / "distres.itp"
    assert output.exists()
    content = output.read_text()
    assert "1 2 1 0 1 0.0 0.3 0.4 DISTRES_FC" in content
    assert '#include "distres.itp"' in topology.read_text()


def test_gen_distres_rejects_mismatched_selection_counts(tmp_path, trajectory_files):
    topology = tmp_path / "system.top"
    topology.write_text("; test topology\n")

    args = _parse_args(
        [
            "gen_distres",
            "-g",
            trajectory_files["pdb"],
            "-p",
            str(topology),
            "-s1",
            "index 0,index 1",
            "-s2",
            "index 2",
        ]
    )

    with pytest.raises(ValueError, match="selection1 and selection2"):
        gen_distres.run(args)
