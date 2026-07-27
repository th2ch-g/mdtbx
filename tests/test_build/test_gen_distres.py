import argparse

import pytest

from mdtbx.build import gen_distres


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


def test_gen_distres_rejects_empty_selections(tmp_path, trajectory_files):
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
            ",",
            "-s2",
            ",",
        ]
    )

    with pytest.raises(ValueError, match="must not be empty"):
        gen_distres.run(args)


@pytest.mark.parametrize(
    ("option", "values"),
    [
        ("lower_bound", [-0.1]),
        ("lower_bound", [float("nan")]),
        ("upper_bound1", [0.5]),
        ("upper_bound2", [0.2]),
    ],
)
def test_gen_distres_rejects_invalid_bounds(tmp_path, trajectory_files, option, values):
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
            "index 0",
            "-s2",
            "index 1",
        ]
    )
    setattr(args, option, values)

    with pytest.raises(ValueError, match="bounds for selection 1"):
        gen_distres.run(args)


def test_gen_distres_uses_topology_relative_include_and_valid_macro(
    tmp_path, trajectory_files
):
    topology = tmp_path / "topology" / "system.top"
    topology.parent.mkdir()
    topology.write_text("; test topology")
    output_prefix = tmp_path / "restraints" / "protein-ca"
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
            "-o",
            str(output_prefix),
        ]
    )

    gen_distres.run(args)
    gen_distres.run(args)

    output = output_prefix.with_suffix(".itp")
    assert output.is_file()
    assert "#ifdef PROTEIN_CA" in output.read_text()
    assert "PROTEIN_CA_FC" in output.read_text()
    include_line = '#include "../restraints/protein-ca.itp"'
    assert topology.read_text().count(include_line) == 1
