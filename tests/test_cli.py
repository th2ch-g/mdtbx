"""
CLI subcommand registration tests

Confirms that every subcommand registers with argparse. No external tool
(PyMOL, Gromacs, ...) is invoked.
"""

import subprocess
import sys

import pytest


def test_cli_importable():
    """The cli module imports without error"""
    from mdtbx.cli import cli, create_parser

    assert callable(cli)
    assert callable(create_parser)


def test_artifact_paths_include_mutated_path_directory():
    """Direct --json runs track --path for mutating commands like agent_run does"""
    import types

    from mdtbx.cli import _artifact_paths

    args = types.SimpleNamespace(_command="run_fep", path="fep", output=None)
    assert "fep" in _artifact_paths(args)

    args = types.SimpleNamespace(_command="fit", path="fep", output=None)
    assert "fep" not in _artifact_paths(args)


def test_mdtraj_pkg_resources_warning_is_hidden():
    result = subprocess.run(
        [sys.executable, "-m", "mdtbx", "fit", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "pkg_resources is deprecated as an API" not in result.stderr


def test_all_subcommands_registered():
    """
    Every subcommand is registered with argparse and its --help works
    (exiting with SystemExit(0))
    """
    from mdtbx.cli import cli

    with pytest.raises(SystemExit) as exc_info:
        # --help raises SystemExit(0)
        import sys

        sys.argv = ["mdtbx", "--help"]
        cli()
    assert exc_info.value.code == 0


@pytest.mark.parametrize(
    "subcmd",
    [
        "addace",
        "addh",
        "mutate",
        "fill_chainname",
        "addnme",
        "add_ndx",
        "mv_crds_mol2",
        "gen_am1bcc",
        "gen_resp",
        "gen_modres_am1bcc",
        "gen_modres_resp",
        "gen_posres",
        "gen_distres",
        "modeling_cf",
        "find_bond",
        "mod_mdp",
        "convert",
        "calc_ion_conc",
        "centering_gro",
        "amb2gro",
        "trjcat",
        "fit",
        "pacs_trjcat",
        "rmfile",
        "extract_ave_str",
        "extract_str",
        "show_mdtraj",
        "show_npy",
        "print_perf",
        "opt_perf",
        "shell_hook",
        "partial_tempering",
        "gen_temperatures",
        "cmd",
        "build_solution",
        "build_vacuum",
        "place_solvent",
        "setup_fep",
        "setup_fep_rest",
        "setup_abfe",
        "run_fep",
        "run_abfe",
        "analyze_fep",
        "optimize_fep_schedule",
        "analyze_fep_rest",
        "analyze_abfe",
        "contactmap",
        "comdist",
        "comvec",
        "mindist",
        "rmsd",
        "rmsf",
        "xyz",
        "pca",
        "densmap",
        "distmap",
        "tica",
        "cluster",
        "msm",
    ],
)
def test_subcommand_help(subcmd):
    """--help of each subcommand exits successfully"""
    import sys

    from mdtbx.cli import cli

    sys.argv = ["mdtbx", subcmd, "--help"]
    with pytest.raises(SystemExit) as exc_info:
        cli()
    assert exc_info.value.code == 0
