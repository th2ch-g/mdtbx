"""
CLI subcommand registration tests

Confirms that every subcommand registers with argparse. No external tool
(PyMOL, Gromacs, ...) is invoked.
"""

import subprocess
import sys
from argparse import _SubParsersAction

import pytest

REMOVED_AGENT_COMMANDS = {
    "agent_cancel",
    "agent_collect",
    "agent_plan",
    "agent_probe",
    "agent_profile_save",
    "agent_run",
    "agent_schema",
    "agent_status",
}
REMOVED_AGENT_OPTIONS = {"--json", "--dry-run", "--cluster-profile"}


def _subcommands(parser):
    return next(
        action for action in parser._actions if isinstance(action, _SubParsersAction)
    )


def test_cli_importable():
    """The cli module imports without error"""
    from mdtbx.cli import cli, create_parser

    assert callable(cli)
    assert callable(create_parser)


def test_agent_commands_and_options_are_removed():
    from mdtbx.cli import create_parser

    parser = create_parser()
    subcommands = _subcommands(parser)
    assert REMOVED_AGENT_COMMANDS.isdisjoint(subcommands.choices)

    option_strings = {
        option
        for command_parser in subcommands.choices.values()
        for action in command_parser._actions
        for option in action.option_strings
    }
    assert REMOVED_AGENT_OPTIONS.isdisjoint(option_strings)


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
