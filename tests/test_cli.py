"""
CLI のサブコマンド登録テスト

全サブコマンドが argparse に正常に登録されることを確認する。
外部ツール(PyMOL, Gromacs 等)の呼び出しは発生しない。
"""

import pytest


def test_cli_importable():
    """cli モジュールが import エラーなくロードできること"""
    from src.cli import cli

    assert callable(cli)


def test_all_subcommands_registered():
    """
    全サブコマンドが argparse に登録されており、
    --help が正常に動作すること（SystemExit(0) で終了）
    """
    from src.cli import cli

    with pytest.raises(SystemExit) as exc_info:
        # --help は SystemExit(0) を発生させる
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
    ],
)
def test_subcommand_help(subcmd):
    """各サブコマンドの --help が正常終了すること"""
    import sys

    from src.cli import cli

    sys.argv = ["mdtbx", subcmd, "--help"]
    with pytest.raises(SystemExit) as exc_info:
        cli()
    assert exc_info.value.code == 0
