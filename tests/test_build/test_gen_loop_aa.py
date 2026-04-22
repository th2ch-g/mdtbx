import argparse

import pytest

from src.build import gen_loop_aa


def _parse_args(add_subcmd, argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    add_subcmd(subparsers)
    return parser.parse_args(argv)


def _fake_run_factory():
    commands = []

    def fake_run(command, shell, check):
        commands.append(command)

    return commands, fake_run


def test_one_letter_sequence(tmp_path, monkeypatch):
    commands, fake_run = _fake_run_factory()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gen_loop_aa.subprocess, "run", fake_run)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "AGS"])
    gen_loop_aa.run(args)

    tleap_in = (tmp_path / "tleap.in").read_text()
    assert "ALA GLY SER" in tleap_in
    assert commands[0] == "tleap -f tleap.in"


def test_three_letter_sequence(tmp_path, monkeypatch):
    commands, fake_run = _fake_run_factory()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gen_loop_aa.subprocess, "run", fake_run)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "ALA GLY SER"])
    gen_loop_aa.run(args)

    tleap_in = (tmp_path / "tleap.in").read_text()
    assert "ALA GLY SER" in tleap_in


def test_cap_adds_ace_nme(tmp_path, monkeypatch):
    commands, fake_run = _fake_run_factory()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gen_loop_aa.subprocess, "run", fake_run)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "--cap", "AGS"])
    gen_loop_aa.run(args)

    tleap_in = (tmp_path / "tleap.in").read_text()
    assert "ACE ALA GLY SER NME" in tleap_in


def test_keepfiles_skips_cleanup(tmp_path, monkeypatch):
    commands, fake_run = _fake_run_factory()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gen_loop_aa.subprocess, "run", fake_run)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "--keepfiles", "AGS"])
    gen_loop_aa.run(args)

    # Only tleap, no rm
    assert len(commands) == 1
    assert commands[0] == "tleap -f tleap.in"


def test_custom_output(tmp_path, monkeypatch):
    commands, fake_run = _fake_run_factory()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gen_loop_aa.subprocess, "run", fake_run)

    args = _parse_args(
        gen_loop_aa.add_subcmd, ["gen_loop_aa", "-o", "mypeptide.pdb", "AGS"]
    )
    gen_loop_aa.run(args)

    tleap_in = (tmp_path / "tleap.in").read_text()
    assert "mypeptide.pdb" in tleap_in


def test_unknown_one_letter_raises():
    with pytest.raises(KeyError):
        gen_loop_aa._to_three_letter("AXG")


def test_all_twenty_amino_acids():
    result = gen_loop_aa._to_three_letter("ARNDCEQGHILKMFPSTWYV")
    assert len(result) == 20
    assert result[0] == "ALA"
    assert result[-1] == "VAL"
