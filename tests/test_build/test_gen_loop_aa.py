import argparse

import pytest

from src.build import gen_loop_aa


def _parse_args(add_subcmd, argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    add_subcmd(subparsers)
    return parser.parse_args(argv)


def _fake_tleap_factory():
    """Capture the tleap input/keepfiles passed to run_tleap (which is mocked)."""
    calls = {}

    def fake_run_tleap(input_text, *, keepfiles=False, extra_cleanup=()):
        calls["input_text"] = input_text
        calls["keepfiles"] = keepfiles

    return calls, fake_run_tleap


def test_one_letter_sequence(monkeypatch):
    calls, fake = _fake_tleap_factory()
    monkeypatch.setattr(gen_loop_aa, "run_tleap", fake)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "AGS"])
    gen_loop_aa.run(args)

    assert "ALA GLY SER" in calls["input_text"]


def test_three_letter_sequence(monkeypatch):
    calls, fake = _fake_tleap_factory()
    monkeypatch.setattr(gen_loop_aa, "run_tleap", fake)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "ALA GLY SER"])
    gen_loop_aa.run(args)

    assert "ALA GLY SER" in calls["input_text"]


def test_cap_adds_ace_nme(monkeypatch):
    calls, fake = _fake_tleap_factory()
    monkeypatch.setattr(gen_loop_aa, "run_tleap", fake)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "--cap", "AGS"])
    gen_loop_aa.run(args)

    assert "ACE ALA GLY SER NME" in calls["input_text"]


def test_keepfiles_passes_through(monkeypatch):
    calls, fake = _fake_tleap_factory()
    monkeypatch.setattr(gen_loop_aa, "run_tleap", fake)

    args = _parse_args(gen_loop_aa.add_subcmd, ["gen_loop_aa", "--keepfiles", "AGS"])
    gen_loop_aa.run(args)

    assert calls["keepfiles"] is True


def test_custom_output(monkeypatch):
    calls, fake = _fake_tleap_factory()
    monkeypatch.setattr(gen_loop_aa, "run_tleap", fake)

    args = _parse_args(
        gen_loop_aa.add_subcmd, ["gen_loop_aa", "-o", "mypeptide.pdb", "AGS"]
    )
    gen_loop_aa.run(args)

    assert "mypeptide.pdb" in calls["input_text"]


def test_unknown_one_letter_raises():
    with pytest.raises(KeyError):
        gen_loop_aa._to_three_letter("AXG")


def test_all_twenty_amino_acids():
    result = gen_loop_aa._to_three_letter("ARNDCEQGHILKMFPSTWYV")
    assert len(result) == 20
    assert result[0] == "ALA"
    assert result[-1] == "VAL"
