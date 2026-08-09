import subprocess

import pytest

from mdtbx.utils import tleap


def test_run_tleap_cleans_working_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fake_run_cmd(_command, **_kwargs):
        (tmp_path / "leap.log").write_text("log\n")

    monkeypatch.setattr(tleap, "run_cmd", fake_run_cmd)

    tleap.run_tleap("quit\n")

    assert not (tmp_path / "tleap.in").exists()
    assert not (tmp_path / "leap.log").exists()


def test_run_tleap_refuses_existing_working_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "leap.log").write_text("user data\n")

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        tleap.run_tleap("quit\n")

    assert (tmp_path / "leap.log").read_text() == "user data\n"


def test_run_tleap_cleans_after_command_failure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fail(_command, **_kwargs):
        (tmp_path / "leap.log").write_text("failed\n")
        raise subprocess.CalledProcessError(1, "tleap")

    monkeypatch.setattr(tleap, "run_cmd", fail)

    with pytest.raises(subprocess.CalledProcessError):
        tleap.run_tleap("quit\n")

    assert not (tmp_path / "tleap.in").exists()
    assert not (tmp_path / "leap.log").exists()
