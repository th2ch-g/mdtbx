from pathlib import Path

from mdtbx.utils.gmx import gmx_index_args, gmx_tempfile, remove_gmx_backups


def test_gmx_index_args():
    assert gmx_index_args(None) == []
    assert gmx_index_args("index with spaces.ndx") == [
        "-n",
        "index with spaces.ndx",
    ]


def test_gmx_tempfile_reserves_unique_missing_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with gmx_tempfile(".xvg") as path:
        temporary_path = Path(path)
        assert not temporary_path.exists()
        temporary_path.write_text("data")

    assert not temporary_path.exists()


def test_remove_gmx_backups_only_removes_backup_pattern(tmp_path):
    backup = tmp_path / "#trajectory.xtc.1#"
    unrelated = tmp_path / "#notes"
    backup.touch()
    unrelated.touch()

    remove_gmx_backups(tmp_path)

    assert not backup.exists()
    assert unrelated.exists()
