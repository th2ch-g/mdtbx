"""
rmfile のユニットテスト

対象パターン（#*#, *cpt, mdout.mdp）のファイルが削除され、
それ以外のファイルは残ることを確認する。
"""

import types


from src.utils.rmfile import run


class TestRmfile:
    def _make_args(self, path):
        return types.SimpleNamespace(path=str(path))

    def test_removes_cpt_files(self, tmp_path):
        cpt = tmp_path / "state.cpt"
        cpt.touch()
        run(self._make_args(tmp_path))
        assert not cpt.exists()

    def test_removes_mdout_mdp(self, tmp_path):
        mdout = tmp_path / "mdout.mdp"
        mdout.touch()
        run(self._make_args(tmp_path))
        assert not mdout.exists()

    def test_removes_backup_files(self, tmp_path):
        backup = tmp_path / "#step.trr.1#"
        backup.touch()
        run(self._make_args(tmp_path))
        assert not backup.exists()

    def test_preserves_other_files(self, tmp_path):
        mdp = tmp_path / "md.mdp"
        top = tmp_path / "topol.top"
        tpr = tmp_path / "topol.tpr"
        for f in [mdp, top, tpr]:
            f.touch()

        run(self._make_args(tmp_path))

        assert mdp.exists()
        assert top.exists()
        assert tpr.exists()

    def test_removes_cpt_recursively(self, tmp_path):
        """サブディレクトリ内の .cpt も再帰的に削除されること"""
        subdir = tmp_path / "run1"
        subdir.mkdir()
        cpt = subdir / "state.cpt"
        cpt.touch()

        run(self._make_args(tmp_path))

        assert not cpt.exists()

    def test_mixed_files(self, tmp_path):
        """削除対象と保持対象が混在していても正しく処理されること"""
        targets = [
            tmp_path / "state.cpt",
            tmp_path / "mdout.mdp",
            tmp_path / "#backup#",
        ]
        keepfiles = [
            tmp_path / "topol.top",
            tmp_path / "md.mdp",
        ]
        for f in targets + keepfiles:
            f.touch()

        run(self._make_args(tmp_path))

        for f in targets:
            assert not f.exists(), f"{f.name} should have been deleted"
        for f in keepfiles:
            assert f.exists(), f"{f.name} should be preserved"
