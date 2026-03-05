"""
mod_mdp のユニットテスト

mod_mdp() 関数（ファイル I/O）と run()（ディレクトリ走査）をテストする。
tmp_path fixture により、テスト後に一時ファイルは自動削除される。
"""

from src.utils.mod_mdp import mod_mdp, run


class TestModMdp:
    def test_replace_existing_value(self, tmp_path):
        """既存キーの値が書き換わること"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps                  = 1000\n")

        mod_mdp("nsteps", "2000", mdp, ljust=23)

        content = mdp.read_text()
        assert "nsteps" in content
        assert "2000" in content
        assert "1000" not in content

    def test_add_new_key(self, tmp_path):
        """存在しないキーが末尾に追加されること"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps                  = 1000\n")

        mod_mdp("dt", "0.002", mdp, ljust=23)

        content = mdp.read_text()
        assert "dt" in content
        assert "0.002" in content
        # 元の内容は保持される
        assert "nsteps" in content

    def test_comment_lines_preserved(self, tmp_path):
        """コメント行（; で始まる）は変更されないこと"""
        mdp = tmp_path / "md.mdp"
        original = "; this is a comment\nnsteps = 1000\n"
        mdp.write_text(original)

        mod_mdp("nsteps", "2000", mdp, ljust=23)

        content = mdp.read_text()
        assert "; this is a comment" in content

    def test_ljust_controls_padding(self, tmp_path):
        """ljust パラメータが新規キーの整形幅を制御すること"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("")

        mod_mdp("dt", "0.002", mdp, ljust=30)

        content = mdp.read_text()
        # フォーマット: f"{key.ljust(30)} = {value}" → "dt" + 28空白 + " = 0.002"
        assert "dt" + " " * 28 + " = 0.002" in content

    def test_multiple_keys(self, tmp_path):
        """複数のキーを順番に書き換えられること"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps = 1000\ndt     = 0.002\n")

        mod_mdp("nsteps", "5000", mdp, ljust=7)
        mod_mdp("dt", "0.001", mdp, ljust=7)

        content = mdp.read_text()
        assert "5000" in content
        assert "0.001" in content


class TestModMdpRun:
    def test_run_modifies_all_mdp_files(self, tmp_path):
        """run() がディレクトリ内の全 .mdp ファイルを更新すること"""
        for name in ["em.mdp", "nvt.mdp", "npt.mdp"]:
            (tmp_path / name).write_text("nsteps = 100\n")

        import types

        args = types.SimpleNamespace(
            path=str(tmp_path),
            target_variable="nsteps",
            new_value="50000",
            exclude=None,
            ljust=23,
        )
        run(args)

        for name in ["em.mdp", "nvt.mdp", "npt.mdp"]:
            content = (tmp_path / name).read_text()
            assert "50000" in content

    def test_run_respects_exclude(self, tmp_path):
        """exclude オプションで指定したファイルは変更されないこと"""
        (tmp_path / "nvt.mdp").write_text("nsteps = 100\n")
        (tmp_path / "npt.mdp").write_text("nsteps = 100\n")

        import types

        args = types.SimpleNamespace(
            path=str(tmp_path),
            target_variable="nsteps",
            new_value="50000",
            exclude=["npt"],
            ljust=23,
        )
        run(args)

        assert "50000" in (tmp_path / "nvt.mdp").read_text()
        assert "100" in (tmp_path / "npt.mdp").read_text()  # 変更されていない

    def test_run_nonexistent_dir(self, tmp_path):
        """存在しないディレクトリを指定した場合でもエラーにならないこと（.mdp ファイルが 0 件）"""
        import types

        args = types.SimpleNamespace(
            path=str(tmp_path / "nodir"),
            target_variable="nsteps",
            new_value="1000",
            exclude=None,
            ljust=23,
        )
        # glob が空を返すだけでエラーにならない
        run(args)
