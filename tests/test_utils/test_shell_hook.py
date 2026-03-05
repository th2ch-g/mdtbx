"""
utils/shell_hook のユニットテスト

run() が適切なシェルスクリプトのテンプレートを標準出力に出力することをテストする。
"""

import types


class TestShellHookRun:
    def _make_args(self):
        return types.SimpleNamespace()

    def test_output_contains_mdtbx_alias(self, capsys):
        """mdtbx の alias / function 定義が含まれること"""
        from src.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "mdtbx" in captured.out

    def test_output_contains_pymol(self, capsys):
        """pymol の設定が含まれること"""
        from src.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "pymol" in captured.out

    def test_output_contains_begin_end_markers(self, capsys):
        """BEGIN / END マーカーが含まれること"""
        from src.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "BEGIN OF MDTBX SHELL HOOK" in captured.out
        assert "END OF MDTBX SHELL HOOK" in captured.out

    def test_output_contains_path_export(self, capsys):
        """PATH 設定が含まれること"""
        from src.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "PATH" in captured.out

    def test_output_is_nonempty(self, capsys):
        """出力が空でないこと"""
        from src.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert len(captured.out.strip()) > 0
