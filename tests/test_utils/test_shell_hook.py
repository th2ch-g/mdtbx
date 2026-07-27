"""
utils/shell_hook unit tests

Tests that run() prints a suitable shell script template to stdout.
"""

import types


class TestShellHookRun:
    def _make_args(self):
        return types.SimpleNamespace()

    def test_output_contains_mdtbx_alias(self, capsys):
        """The mdtbx alias / function definitions are present"""
        from mdtbx.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "mdtbx" in captured.out

    def test_output_contains_pymol(self, capsys):
        """The pymol settings are present"""
        from mdtbx.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "pymol" in captured.out

    def test_output_contains_begin_end_markers(self, capsys):
        """The BEGIN / END markers are present"""
        from mdtbx.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "BEGIN OF MDTBX SHELL HOOK" in captured.out
        assert "END OF MDTBX SHELL HOOK" in captured.out

    def test_output_contains_path_export(self, capsys):
        """The PATH export is present"""
        from mdtbx.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert "PATH" in captured.out

    def test_output_is_nonempty(self, capsys):
        """The output is not empty"""
        from mdtbx.utils.shell_hook import run

        run(self._make_args())
        captured = capsys.readouterr()
        assert len(captured.out.strip()) > 0
