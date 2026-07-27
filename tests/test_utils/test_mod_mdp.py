"""
mod_mdp unit tests

Tests the mod_mdp() function (file I/O) and run() (directory traversal).
The tmp_path fixture removes the temporary files after each test.
"""

from src.utils.mod_mdp import mod_mdp, run


class TestModMdp:
    def test_replace_existing_value(self, tmp_path):
        """The value of an existing key is rewritten"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps                  = 1000\n")

        mod_mdp("nsteps", "2000", mdp, ljust=23)

        content = mdp.read_text()
        assert "nsteps" in content
        assert "2000" in content
        assert "1000" not in content

    def test_add_new_key(self, tmp_path):
        """A key that does not exist is appended at the end"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps                  = 1000\n")

        mod_mdp("dt", "0.002", mdp, ljust=23)

        content = mdp.read_text()
        assert "dt" in content
        assert "0.002" in content
        # The original content is preserved
        assert "nsteps" in content

    def test_comment_lines_preserved(self, tmp_path):
        """Comment lines (starting with ;) are left unchanged"""
        mdp = tmp_path / "md.mdp"
        original = "; this is a comment\nnsteps = 1000\n"
        mdp.write_text(original)

        mod_mdp("nsteps", "2000", mdp, ljust=23)

        content = mdp.read_text()
        assert "; this is a comment" in content

    def test_ljust_controls_padding(self, tmp_path):
        """The ljust parameter controls the padding width of a new key"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("")

        mod_mdp("dt", "0.002", mdp, ljust=30)

        content = mdp.read_text()
        # Format: f"{key.ljust(30)} = {value}" -> "dt" + 28 spaces + " = 0.002"
        assert "dt" + " " * 28 + " = 0.002" in content

    def test_multiple_keys(self, tmp_path):
        """Several keys can be rewritten in sequence"""
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps = 1000\ndt     = 0.002\n")

        mod_mdp("nsteps", "5000", mdp, ljust=7)
        mod_mdp("dt", "0.001", mdp, ljust=7)

        content = mdp.read_text()
        assert "5000" in content
        assert "0.001" in content

    def test_add_key_after_line_without_newline(self, tmp_path):
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps = 1000")

        mod_mdp("dt", "0.002", mdp, ljust=7)

        assert mdp.read_text() == "nsteps = 1000\ndt      = 0.002\n"

    def test_preserves_inline_comment(self, tmp_path):
        mdp = tmp_path / "md.mdp"
        mdp.write_text("nsteps = 1000 ; production length\n")

        mod_mdp("nsteps", "2000", mdp, ljust=23)

        assert mdp.read_text() == "nsteps = 2000 ; production length\n"

    def test_rejects_invalid_key(self, tmp_path):
        import pytest

        mdp = tmp_path / "md.mdp"
        mdp.write_text("")

        with pytest.raises(ValueError, match="MDP key"):
            mod_mdp("bad=key", "1", mdp, ljust=23)

    def test_rejects_multiline_value(self, tmp_path):
        import pytest

        mdp = tmp_path / "md.mdp"
        mdp.write_text("")

        with pytest.raises(ValueError, match="single line"):
            mod_mdp("nsteps", "1\nintegrator = md", mdp, ljust=23)


class TestModMdpRun:
    def test_run_modifies_all_mdp_files(self, tmp_path):
        """run() updates every .mdp file in the directory"""
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
        """Files named by the exclude option are left unchanged"""
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
        assert "100" in (tmp_path / "npt.mdp").read_text()  # unchanged

    def test_run_nonexistent_dir(self, tmp_path):
        """A directory that does not exist is not an error (it holds 0 .mdp files)"""
        import types

        args = types.SimpleNamespace(
            path=str(tmp_path / "nodir"),
            target_variable="nsteps",
            new_value="1000",
            exclude=None,
            ljust=23,
        )
        # glob just returns empty rather than raising
        run(args)
