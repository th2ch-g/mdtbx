"""
utils/partial_tempering unit tests

Uses sample.top to test appending an underscore to atom_type.
"""

import types


class TestPartialTemperingRun:
    def _make_args(self, topology_path, output_path, selection):
        return types.SimpleNamespace(
            topology=str(topology_path),
            selection=selection,
            output=str(output_path),
        )

    def test_output_file_created(self, sample_top_path, tmp_path):
        """The output file is generated"""
        from mdtbx.utils.partial_tempering import run

        out = tmp_path / "output.top"
        run(self._make_args(sample_top_path, out, selection="resname ALA"))
        assert out.exists()

    def test_selected_atoms_get_underscore(self, sample_top_path, tmp_path):
        """The atom_type of a selected ALA atom gains a trailing _"""
        from mdtbx.utils.partial_tempering import run

        out = tmp_path / "output.top"
        run(self._make_args(sample_top_path, out, selection="resname ALA"))

        content = out.read_text()
        # The ALA atom_type CT must have become CT_
        assert "CT_" in content

    def test_unselected_atoms_unchanged(self, sample_top_path, tmp_path):
        """The atom_type of an unselected GLY atom is unchanged"""
        from mdtbx.utils.partial_tempering import run

        out = tmp_path / "output.top"
        # Select ALA only
        run(self._make_args(sample_top_path, out, selection="resname ALA"))

        content = out.read_text()
        # The GLY lines (residue 2) stay as CT
        lines = content.splitlines()
        gly_lines = [
            line
            for line in lines
            if "GLY" in line and line.strip() and line.strip()[0].isdigit()
        ]
        for line in gly_lines:
            tokens = line.split()
            if len(tokens) >= 2:
                atom_type = tokens[1]
                assert not atom_type.endswith("_"), (
                    f"GLY atom_type should not have underscore: {line}"
                )

    def test_original_file_not_modified(self, sample_top_path, tmp_path):
        """The input file is left unchanged"""
        from mdtbx.utils.partial_tempering import run

        original_content = sample_top_path.read_text()
        out = tmp_path / "output.top"
        run(self._make_args(sample_top_path, out, selection="resname ALA"))

        assert sample_top_path.read_text() == original_content

    def test_no_match_selection_produces_unchanged_output(
        self, sample_top_path, tmp_path
    ):
        """A selection that matches nothing appends no underscore"""
        from mdtbx.utils.partial_tempering import run

        out = tmp_path / "output_nomatch.top"
        run(self._make_args(sample_top_path, out, selection="resname LIG"))

        output = out.read_text()
        # CT_ must not appear
        assert "CT_" not in output

    def test_second_application_is_idempotent(self, sample_top_path, tmp_path):
        from mdtbx.utils.partial_tempering import run

        first = tmp_path / "first.top"
        second = tmp_path / "second.top"
        run(self._make_args(sample_top_path, first, selection="resname ALA"))
        run(self._make_args(first, second, selection="resname ALA"))

        assert second.read_text() == first.read_text()
        assert "CT__" not in second.read_text()

    def test_nested_output_directory_is_created(self, sample_top_path, tmp_path):
        from mdtbx.utils.partial_tempering import run

        out = tmp_path / "nested" / "output.top"
        run(self._make_args(sample_top_path, out, selection="resname ALA"))

        assert out.exists()
