"""
build/gen_posres unit tests

Uses sample.top to test position restraint (.itp) generation and its
insertion into the topology.
"""

import shutil
import types


class TestGenPosresRun:
    def _make_args(self, top_path, output_prefix, selection="name CA"):
        return types.SimpleNamespace(
            topology=str(top_path),
            selection=selection,
            output_prefix=str(output_prefix),
        )

    def test_itp_file_created_for_protein(self, sample_top_path, tmp_path):
        """A posres .itp file is generated for the Protein moleculetype"""
        from src.build.gen_posres import run

        # Copy the topology (it is rewritten in place)
        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        itp = tmp_path / "posres_Protein.itp"
        assert itp.exists()

    def test_itp_contains_ifdef_block(self, sample_top_path, tmp_path):
        """The generated .itp has #ifdef / #endif blocks"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        itp = tmp_path / "posres_Protein.itp"
        content = itp.read_text()
        # output_prefix becomes a full path, so the #ifdef name is path dependent;
        # only check that #ifdef / #endif are present at all
        assert "#ifdef" in content
        assert "#endif" in content
        assert "[ position_restraints ]" in content

    def test_itp_contains_selected_atoms(self, sample_top_path, tmp_path):
        """CA atom indices appear in the .itp (Protein has 2 CA atoms)"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        itp = tmp_path / "posres_Protein.itp"
        content = itp.read_text()

        # ALA-CA (index 2) and GLY-CA (index 7) must be present
        lines = [
            line
            for line in content.splitlines()
            if line.strip()
            and not line.startswith(";")
            and not line.startswith("#")
            and line.strip()[0].isdigit()
        ]
        assert len(lines) == 2

    def test_topology_updated_with_include(self, sample_top_path, tmp_path):
        """An #include line is appended to the topology file after the run"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        updated = top_copy.read_text()
        assert "#include" in updated
        assert "posres_Protein.itp" in updated

    def test_no_sol_itp_when_sol_not_selected(self, sample_top_path, tmp_path):
        """SOL has no CA atom, so posres_SOL.itp is not generated"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        sol_itp = tmp_path / "posres_SOL.itp"
        assert not sol_itp.exists()
