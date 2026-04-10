"""
utils/partial_tempering のユニットテスト

sample.top を使って atom_type へのアンダースコア付与をテストする。
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
        """出力ファイルが生成されること"""
        from src.utils.partial_tempering import run

        out = tmp_path / "output.top"
        run(self._make_args(sample_top_path, out, selection="resname ALA"))
        assert out.exists()

    def test_selected_atoms_get_underscore(self, sample_top_path, tmp_path):
        """選択された ALA 原子の atom_type に _ が付加されること"""
        from src.utils.partial_tempering import run

        out = tmp_path / "output.top"
        run(self._make_args(sample_top_path, out, selection="resname ALA"))

        content = out.read_text()
        # ALA の atom_type CT が CT_ に変更されているはず
        assert "CT_" in content

    def test_unselected_atoms_unchanged(self, sample_top_path, tmp_path):
        """選択されていない GLY 原子の atom_type は変更されないこと"""
        from src.utils.partial_tempering import run

        out = tmp_path / "output.top"
        # ALA のみを選択
        run(self._make_args(sample_top_path, out, selection="resname ALA"))

        content = out.read_text()
        # GLY の行（residue 2）は CT のまま残る
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
        """入力ファイルが変更されないこと"""
        from src.utils.partial_tempering import run

        original_content = sample_top_path.read_text()
        out = tmp_path / "output.top"
        run(self._make_args(sample_top_path, out, selection="resname ALA"))

        assert sample_top_path.read_text() == original_content

    def test_no_match_selection_produces_unchanged_output(
        self, sample_top_path, tmp_path
    ):
        """マッチしない選択ではアンダースコアが付かないこと"""
        from src.utils.partial_tempering import run

        out = tmp_path / "output_nomatch.top"
        run(self._make_args(sample_top_path, out, selection="resname LIG"))

        output = out.read_text()
        # CT_ は含まれないはず
        assert "CT_" not in output
