"""
build/gen_posres のユニットテスト

sample.top を使って位置拘束 (.itp) の生成と topology への挿入をテストする。
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
        """Protein モジュールの posres .itp ファイルが生成されること"""
        from src.build.gen_posres import run

        # topology をコピー（in-place で書き換えられるため）
        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        itp = tmp_path / "posres_Protein.itp"
        assert itp.exists()

    def test_itp_contains_ifdef_block(self, sample_top_path, tmp_path):
        """生成された .itp が #ifdef / #endif ブロックを持つこと"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        itp = tmp_path / "posres_Protein.itp"
        content = itp.read_text()
        # output_prefix がフルパスになるため #ifdef の名前はパス依存
        # 形式的な #ifdef / #endif の存在を確認する
        assert "#ifdef" in content
        assert "#endif" in content
        assert "[ position_restraints ]" in content

    def test_itp_contains_selected_atoms(self, sample_top_path, tmp_path):
        """CA 原子のインデックスが .itp に含まれること（Protein に CA は 2 個）"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        itp = tmp_path / "posres_Protein.itp"
        content = itp.read_text()

        # ALA-CA (index 2) と GLY-CA (index 7) が含まれること
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
        """実行後の topology ファイルに #include 行が追加されること"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        updated = top_copy.read_text()
        assert "#include" in updated
        assert "posres_Protein.itp" in updated

    def test_no_sol_itp_when_sol_not_selected(self, sample_top_path, tmp_path):
        """SOL に CA 原子がないため posres_SOL.itp は生成されないこと"""
        from src.build.gen_posres import run

        top_copy = tmp_path / "test.top"
        shutil.copy(str(sample_top_path), str(top_copy))

        prefix = tmp_path / "posres"
        run(self._make_args(top_copy, prefix, selection="name CA"))

        sol_itp = tmp_path / "posres_SOL.itp"
        assert not sol_itp.exists()
