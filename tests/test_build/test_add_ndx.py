"""
build/add_ndx のユニットテスト

gmx make_ndx に依存しない純粋な count_index_group() ヘルパーをテストする。
"""

import types
from pathlib import Path


class TestCountIndexGroup:
    def _make_args(self, index_path):
        return types.SimpleNamespace(index=str(index_path))

    def _write_ndx(self, path: Path, n_groups: int) -> Path:
        """n_groups 個のインデックスグループを持つ .ndx ファイルを作成する"""
        lines = []
        for i in range(n_groups):
            lines.append(f"[ Group{i} ]\n")
            lines.append(f"{i + 1}\n")
        path.write_text("".join(lines))
        return path

    def test_counts_single_group(self, tmp_path):
        """グループが 1 つの .ndx を正しくカウントすること"""
        from src.build.add_ndx import count_index_group

        ndx = self._write_ndx(tmp_path / "single.ndx", n_groups=1)
        assert count_index_group(self._make_args(ndx)) == 1

    def test_counts_multiple_groups(self, tmp_path):
        """複数グループの .ndx を正しくカウントすること"""
        from src.build.add_ndx import count_index_group

        ndx = self._write_ndx(tmp_path / "multi.ndx", n_groups=5)
        assert count_index_group(self._make_args(ndx)) == 5

    def test_ignores_non_bracket_lines(self, tmp_path):
        """[ で始まらない行はカウントしないこと"""
        from src.build.add_ndx import count_index_group

        ndx = tmp_path / "mixed.ndx"
        ndx.write_text(
            "[ System ]\n"
            "1 2 3 4 5\n"
            "[ Protein ]\n"
            "1 2 3\n"
            "Some comment line\n"
            "[ Water ]\n"
            "4 5\n"
        )
        assert count_index_group(self._make_args(ndx)) == 3

    def test_empty_file_returns_zero(self, tmp_path):
        """空ファイルは 0 を返すこと"""
        from src.build.add_ndx import count_index_group

        ndx = tmp_path / "empty.ndx"
        ndx.write_text("")
        assert count_index_group(self._make_args(ndx)) == 0

    def test_real_gromacs_style_ndx(self, tmp_path):
        """GROMACS 形式の .ndx を正しくカウントすること"""
        from src.build.add_ndx import count_index_group

        ndx = tmp_path / "gromacs.ndx"
        ndx.write_text(
            "[ System ]\n"
            "   1    2    3    4    5    6    7    8    9\n"
            "[ Protein ]\n"
            "   1    2    3    4    5\n"
            "[ non-Protein ]\n"
            "   6    7    8    9\n"
        )
        assert count_index_group(self._make_args(ndx)) == 3
