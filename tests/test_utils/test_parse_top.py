"""
parse_top (GromacsTopologyParser) のユニットテスト

tests/fixtures/sample.top を使ってトポロジーパーサーをテストする。
"""

import pytest

from src.utils.parse_top import GromacsTopologyParser


@pytest.fixture(scope="module")
def parser(sample_top_path):
    return GromacsTopologyParser(str(sample_top_path))


class TestGromacsTopologyParser:
    def test_get_all_moleculetypes(self, parser):
        """モジュール名が正しく取得できること"""
        moltypes = parser.get_all_moleculetypes()
        assert "Protein" in moltypes
        assert "SOL" in moltypes

    def test_moleculetype_order(self, parser):
        """定義順が保持されること"""
        moltypes = parser.get_all_moleculetypes()
        assert moltypes.index("Protein") < moltypes.index("SOL")

    def test_get_atoms_in_protein(self, parser):
        """Protein の原子リストが取得できること"""
        atoms = parser.get_atoms_in("Protein")
        assert len(atoms) > 0

    def test_protein_atom_fields(self, parser):
        """原子辞書に必須フィールドが含まれること"""
        atoms = parser.get_atoms_in("Protein")
        first = atoms[0]
        assert "atom_type" in first
        assert "index" in first
        assert "resid" in first
        assert "resname" in first
        assert "name" in first

    def test_protein_residue_names(self, parser):
        """ALA と GLY の残基が含まれること（fixture の内容に対応）"""
        atoms = parser.get_atoms_in("Protein")
        resnames = {a["resname"] for a in atoms}
        assert "ALA" in resnames
        assert "GLY" in resnames

    def test_get_atoms_in_sol(self, parser):
        """SOL の原子リストが取得できること"""
        atoms = parser.get_atoms_in("SOL")
        assert len(atoms) == 3  # OW, HW1, HW2

    def test_sol_atom_names(self, parser):
        """SOL の原子名が正しいこと"""
        atoms = parser.get_atoms_in("SOL")
        names = [a["name"] for a in atoms]
        assert "OW" in names
        assert "HW1" in names
        assert "HW2" in names

    def test_get_insert_linenumber(self, parser):
        """挿入行番号が整数で返ること"""
        lineno = parser.get_insert_linenumber_in("Protein")
        assert isinstance(lineno, int)
        assert lineno > 0

    def test_atom_index_sequential(self, parser):
        """原子のインデックスが連続していること"""
        atoms = parser.get_atoms_in("Protein")
        indices = [a["index"] for a in atoms]
        assert indices == list(range(1, len(atoms) + 1))

    def test_invalid_moleculetype_raises(self, parser):
        """存在しないモジュール名は KeyError になること"""
        with pytest.raises(KeyError):
            parser.get_atoms_in("NONEXISTENT")
