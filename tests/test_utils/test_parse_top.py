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
        assert "moleculetype" in first

    def test_atoms_carry_their_moleculetype(self, parser):
        """各原子が所属する moleculetype 名を保持すること"""
        assert {a["moleculetype"] for a in parser.get_atoms_in("Protein")} == {
            "Protein"
        }
        assert {a["moleculetype"] for a in parser.get_atoms_in("SOL")} == {"SOL"}

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

    def test_section_headers_allow_inline_comments(self, tmp_path):
        topology = tmp_path / "inline-comment.top"
        topology.write_text(
            "[ moleculetype ] ; molecule definition\n"
            "Protein 3\n"
            "[ atoms ] ; atom records\n"
            "1 CT 1 ALA CA 1 0.0 12.0 ; alpha carbon\n"
            "[ system ] ; system name\n"
            "Example\n"
        )

        parsed = GromacsTopologyParser(str(topology))

        assert parsed.get_all_moleculetypes() == ["Protein"]
        assert parsed.get_atoms_in("Protein")[0]["name"] == "CA"

    def test_empty_topology_is_supported(self, tmp_path):
        topology = tmp_path / "empty.top"
        topology.write_text("")

        parsed = GromacsTopologyParser(str(topology))

        assert parsed.get_all_moleculetypes() == []

    def test_section_names_are_case_insensitive(self, tmp_path):
        topology = tmp_path / "uppercase.top"
        topology.write_text(
            "[ MOLECULETYPE ]\nProtein 3\n[ ATOMS ]\n1 CT 1 ALA CA 1 0.0 12.0\n"
        )

        parsed = GromacsTopologyParser(str(topology))

        assert parsed.get_atoms_in("Protein")[0]["name"] == "CA"

    def test_duplicate_moleculetype_raises(self, tmp_path):
        topology = tmp_path / "duplicate.top"
        topology.write_text(
            "[ moleculetype ]\nProtein 3\n[ moleculetype ]\nProtein 3\n"
        )

        with pytest.raises(ValueError, match="Duplicate moleculetype"):
            GromacsTopologyParser(str(topology))

    def test_malformed_atom_record_reports_line(self, tmp_path):
        topology = tmp_path / "malformed.top"
        topology.write_text("[ moleculetype ]\nProtein 3\n[ atoms ]\n1 CT\n")

        with pytest.raises(ValueError, match=r"malformed\.top:4"):
            GromacsTopologyParser(str(topology))

    def test_atoms_without_moleculetype_raises(self, tmp_path):
        topology = tmp_path / "orphan-atoms.top"
        topology.write_text("[ atoms ]\n1 CT 1 ALA CA\n")

        with pytest.raises(ValueError, match="no preceding moleculetype"):
            GromacsTopologyParser(str(topology))

    def test_moleculetype_without_declaration_raises(self, tmp_path):
        topology = tmp_path / "missing-name.top"
        topology.write_text("[ moleculetype ]\n")

        with pytest.raises(ValueError, match="declaration is missing"):
            GromacsTopologyParser(str(topology))
