"""
atom_selection_parser のユニットテスト

SelectionParser / AtomSelector の純粋なパースロジックをテストする。
外部ツール・ファイル I/O に依存しない。
"""

import pytest

from src.utils.atom_selection_parser import (
    All,
    And,
    AtomSelector,
    Backbone,
    Bracket,
    Chain,
    Index,
    Name,
    Not,
    Or,
    Protein,
    ResId,
    ResName,
    SelectionParser,
    Sidechain,
    Water,
    parse_selection,
)

# ---------------------------------------------------------------------------
# ヘルパー: テスト用の原子辞書
# ---------------------------------------------------------------------------

ALA_CA = {"resname": "ALA", "name": "CA", "resid": 1, "index": 1}
ALA_CB = {"resname": "ALA", "name": "CB", "resid": 1, "index": 2}
GLY_N = {"resname": "GLY", "name": "N", "resid": 2, "index": 3}
HOH_O = {"resname": "HOH", "name": "O", "resid": 3, "index": 4}
NA_ION = {"resname": "NA", "name": "NA", "resid": 4, "index": 5}

# ---------------------------------------------------------------------------
# SelectionNode の基本評価
# ---------------------------------------------------------------------------


class TestAll:
    def test_always_true(self):
        node = All()
        assert node.eval(ALA_CA)
        assert node.eval(HOH_O)
        assert node.eval({})


class TestProtein:
    def test_ala_is_protein(self):
        assert Protein().eval(ALA_CA)

    def test_water_is_not_protein(self):
        assert not Protein().eval(HOH_O)

    def test_unknown_resname_is_not_protein(self):
        assert not Protein().eval({"resname": "LIG", "name": "C1"})


class TestWater:
    def test_hoh_is_water(self):
        assert Water().eval(HOH_O)

    def test_sol_is_water(self):
        assert Water().eval({"resname": "SOL", "name": "OW"})

    def test_protein_is_not_water(self):
        assert not Water().eval(ALA_CA)


class TestBackbone:
    def test_ca_is_backbone(self):
        assert Backbone().eval(ALA_CA)

    def test_cb_is_not_backbone(self):
        assert not Backbone().eval(ALA_CB)

    def test_n_is_backbone(self):
        assert Backbone().eval({"resname": "ALA", "name": "N"})

    def test_water_o_is_not_backbone(self):
        assert not Backbone().eval(HOH_O)


class TestSidechain:
    def test_cb_is_sidechain(self):
        assert Sidechain().eval(ALA_CB)

    def test_ca_is_not_sidechain(self):
        assert not Sidechain().eval(ALA_CA)

    def test_water_is_not_sidechain(self):
        assert not Sidechain().eval(HOH_O)


class TestResName:
    def test_match_single(self):
        assert ResName(["ALA"]).eval(ALA_CA)

    def test_match_multiple(self):
        assert ResName(["ALA", "GLY"]).eval(GLY_N)

    def test_no_match(self):
        assert not ResName(["ALA"]).eval(HOH_O)

    def test_wildcard(self):
        assert ResName(["A*"]).eval(ALA_CA)
        assert not ResName(["A*"]).eval(GLY_N)


class TestResId:
    def test_match(self):
        assert ResId([1]).eval(ALA_CA)

    def test_range(self):
        assert ResId([1, 2, 3]).eval(GLY_N)

    def test_no_match(self):
        assert not ResId([5, 6]).eval(ALA_CA)


class TestName:
    def test_match(self):
        assert Name(["CA"]).eval(ALA_CA)

    def test_wildcard_prefix(self):
        assert Name(["C*"]).eval(ALA_CA)
        assert Name(["C*"]).eval(ALA_CB)
        assert not Name(["C*"]).eval(GLY_N)

    def test_wildcard_all(self):
        assert Name(["*"]).eval(ALA_CA)


class TestIndex:
    def test_match(self):
        assert Index([1]).eval(ALA_CA)

    def test_no_match(self):
        assert not Index([99]).eval(ALA_CA)


class TestNot:
    def test_inverts_protein(self):
        assert Not(Protein()).eval(HOH_O)
        assert not Not(Protein()).eval(ALA_CA)


class TestAnd:
    def test_both_true(self):
        assert And([Protein(), Backbone()]).eval(ALA_CA)

    def test_one_false(self):
        assert not And([Protein(), Backbone()]).eval(ALA_CB)

    def test_empty_is_true(self):
        assert And([]).eval(ALA_CA)


class TestOr:
    def test_one_true(self):
        assert Or([Protein(), Water()]).eval(HOH_O)
        assert Or([Protein(), Water()]).eval(ALA_CA)

    def test_both_false(self):
        assert not Or([Protein(), Water()]).eval(NA_ION)

    def test_empty_is_false(self):
        assert not Or([]).eval(ALA_CA)


# ---------------------------------------------------------------------------
# SelectionParser パース結果の検証
# ---------------------------------------------------------------------------


class TestSelectionParser:
    def test_protein(self):
        result = SelectionParser("protein").parse()
        assert isinstance(result, Protein)

    def test_all(self):
        result = SelectionParser("all").parse()
        assert isinstance(result, All)

    def test_resname_single(self):
        result = SelectionParser("resname ALA").parse()
        assert isinstance(result, ResName)
        assert result.names == ["ALA"]

    def test_resname_multiple(self):
        result = SelectionParser("resname ALA GLY").parse()
        assert isinstance(result, ResName)
        assert result.names == ["ALA", "GLY"]

    def test_resid_single(self):
        result = SelectionParser("resid 5").parse()
        assert isinstance(result, ResId)
        assert result.ids == [5]

    def test_resid_range(self):
        result = SelectionParser("resid 1 to 5").parse()
        assert isinstance(result, ResId)
        assert result.ids == [1, 2, 3, 4, 5]

    def test_and_expression(self):
        result = SelectionParser("protein and backbone").parse()
        assert isinstance(result, And)

    def test_or_expression(self):
        result = SelectionParser("protein or water").parse()
        assert isinstance(result, Or)

    def test_not_expression(self):
        result = SelectionParser("not water").parse()
        assert isinstance(result, Not)

    def test_bracket(self):
        result = SelectionParser("(protein and backbone)").parse()
        assert isinstance(result, Bracket)

    def test_chain(self):
        result = SelectionParser("chain A").parse()
        assert isinstance(result, Chain)

    def test_name_wildcard(self):
        result = SelectionParser("name C*").parse()
        assert isinstance(result, Name)

    def test_index_range(self):
        result = SelectionParser("index 1 to 3").parse()
        assert isinstance(result, Index)
        assert result.indices == [1, 2, 3]


# ---------------------------------------------------------------------------
# AtomSelector の end-to-end テスト
# ---------------------------------------------------------------------------


class TestAtomSelector:
    @pytest.mark.parametrize(
        "selection, mol, expected",
        [
            ("protein and backbone", ALA_CA, True),
            ("protein and backbone", ALA_CB, False),
            ("protein and sidechain", ALA_CB, True),
            ("protein and sidechain", HOH_O, False),
            ("water", HOH_O, True),
            ("not water", HOH_O, False),
            ("resid 1 to 3 and name CA", ALA_CA, True),
            ("resid 1 to 3 and name CA", GLY_N, False),
            ("all", {"resname": "XYZ"}, True),
            ("name C*", ALA_CA, True),
            ("name C*", ALA_CB, True),
            ("name C*", GLY_N, False),
            ("resname A*", ALA_CA, True),
            ("resname A*", HOH_O, False),
            ("(resname ALA GLY) or (resname HOH and name O)", ALA_CA, True),
            ("(resname ALA GLY) or (resname HOH and name O)", HOH_O, True),
            ("(resname ALA GLY) or (resname HOH and name O)", NA_ION, False),
        ],
    )
    def test_selection_result(self, selection, mol, expected):
        selector = AtomSelector(selection)
        assert selector.eval(mol) == expected

    def test_double_not(self):
        selector = AtomSelector("not not protein")
        assert selector.eval(ALA_CA)

    def test_chain_selection(self):
        mol_a = {"resname": "ALA", "name": "CA", "chain": "A"}
        mol_b = {"resname": "ALA", "name": "CA", "chain": "B"}
        assert AtomSelector("chain A and protein").eval(mol_a)
        assert not AtomSelector("chain B and protein").eval(mol_a)
        assert AtomSelector("chain B and protein").eval(mol_b)


# ---------------------------------------------------------------------------
# エラー処理
# ---------------------------------------------------------------------------


class TestParseErrors:
    @pytest.mark.parametrize(
        "bad_selection",
        [
            "resid 10 to",
            "resname and",
            "( resid 1",
            "name ca and",
        ],
    )
    def test_invalid_selection_raises(self, bad_selection):
        with pytest.raises(ValueError):
            AtomSelector(bad_selection)

    def test_parse_selection_returns_error_string_on_failure(self):
        result = parse_selection("resid 10 to")
        assert isinstance(result, str)
