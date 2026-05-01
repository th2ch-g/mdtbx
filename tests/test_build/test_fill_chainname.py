"""
build/fill_chainname のユニットテスト

PyMOL の load/iterate/alter/save をモックして、
ブロック分割と chain ID 割当ロジックを検証する。
"""

import types
from unittest.mock import MagicMock

import pytest


class TestSplitIntoBlocks:
    def test_single_continuous_block(self):
        from src.build.fill_chainname import _split_into_blocks

        atoms = [(1, 1, ""), (2, 1, ""), (3, 2, ""), (4, 2, "")]
        blocks = _split_into_blocks(atoms, max_resi_gap=1, use_segi=True)

        assert len(blocks) == 1
        assert blocks[0] == [1, 2, 3, 4]

    def test_split_on_residue_gap(self):
        from src.build.fill_chainname import _split_into_blocks

        # resv jumps 2 -> 100, gap of 98 should start a new block
        atoms = [(1, 1, ""), (2, 2, ""), (3, 100, ""), (4, 101, "")]
        blocks = _split_into_blocks(atoms, max_resi_gap=1, use_segi=True)

        assert len(blocks) == 2
        assert blocks[0] == [1, 2]
        assert blocks[1] == [3, 4]

    def test_split_on_segi_change(self):
        from src.build.fill_chainname import _split_into_blocks

        atoms = [(1, 1, "P1"), (2, 2, "P1"), (3, 3, "P2"), (4, 4, "P2")]
        blocks = _split_into_blocks(atoms, max_resi_gap=1, use_segi=True)

        assert len(blocks) == 2
        assert blocks[0] == [1, 2]
        assert blocks[1] == [3, 4]

    def test_ignore_segi_when_use_segi_false(self):
        from src.build.fill_chainname import _split_into_blocks

        # Same resv flow but segi differs; with use_segi=False stay one block
        atoms = [(1, 1, "P1"), (2, 2, "P2"), (3, 3, "P3")]
        blocks = _split_into_blocks(atoms, max_resi_gap=1, use_segi=False)

        assert len(blocks) == 1

    def test_empty_input(self):
        from src.build.fill_chainname import _split_into_blocks

        assert _split_into_blocks([], max_resi_gap=1, use_segi=True) == []


class TestNextChainId:
    def test_returns_first_unused_uppercase(self):
        from src.build.fill_chainname import _next_chain_id

        assert _next_chain_id(set()) == "A"
        assert _next_chain_id({"A"}) == "B"
        assert _next_chain_id({"A", "B", "C"}) == "D"

    def test_skips_used_then_falls_back_to_lowercase(self):
        from src.build.fill_chainname import _next_chain_id

        used = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        assert _next_chain_id(used) == "a"

    def test_raises_when_pool_exhausted(self):
        from src.build.fill_chainname import _next_chain_id, CHAIN_POOL

        with pytest.raises(RuntimeError, match="No single-letter chain ID"):
            _next_chain_id(set(CHAIN_POOL))


class TestRun:
    def _make_args(self, structure, output, **overrides):
        args = types.SimpleNamespace(
            structure=str(structure),
            output=str(output),
            max_resi_gap=1,
            ignore_segi=False,
            chainname=None,
        )
        for k, v in overrides.items():
            setattr(args, k, v)
        return args

    def test_no_blank_chain_saves_as_is(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")
        output = tmp_path / "out.pdb"

        pymol_cmd = MagicMock()
        pymol_cmd.get_chains.return_value = ["A", "B"]
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(fill_chainname, "_collect_blank_atoms", lambda _: [])

        fill_chainname.run(self._make_args(structure, output))

        pymol_cmd.alter.assert_not_called()
        pymol_cmd.save.assert_called_once_with(str(output), "target")

    def test_assigns_next_unused_chain_id(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")
        output = tmp_path / "out.pdb"

        pymol_cmd = MagicMock()
        # existing: A and C, so blanks should become B then D
        pymol_cmd.get_chains.return_value = ["A", "C", ""]
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(
            fill_chainname,
            "_collect_blank_atoms",
            lambda _: [(1, 1, "S1"), (2, 2, "S1"), (3, 100, "S2")],
        )

        fill_chainname.run(self._make_args(structure, output))

        # Two blocks => two alter calls with chain='B' then chain='D'
        assert pymol_cmd.alter.call_count == 2
        first_expr = pymol_cmd.alter.call_args_list[0].args[1]
        second_expr = pymol_cmd.alter.call_args_list[1].args[1]
        assert first_expr == "chain='B'"
        assert second_expr == "chain='D'"
        pymol_cmd.sort.assert_called_once()
        pymol_cmd.save.assert_called_once_with(str(output), "target")

    def test_user_specified_chainnames_are_used(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")
        output = tmp_path / "out.pdb"

        pymol_cmd = MagicMock()
        pymol_cmd.get_chains.return_value = ["A"]
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(
            fill_chainname,
            "_collect_blank_atoms",
            lambda _: [(1, 1, "S1"), (2, 100, "S2")],
        )

        args = self._make_args(structure, output, chainname=["X", "Y"])
        fill_chainname.run(args)

        exprs = [c.args[1] for c in pymol_cmd.alter.call_args_list]
        assert exprs == ["chain='X'", "chain='Y'"]

    def test_user_chainname_count_mismatch_raises(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")

        pymol_cmd = MagicMock()
        pymol_cmd.get_chains.return_value = ["A"]
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(
            fill_chainname,
            "_collect_blank_atoms",
            lambda _: [(1, 1, "S1"), (2, 100, "S2")],
        )

        args = self._make_args(structure, tmp_path / "out.pdb", chainname=["B"])
        with pytest.raises(ValueError, match="2 blank-chain block"):
            fill_chainname.run(args)

    def test_user_chainname_collision_with_existing_raises(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")

        pymol_cmd = MagicMock()
        pymol_cmd.get_chains.return_value = ["A", "B"]
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(
            fill_chainname, "_collect_blank_atoms", lambda _: [(1, 1, "S1")]
        )

        args = self._make_args(structure, tmp_path / "out.pdb", chainname=["B"])
        with pytest.raises(ValueError, match="already exists"):
            fill_chainname.run(args)

    def test_user_chainname_invalid_letter_raises(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")

        pymol_cmd = MagicMock()
        pymol_cmd.get_chains.return_value = []
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(
            fill_chainname, "_collect_blank_atoms", lambda _: [(1, 1, "S1")]
        )

        args = self._make_args(structure, tmp_path / "out.pdb", chainname=["XX"])
        with pytest.raises(ValueError, match="single A-Za-z"):
            fill_chainname.run(args)

    def test_user_chainname_duplicates_raises(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")

        pymol_cmd = MagicMock()
        pymol_cmd.get_chains.return_value = []
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(
            fill_chainname,
            "_collect_blank_atoms",
            lambda _: [(1, 1, "S1"), (2, 100, "S2")],
        )

        args = self._make_args(structure, tmp_path / "out.pdb", chainname=["B", "B"])
        with pytest.raises(ValueError, match="duplicates"):
            fill_chainname.run(args)

    def test_uses_index_selection(self, tmp_path, monkeypatch):
        from src.build import fill_chainname

        structure = tmp_path / "input.pdb"
        structure.write_text("ATOM\n")
        output = tmp_path / "out.pdb"

        pymol_cmd = MagicMock()
        pymol_cmd.get_chains.return_value = []
        monkeypatch.setattr(fill_chainname, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(
            fill_chainname,
            "_collect_blank_atoms",
            lambda _: [(7, 1, ""), (8, 2, "")],
        )

        fill_chainname.run(self._make_args(structure, output))

        sel_arg = pymol_cmd.alter.call_args_list[0].args[0]
        assert "index 7+8" in sel_arg
