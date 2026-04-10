"""
build/mutate のユニットテスト

PyMOL Mutagenesis Wizard の呼び出しをモックして検証する。
"""

import types
from unittest.mock import MagicMock

import pytest


class TestMutateRun:
    def _make_args(self, structure, output, selection="chain A and resi 42"):
        return types.SimpleNamespace(
            structure=str(structure),
            selection=selection,
            mutant="VAL",
            output=str(output),
            output_prefix=None,
            rotamer=1,
            hydrogens="auto",
        )

    def test_run_invokes_mutagenesis_wizard(self, tmp_path, monkeypatch):
        from src.build import mutate

        structure = tmp_path / "input.pdb"
        structure.write_text("MODEL\nENDMDL\n")

        wizard = MagicMock()
        pymol_cmd = MagicMock()
        pymol_cmd.get_wizard.return_value = wizard
        monkeypatch.setattr(mutate, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(mutate, "_count_selected_residues", lambda _: 1)

        mutate.run(self._make_args(structure, tmp_path / "mutated"))

        pymol_cmd.wizard.assert_called_once_with("mutagenesis")
        wizard.set_mode.assert_called_once_with("VAL")
        wizard.do_select.assert_called_once_with("_mutate_input")
        wizard.apply.assert_called_once()
        pymol_cmd.save.assert_called_once()

    def test_run_rejects_multi_residue_selection(self, tmp_path, monkeypatch):
        from src.build import mutate

        structure = tmp_path / "input.pdb"
        structure.write_text("MODEL\nENDMDL\n")

        pymol_cmd = MagicMock()
        monkeypatch.setattr(mutate, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(mutate, "_count_selected_residues", lambda _: 2)

        with pytest.raises(ValueError, match="exactly one residue"):
            mutate.run(self._make_args(structure, tmp_path / "mutated"))

    def test_output_prefix_alias_is_supported(self, tmp_path, monkeypatch):
        from src.build import mutate

        structure = tmp_path / "input.pdb"
        structure.write_text("MODEL\nENDMDL\n")

        wizard = MagicMock()
        pymol_cmd = MagicMock()
        pymol_cmd.get_wizard.return_value = wizard
        monkeypatch.setattr(mutate, "pymol_cmd", pymol_cmd)
        monkeypatch.setattr(mutate, "_count_selected_residues", lambda _: 1)

        args = self._make_args(structure, tmp_path / "mutated.pdb")
        args.output_prefix = str(tmp_path / "prefix")
        mutate.run(args)

        saved_path = pymol_cmd.save.call_args.args[0]
        assert saved_path.endswith("prefix.pdb")
