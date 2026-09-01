import types

import pytest

from mdtbx.build import modeling_cf


def _args(tmp_path, **overrides):
    values = {
        "input": None,
        "sequence": "ACDE",
        "output": str(tmp_path / "results"),
    }
    values.update(overrides)
    return types.SimpleNamespace(**values)


def test_run_uses_argv_and_preserves_unrelated_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(modeling_cf.shutil, "which", lambda _name: "/bin/tool")
    existing_fasta = tmp_path / "input.fasta"
    existing_fasta.write_text("keep\n")
    existing_template = tmp_path / "tmp_template"
    existing_template.mkdir()
    (existing_template / "keep.pdb").write_text("keep\n")
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        output = modeling_cf.Path(command[-2])
        output.mkdir()
        (output / "input_unrelaxed_rank_001_model_1.pdb").write_text("ATOM\n")

    monkeypatch.setattr(modeling_cf, "run_cmd", fake_run)

    modeling_cf.run(_args(tmp_path))

    assert isinstance(captured["command"], list)
    assert captured["command"][0] == "colabfold_batch"
    assert captured["command"][-2:] == [str(tmp_path / "results"), "--amber"]
    assert existing_fasta.read_text() == "keep\n"
    assert (existing_template / "keep.pdb").read_text() == "keep\n"
    assert not list(tmp_path.glob(".mdtbx_colabfold_*"))


def test_run_adds_template_arguments(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(modeling_cf.shutil, "which", lambda _name: "/bin/tool")
    template = tmp_path / "reference.pdb"
    template.write_text("ATOM\n")
    captured = {}

    def fake_run(command, **_kwargs):
        captured["command"] = command
        template_arg = command[command.index("--custom-template-path") + 1]
        copied_template = modeling_cf.Path(template_arg) / "tmpl.pdb"
        assert copied_template.is_file()
        assert len(copied_template.stem) == 4
        output = modeling_cf.Path(command[-2])
        output.mkdir()
        (output / "input_unrelaxed_rank_001_model_1.pdb").write_text("ATOM\n")

    monkeypatch.setattr(modeling_cf, "run_cmd", fake_run)

    modeling_cf.run(_args(tmp_path, input=str(template)))

    assert "--templates" in captured["command"]


def test_missing_command_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(modeling_cf.shutil, "which", lambda _name: None)

    with pytest.raises(RuntimeError, match="not installed"):
        modeling_cf.run(_args(tmp_path))


def test_empty_sequence_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(modeling_cf.shutil, "which", lambda _name: "/bin/tool")

    with pytest.raises(ValueError, match="must not be empty"):
        modeling_cf.run(_args(tmp_path, sequence=" \n"))


def test_missing_ranked_model_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(modeling_cf.shutil, "which", lambda _name: "/bin/tool")
    monkeypatch.setattr(modeling_cf, "run_cmd", lambda *_args, **_kwargs: None)

    with pytest.raises(RuntimeError, match="rank-001 model"):
        modeling_cf.run(_args(tmp_path))
