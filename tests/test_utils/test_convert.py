import types


class FakeCmd:
    def __init__(self):
        self.calls = []

    def reinitialize(self):
        self.calls.append(("reinitialize",))

    def load(self, structure, name):
        self.calls.append(("load", structure, name))

    def save(self, output, selection):
        self.calls.append(("save", output, selection))


def test_pymol_conversion_uses_isolated_target_and_creates_parent(
    tmp_path, monkeypatch
):
    from mdtbx.utils import convert

    fake_cmd = FakeCmd()
    monkeypatch.setattr(convert, "cmd", fake_cmd)
    output = tmp_path / "nested" / "output.pdb"
    args = types.SimpleNamespace(
        structure="input.pdb",
        output=str(output),
        type="pymol",
    )

    convert.run(args)

    assert output.parent.is_dir()
    assert fake_cmd.calls == [
        ("reinitialize",),
        ("load", "input.pdb", "target"),
        ("save", str(output), "target"),
    ]


def test_mdtraj_conversion_creates_parent(trajectory_files, tmp_path):
    import mdtraj as md

    from mdtbx.utils.convert import run

    output = tmp_path / "nested" / "output.pdb"
    args = types.SimpleNamespace(
        structure=trajectory_files["pdb"],
        output=str(output),
        type="mdtraj",
    )

    run(args)

    assert md.load(output).n_atoms == trajectory_files["traj"].n_atoms


def test_conversion_rejects_unknown_backend(tmp_path):
    import pytest

    from mdtbx.utils.convert import run

    args = types.SimpleNamespace(
        structure="input.pdb",
        output=str(tmp_path / "output.pdb"),
        type="unknown",
    )

    with pytest.raises(ValueError, match="Unsupported conversion type"):
        run(args)
