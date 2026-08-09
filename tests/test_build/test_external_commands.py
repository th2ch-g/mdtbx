import types
from pathlib import Path


def test_make_default_index_uses_argv_and_stdin(monkeypatch):
    from mdtbx.build import add_ndx

    captured = {}

    def fake_run_cmd(command, **kwargs):
        captured["command"] = command
        captured["input"] = kwargs["input"]

    monkeypatch.setattr(add_ndx, "run_cmd", fake_run_cmd)
    args = types.SimpleNamespace(
        gro="structure with spaces.gro",
        output="index with spaces.ndx",
    )

    add_ndx.make_default_index(args)

    assert captured["command"] == [
        "gmx",
        "make_ndx",
        "-f",
        "structure with spaces.gro",
        "-o",
        "index with spaces.ndx",
    ]
    assert captured["input"] == "q\n"


def test_centering_gro_uses_unique_topology_and_argv(tmp_path, monkeypatch):
    from mdtbx.build import centering_gro

    calls = []

    def fake_run_cmd(command, **kwargs):
        calls.append((command, kwargs.get("input")))

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(centering_gro, "run_cmd", fake_run_cmd)
    args = types.SimpleNamespace(
        gmx="gmx_mpi",
        structure="structure with spaces.gro",
        topology="topology with spaces.top",
        index="index with spaces.ndx",
        output="centered structure.gro",
        centering_selection="Protein",
        no_editconf=False,
    )

    centering_gro.run(args)

    assert [command[1] for command, _input in calls] == [
        "grompp",
        "trjconv",
        "editconf",
    ]
    assert calls[1][1] == "Protein\nSystem\n"
    temporary_topology = Path(calls[0][0][calls[0][0].index("-o") + 1])
    centered_path = Path(calls[1][0][calls[1][0].index("-o") + 1])
    assert centered_path != Path(args.output)
    assert calls[2][0][calls[2][0].index("-f") + 1] == str(centered_path)
    assert calls[2][0][calls[2][0].index("-o") + 1] == args.output
    assert not temporary_topology.exists()
    assert not centered_path.exists()


def test_amb2gro_parmed_preserves_paths_as_single_arguments(monkeypatch):
    from mdtbx.build import amb2gro

    calls = []
    monkeypatch.setattr(
        amb2gro,
        "run_cmd",
        lambda command, **_kwargs: calls.append(command),
    )
    args = types.SimpleNamespace(
        type="parmed",
        parm="input files/leap.parm7",
        rst="input files/leap.rst7",
        no_editconf=True,
    )

    amb2gro.run(args)

    assert calls[0][calls[0].index("-p") + 1] == "input files/leap.parm7"
    assert calls[0][calls[0].index("-c") + 1] == "input files/leap.rst7"


def test_amb2gro_acpype_uses_isolated_temporary_directory(tmp_path, monkeypatch):
    from mdtbx.build import amb2gro

    workdirs = []

    def fake_run_cmd(_command, **kwargs):
        workdir = Path(kwargs["cwd"])
        workdirs.append(workdir)
        generated = workdir / "leap.amb2gmx"
        generated.mkdir()
        (generated / "leap_GMX.gro").write_text("gro\n")
        (generated / "leap_GMX.top").write_text("top\n")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(amb2gro, "run_cmd", fake_run_cmd)
    args = types.SimpleNamespace(
        type="acpype",
        parm="leap.parm7",
        rst="leap.rst7",
        no_editconf=True,
    )

    amb2gro.run(args)

    assert (tmp_path / "gmx.gro").read_text() == "gro\n"
    assert (tmp_path / "gmx.top").read_text() == "top\n"
    assert not workdirs[0].exists()


def test_gen_am1bcc_reuses_tleap_helper(monkeypatch):
    from mdtbx.build import gen_am1bcc

    commands = []
    tleap_inputs = []
    monkeypatch.setattr(
        gen_am1bcc,
        "run_cmd",
        lambda command, **_kwargs: commands.append(command),
    )
    monkeypatch.setattr(
        gen_am1bcc,
        "run_tleap",
        lambda input_text: tleap_inputs.append(input_text),
    )
    args = types.SimpleNamespace(
        structure="ligand with spaces.mol",
        resname="LIG",
        charge=-1,
        multiplicity=1,
    )

    gen_am1bcc.run(args)

    assert all(isinstance(command, list) for command in commands)
    assert commands[0][commands[0].index("-fi") + 1] == "mdl"
    assert commands[0][commands[0].index("-i") + 1] == "ligand with spaces.mol"
    assert len(tleap_inputs) == 1
    assert "saveoff LIG LIG.lib" in tleap_inputs[0]
