import types

import pytest

from src.build import gaussian


def test_configure_gaussian_input_updates_resources_and_charge(tmp_path):
    input_path = tmp_path / "input.gjf"
    input_path.write_text("%oldchk\n%oldmem\n\nTitle\n\n0 1\nC 0.0 0.0 0.0\n")

    gaussian.configure_gaussian_input(
        input_path,
        checkpoint="job.chk",
        memory_gb=32,
        threads=8,
        route="#p hf/6-31g(d)",
        charge=-1,
        multiplicity=2,
    )

    text = input_path.read_text()
    assert text.startswith("%chk=job.chk\n%mem=32GB\n%nprocshared=8\n#p hf/6-31g(d)\n")
    assert "-1 2\n" in text


def test_configure_gaussian_input_requires_charge_line(tmp_path):
    input_path = tmp_path / "input.gjf"
    input_path.write_text("%chk=old\n%mem=1GB\nC 0.0 0.0 0.0\n")

    with pytest.raises(ValueError, match="Charge and multiplicity"):
        gaussian.configure_gaussian_input(
            input_path,
            checkpoint="job.chk",
            memory_gb=32,
            threads=8,
            route="#p hf/6-31g(d)",
            charge=0,
            multiplicity=1,
        )


def test_convert_to_gaussian_input_uses_argv(monkeypatch, tmp_path):
    captured = {}

    def fake_run_cmd(command, **kwargs):
        captured["command"] = command
        kwargs["stdout"].write("generated\n")

    monkeypatch.setattr(gaussian, "run_cmd", fake_run_cmd)
    output = tmp_path / "output.gjf"

    gaussian.convert_to_gaussian_input(
        "input with spaces.mol2",
        "mol2",
        output,
    )

    assert captured["command"] == [
        "obabel",
        "-i",
        "mol2",
        "input with spaces.mol2",
        "-o",
        "gjf",
    ]
    assert output.read_text() == "generated\n"


def test_run_gaussian_uses_file_streams(monkeypatch, tmp_path):
    captured = {}
    input_path = tmp_path / "input.gjf"
    output_path = tmp_path / "output.log"
    input_path.write_bytes(b"input")

    def fake_run_cmd(command, **kwargs):
        captured["command"] = command
        captured["input"] = kwargs["stdin"].read()
        kwargs["stdout"].write(b"output")

    monkeypatch.setattr(gaussian, "run_cmd", fake_run_cmd)

    gaussian.run_gaussian("g16 -p", input_path, output_path)

    assert captured == {"command": ["g16", "-p"], "input": b"input"}
    assert output_path.read_bytes() == b"output"


def test_gen_resp_no_opt_uses_shared_helpers(monkeypatch, tmp_path):
    from src.build import gen_resp

    conversions = []
    configurations = []
    gaussian_runs = []
    commands = []
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        gen_resp,
        "convert_to_gaussian_input",
        lambda *args: conversions.append(args),
    )
    monkeypatch.setattr(
        gen_resp,
        "configure_gaussian_input",
        lambda *args, **kwargs: configurations.append((args, kwargs)),
    )
    monkeypatch.setattr(
        gen_resp,
        "run_gaussian",
        lambda *args: gaussian_runs.append(args),
    )
    monkeypatch.setattr(
        gen_resp,
        "run_cmd",
        lambda command, **_kwargs: commands.append(command),
    )
    monkeypatch.setattr(gen_resp, "run_tleap", lambda _text: None)
    args = types.SimpleNamespace(
        structure="ligand with spaces.mol2",
        resname="LIG",
        multiplicity=1,
        charge=-1,
        memory=32,
        threads=8,
        no_opt=True,
    )

    gen_resp.run(args)

    assert conversions == [
        (
            "ligand with spaces.mol2",
            "mol2",
            "single_point_calculation.gjf",
        )
    ]
    assert len(configurations) == 1
    assert len(gaussian_runs) == 1
    assert all(isinstance(command, list) for command in commands)


def test_gen_modres_resp_validates_atom_charge_format(monkeypatch, tmp_path):
    from src.build import gen_modres_resp

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        gen_modres_resp,
        "convert_to_gaussian_input",
        lambda *_args: None,
    )
    monkeypatch.setattr(
        gen_modres_resp,
        "configure_gaussian_input",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        gen_modres_resp,
        "run_gaussian",
        lambda *_args: None,
    )
    monkeypatch.setattr(
        gen_modres_resp,
        "run_cmd",
        lambda *_args, **_kwargs: None,
    )
    args = types.SimpleNamespace(
        structure="residue.mol2",
        resname="MOD",
        sepbond1=None,
        sepbond2=None,
        atomcharge=["invalid"],
        multiplicity=1,
        charge=0,
        memory=32,
        threads=8,
    )

    with pytest.raises(ValueError, match="ATOM:CHARGE"):
        gen_modres_resp.run(args)
