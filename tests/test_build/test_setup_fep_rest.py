import json
from pathlib import Path
from types import SimpleNamespace

from mdtbx.build import setup_fep_rest


DUAL_TOPOLOGY = """\
[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
C 6 12.0 0.0 A 0.34 0.45
DUM 0 12.0 0.0 A 0.0 0.0
[ moleculetype ]
LIG 3
[ atoms ]
1 C 1 LIG C1 1 0.2 12.0 DUM 0.0 12.0
2 C 1 LIG C2 2 -0.2 12.0 C -0.2 12.0
[ system ]
Test
[ molecules ]
LIG 1
"""


def _args(tmp_path, **overrides):
    mdp = tmp_path / "base.mdp"
    topology = tmp_path / "dual.top"
    structure = tmp_path / "system.pdb"
    mdp.write_text("integrator = md\nnsteps = 10\n")
    topology.write_text(DUAL_TOPOLOGY)
    structure.write_text(
        "ATOM      1  C1  LIG A   1       0.000   0.000   0.000"
        "  1.00  0.00           C\n"
        "ATOM      2  C2  LIG A   1       1.500   0.000   0.000"
        "  1.00  0.00           C\nEND\n"
    )
    values = {
        "mdp": str(mdp),
        "topology": str(topology),
        "structure": str(structure),
        "b_structure": None,
        "outdir": str(tmp_path / "fep_rest"),
        "replicas": 4,
        "temperature": 300.0,
        "max_temperature": 900.0,
        "hot_selection": "resname LIG",
        "hot_distance": 0.4,
        "hot_region": "not water",
        "nstdhdl": 5,
        "checkpoint": None,
        "b_checkpoint": None,
        "reference": None,
        "b_reference": None,
        "index": None,
        "deffnm": "rest",
        "gmx": "gmx_mpi",
        "plumed": "plumed",
        "maxwarn": 0,
        "no_tpr": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_run_generates_plumed_rest_replicas(tmp_path, monkeypatch):
    calls = []

    def fake_run_cmd(command, **kwargs):
        calls.append((command, kwargs))
        if command[1] == "grompp":
            Path(command[command.index("-o") + 1]).write_text("tpr")
            if "-pp" in command:
                source = Path(command[command.index("-p") + 1])
                Path(command[command.index("-pp") + 1]).write_text(source.read_text())
        elif command[1] == "partial_tempering":
            kwargs["stdout"].write(kwargs["stdin"].read())

    monkeypatch.setattr(setup_fep_rest, "run_cmd", fake_run_cmd)

    setup_fep_rest.run(_args(tmp_path))

    base = tmp_path / "fep_rest"
    manifest = json.loads((base / "fep_manifest.json").read_text())
    assert manifest["workflow"] == "fep-rest"
    assert len(manifest["windows"]) == 4
    assert manifest["hot_atom_indices"] == [1, 2]
    assert manifest["hrex"] is True
    assert manifest["gmx_executable"] == "gmx_mpi"
    assert manifest["plumed_executable"] == "plumed"
    assert (base / "plumed.dat").is_file()
    assert sum(call[0][1] == "partial_tempering" for call in calls) == 4
    assert sum(call[0][1] == "grompp" for call in calls) == 5

    topology = (base / "lambda_000" / "fep_rest.top").read_text()
    assert "DUM_" in topology
    mdp = (base / "lambda_003" / "rest.mdp").read_text()
    assert "init-lambda-state" in mdp
    assert "= 3" in mdp
    assert "nstxout-compressed" in mdp
    assert "sc-alpha                 = 0.5" in mdp
    assert "sc-coul                  = no" in mdp


def test_run_uses_nearest_endpoint_starting_state(tmp_path, monkeypatch):
    b_structure = tmp_path / "system_b.pdb"
    b_structure.write_text(
        "ATOM      1  C1  LIG A   1       0.100   0.000   0.000"
        "  1.00  0.00           C\n"
        "ATOM      2  C2  LIG A   1       1.600   0.000   0.000"
        "  1.00  0.00           C\nEND\n"
    )
    a_checkpoint = tmp_path / "state_a.cpt"
    b_checkpoint = tmp_path / "state_b.cpt"
    a_checkpoint.write_text("a")
    b_checkpoint.write_text("b")
    calls = []

    def fake_run_cmd(command, **kwargs):
        calls.append((command, kwargs))
        if command[1] == "grompp":
            Path(command[command.index("-o") + 1]).write_text("tpr")
            if "-pp" in command:
                source = Path(command[command.index("-p") + 1])
                Path(command[command.index("-pp") + 1]).write_text(source.read_text())
        elif command[1] == "partial_tempering":
            kwargs["stdout"].write(kwargs["stdin"].read())

    monkeypatch.setattr(setup_fep_rest, "run_cmd", fake_run_cmd)
    args = _args(
        tmp_path,
        b_structure=str(b_structure),
        checkpoint=str(a_checkpoint),
        b_checkpoint=str(b_checkpoint),
    )
    setup_fep_rest.run(args)

    grompp_calls = [call[0] for call in calls if call[0][1] == "grompp"]
    first_window = grompp_calls[1]
    last_window = grompp_calls[-1]
    assert first_window[first_window.index("-c") + 1] == args.structure
    assert first_window[first_window.index("-t") + 1] == str(a_checkpoint.resolve())
    assert last_window[last_window.index("-c") + 1] == str(b_structure.resolve())
    assert last_window[last_window.index("-t") + 1] == str(b_checkpoint.resolve())

    manifest = json.loads((tmp_path / "fep_rest" / "fep_manifest.json").read_text())
    assert manifest["windows"][0]["starting_structure"] == args.structure
    assert manifest["windows"][-1]["starting_structure"] == str(b_structure.resolve())
