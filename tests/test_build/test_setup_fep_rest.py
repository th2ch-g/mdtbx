import json
from pathlib import Path
from types import SimpleNamespace

from src.build import setup_fep_rest


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
        "outdir": str(tmp_path / "fep_rest"),
        "replicas": 4,
        "temperature": 300.0,
        "max_temperature": 900.0,
        "hot_selection": "resname LIG",
        "hot_distance": 0.4,
        "hot_region": "not water",
        "nstdhdl": 5,
        "checkpoint": None,
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
