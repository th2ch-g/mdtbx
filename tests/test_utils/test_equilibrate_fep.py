from pathlib import Path
from types import SimpleNamespace

from mdtbx.utils import equilibrate_fep


def test_run_executes_stages_in_order(tmp_path, monkeypatch):
    topology = tmp_path / "dual.top"
    structure = tmp_path / "dual.gro"
    mdp1 = tmp_path / "mini.mdp"
    mdp2 = tmp_path / "eq1.mdp"
    topology.write_text("[ system ]\nTest\n")
    structure.write_text("Test\n0\n1 1 1\n")
    mdp1.write_text("integrator = steep\nnsteps = 10\n")
    mdp2.write_text("integrator = md\nnsteps = 10\nref_t = 310\ngen_temp = 310\n")
    calls = []

    def fake_run(command, cwd=None, **_kwargs):
        calls.append((command, cwd))
        if command[1] == "mdrun":
            Path(cwd, "stage.gro").write_text("Test\n0\n1 1 1\n")
            Path(cwd, "stage.cpt").write_text("checkpoint")

    monkeypatch.setattr(equilibrate_fep, "run_cmd", fake_run)
    outdir = tmp_path / "equal"

    equilibrate_fep.run(
        SimpleNamespace(
            topology=str(topology),
            structure=str(structure),
            mdps=[str(mdp1), str(mdp2)],
            outdir=str(outdir),
            reference=None,
            index=None,
            lambda_state=0,
            temperature=300.0,
            gmx="gmx_mpi",
            ntomp=1,
            maxwarn=0,
            gpu_offload=True,
        )
    )

    assert [command[1] for command, _cwd in calls] == [
        "grompp",
        "mdrun",
        "grompp",
        "mdrun",
    ]
    assert "-t" in calls[2][0]
    assert "-nb" not in calls[1][0]
    assert "-pme" not in calls[1][0]
    assert "-nb" in calls[3][0]
    assert "-pme" in calls[3][0]
    assert (
        "ref-t                    = 300"
        in (outdir / "01_eq1" / "stage.mdp").read_text()
    )
    assert (outdir / "equilibration_manifest.json").is_file()
