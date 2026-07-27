import json
from types import SimpleNamespace

from mdtbx.utils import run_abfe


def _setup(tmp_path):
    legs = {
        "solvent_charge": "solvent_charge",
        "solvent_vdw": "solvent_vdw",
        "complex_charge": "complex_charge",
        "complex_vdw": "complex_vdw",
        "complex_restraint": "complex_restraint",
    }
    manifest = {
        "schema_version": 1,
        "workflow": "abfe",
        "legs": legs,
    }
    (tmp_path / "abfe_manifest.json").write_text(json.dumps(manifest))
    return SimpleNamespace(
        path=str(tmp_path),
        legs=["solvent_charge", "complex_restraint"],
        start_window=1,
        end_window=2,
        multidir=True,
        replex=100,
        ntmpi=4,
        ntomp=2,
        maxh=1.0,
        nsteps=10,
        no_continue=True,
        gmx="gmx",
        mpi_launcher="mpirun -np {ntmpi}",
    )


def test_run_forwards_options_to_selected_abfe_legs(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(run_abfe.run_fep, "run", calls.append)

    run_abfe.run(_setup(tmp_path))

    assert [call.path for call in calls] == [
        str(tmp_path / "solvent_charge"),
        str(tmp_path / "complex_restraint"),
    ]
    assert all(call.start_window == 1 for call in calls)
    assert all(call.mpi_launcher == "mpirun -np {ntmpi}" for call in calls)
    assert all(call.no_continue is True for call in calls)
