import json
from types import SimpleNamespace

from mdtbx.build import setup_abfe


def _anchor_pdb(path):
    coordinates = [
        (0.0, 10.0, 0.0),
        (0.0, 0.0, 0.0),
        (10.0, 0.0, 0.0),
        (10.0, 10.0, 0.0),
        (20.0, 10.0, 0.0),
        (20.0, 10.0, 10.0),
    ]
    lines = []
    for index, (x, y, z) in enumerate(coordinates, start=1):
        residue = "REC" if index <= 3 else "LIG"
        lines.append(
            f"ATOM  {index:5d}  C{index:<2} {residue:>3} A{index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
    path.write_text("".join(lines) + "END\n")


def _args(tmp_path):
    mdp = tmp_path / "base.mdp"
    topology = tmp_path / "system.top"
    structure = tmp_path / "system.pdb"
    mdp.write_text(
        "integrator = md\n"
        "nsteps = 10\n"
        "vdwtype = Cut-off\n"
        "pull = yes\n"
        "pull-coord7-k = 1\n"
    )
    topology.write_text("topology\n")
    _anchor_pdb(structure)
    return SimpleNamespace(
        mdp=str(mdp),
        complex_mdp=None,
        solvent_mdp=None,
        complex_topology=str(topology),
        complex_structure=str(structure),
        solvent_topology=str(topology),
        solvent_structure=str(structure),
        moltype="LIG",
        outdir=str(tmp_path / "abfe"),
        anchors=[1, 2, 3, 4, 5, 6],
        receptor_selection="protein",
        ligand_selection=None,
        anchor_search_distance=0.5,
        distance_spring=4184.0,
        angle_spring=41.84,
        dihedral_spring=41.84,
        temperature=300.0,
        lj_pme_comb_rule="Lorentz-Berthelot",
        charge_windows=2,
        vdw_windows=2,
        restraint_windows=2,
        nstdhdl=5,
        complex_checkpoint=None,
        solvent_checkpoint=None,
        complex_reference=None,
        solvent_reference=None,
        complex_index=None,
        solvent_index=None,
        deffnm="abfe",
        gmx="gmx",
        maxwarn=0,
        no_tpr=True,
    )


def test_run_generates_five_leg_abfe_cycle(tmp_path, monkeypatch):
    def fake_run_cmd(command, **_kwargs):
        output = command[command.index("-o") + 1]
        setup_abfe.Path(output).write_text("[ System ]\n1 2 3 4 5 6\n")

    monkeypatch.setattr(setup_abfe, "run_cmd", fake_run_cmd)
    setup_abfe.run(_args(tmp_path))

    base = tmp_path / "abfe"
    manifest = json.loads((base / "abfe_manifest.json").read_text())
    assert set(manifest["legs"]) == {
        "solvent_charge",
        "solvent_vdw",
        "complex_charge",
        "complex_vdw",
        "complex_restraint",
    }
    assert manifest["long_range_method"] == "LJ-PME"
    assert manifest["lj_pme_comb_rule"] == "Lorentz-Berthelot"
    assert manifest["gmx_executable"] == "gmx"
    for leg in manifest["legs"]:
        assert (base / leg / "fep_manifest.json").is_file()

    restrained = (base / "inputs" / "complex_restrained.mdp").read_text()
    release = (base / "inputs" / "complex_restraint_release.mdp").read_text()
    solvent = (base / "inputs" / "solvent.mdp").read_text()
    assert "pull-coord6-kB" in restrained
    assert "pull-coord6-kB" in release
    assert "pull-coord7-k" not in restrained
    assert "pull" not in solvent
    assert "vdwtype" in solvent
    assert "PME" in solvent
    assert "Lorentz-Berthelot" in solvent
    index = (base / "inputs" / "complex_abfe.ndx").read_text()
    assert "[ System ]" in index
    assert "[ ABFE_6 ]" in index


def test_run_records_trajectory_anchor_statistics(tmp_path, monkeypatch):
    import mdtraj as md
    import numpy as np

    args = _args(tmp_path)
    frame = md.load(args.complex_structure)
    trajectory = tmp_path / "complex.xtc"
    md.Trajectory(np.repeat(frame.xyz, 2, axis=0), frame.topology).save_xtc(
        str(trajectory)
    )
    args.anchor_trajectory = str(trajectory)
    args.anchor_stride = 1

    def fake_run_cmd(command, **_kwargs):
        output = command[command.index("-o") + 1]
        setup_abfe.Path(output).write_text("[ System ]\n1 2 3 4 5 6\n")

    monkeypatch.setattr(setup_abfe, "run_cmd", fake_run_cmd)
    setup_abfe.run(args)

    manifest = json.loads((tmp_path / "abfe" / "abfe_manifest.json").read_text())
    selection = manifest["anchor_selection"]
    assert selection["method"] == "explicit_trajectory"
    assert selection["frame_count"] == 2
    assert selection["trajectory"] == str(trajectory.resolve())


def test_run_accepts_per_leg_optimized_schedule(tmp_path, monkeypatch):
    args = _args(tmp_path)
    schedule = tmp_path / "charge_schedule.json"
    schedule.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "workflow": "optimized-fep-schedule",
                "source_workflow": "fep",
                "mode": "charge",
                "state_count": 3,
                "coordinates": [0.0, 0.25, 1.0],
                "lambda_components": {
                    "coul": [0.0, 0.25, 1.0],
                    "vdw": [0.0, 0.0, 0.0],
                },
            }
        )
    )
    args.solvent_charge_schedule = str(schedule)

    def fake_run_cmd(command, **_kwargs):
        output = command[command.index("-o") + 1]
        setup_abfe.Path(output).write_text("[ System ]\n1 2 3 4 5 6\n")

    monkeypatch.setattr(setup_abfe, "run_cmd", fake_run_cmd)
    setup_abfe.run(args)

    base = tmp_path / "abfe"
    leg_manifest = json.loads(
        (base / "solvent_charge" / "fep_manifest.json").read_text()
    )
    manifest = json.loads((base / "abfe_manifest.json").read_text())
    assert len(leg_manifest["windows"]) == 3
    assert manifest["schedule_sources"]["solvent_charge"] == str(schedule.resolve())
