from argparse import Namespace

from src.trajectory import pacs_trjcat


def test_check_replica_reads_first_cycle_directory(tmp_path):
    trial_dir = tmp_path / "trial"
    (trial_dir / "cycle000" / "replica001").mkdir(parents=True)
    (trial_dir / "cycle000" / "replica002").mkdir(parents=True)

    args = Namespace(trial_dir=str(trial_dir))

    assert pacs_trjcat.check_cycle(args) == 1
    assert pacs_trjcat.check_replica(args) == 2


def test_run_uses_trajectory_extension_for_cleanup(tmp_path, monkeypatch):
    trial_dir = tmp_path / "trial"
    (trial_dir / "cycle000" / "replica001").mkdir(parents=True)

    commands = []

    def fake_run(command, shell, check):
        commands.append(command)

    monkeypatch.setattr(pacs_trjcat.subprocess, "run", fake_run)

    args = Namespace(
        trial_dir=str(trial_dir),
        ref_structure=pacs_trjcat.DEFAULT_TOPOLOGY,
        fit_selection="Protein",
        skip=1,
        trjname="prd.trr",
        keep_selection="System",
        centering_selection="Protein",
        index=None,
        pbc="mol",
        keep_cycle_trj=False,
    )

    pacs_trjcat.run(args)

    assert any("tmp_all.trr" in command for command in commands)
    assert not any("tmp_all.xtc" in command for command in commands)
