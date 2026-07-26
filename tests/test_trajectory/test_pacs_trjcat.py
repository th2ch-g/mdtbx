from argparse import Namespace

import pytest

from src.trajectory import pacs_trjcat


def test_check_replica_reads_first_cycle_directory(tmp_path):
    trial_dir = tmp_path / "trial"
    (trial_dir / "cycle000" / "replica001").mkdir(parents=True)
    (trial_dir / "cycle000" / "replica002").mkdir(parents=True)

    args = Namespace(trial_dir=str(trial_dir))

    assert pacs_trjcat.check_cycle(args) == 1
    assert pacs_trjcat.check_replica(args) == 2


def test_check_replica_uses_largest_cycle(tmp_path):
    trial_dir = tmp_path / "trial"
    (trial_dir / "cycle000" / "replica001").mkdir(parents=True)
    for replica in range(1, 4):
        (trial_dir / "cycle001" / f"replica{replica:03}").mkdir(parents=True)

    args = Namespace(trial_dir=str(trial_dir))

    assert pacs_trjcat.check_replica(args) == 3


def test_run_uses_trajectory_extension_for_cleanup(tmp_path, monkeypatch):
    trial_dir = tmp_path / "trial"
    replica_dir = trial_dir / "cycle000" / "replica001"
    replica_dir.mkdir(parents=True)
    (replica_dir / "prd.trr").write_text("trajectory")
    (replica_dir / "rmmol_top.tpr").write_text("topology")

    commands = []

    def fake_run_cmd(command, *, log=None, check=True, **kwargs):
        commands.append(command)

    monkeypatch.setattr(pacs_trjcat, "run_cmd", fake_run_cmd)

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

    rendered_commands = [" ".join(map(str, command)) for command in commands]
    assert any("tmp_all.trr" in command for command in rendered_commands)
    assert not any("tmp_all.xtc" in command for command in rendered_commands)
    assert all(isinstance(command, list) for command in commands)


def test_run_uses_actual_replica_directories_for_each_cycle(tmp_path, monkeypatch):
    trial_dir = tmp_path / "trial"
    initial = trial_dir / "cycle000" / "replica001"
    initial.mkdir(parents=True)
    (initial / "prd.xtc").write_text("initial")
    (initial / "rmmol_top.tpr").write_text("topology")
    for replica in range(1, 4):
        replica_dir = trial_dir / "cycle001" / f"replica{replica:03}"
        replica_dir.mkdir(parents=True)
        (replica_dir / "prd.xtc").write_text(f"replica {replica}")

    calls = []

    def fake_run_cmd(command, **kwargs):
        calls.append((command, kwargs))

    monkeypatch.setattr(pacs_trjcat, "run_cmd", fake_run_cmd)
    args = Namespace(
        trial_dir=str(trial_dir),
        ref_structure=None,
        fit_selection="Protein",
        skip=1,
        trjname="prd.xtc",
        keep_selection="System",
        centering_selection="Protein",
        index=None,
        pbc="mol",
        keep_cycle_trj=False,
    )

    pacs_trjcat.run(args)

    cycle_command, kwargs = next(
        (command, kwargs)
        for command, kwargs in calls
        if command[:2] == ["gmx", "trjcat"]
        and str(trial_dir / "cycle001") in " ".join(command)
    )
    for replica in range(1, 4):
        expected = trial_dir / "cycle001" / f"replica{replica:03}" / "prd.xtc"
        assert str(expected) in cycle_command
    assert kwargs["input"] == "c\n" * 3


def test_cluster_pbc_supplies_all_three_selections(tmp_path, monkeypatch):
    trial_dir = tmp_path / "trial"
    replica_dir = trial_dir / "cycle000" / "replica001"
    replica_dir.mkdir(parents=True)
    (replica_dir / "prd.xtc").write_text("trajectory")
    (replica_dir / "rmmol_top.tpr").write_text("topology")
    calls = []
    monkeypatch.setattr(
        pacs_trjcat,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs)),
    )
    args = Namespace(
        trial_dir=str(trial_dir),
        ref_structure=None,
        fit_selection="Backbone",
        skip=1,
        trjname="prd.xtc",
        keep_selection="Protein",
        centering_selection="Protein",
        index=None,
        pbc="cluster",
        keep_cycle_trj=False,
    )

    pacs_trjcat.run(args)

    pbc_calls = [
        kwargs
        for command, kwargs in calls
        if command[:2] == ["gmx", "trjconv"] and "-pbc" in command
    ]
    assert len(pbc_calls) == 2
    assert all(kwargs["input"] == "Protein\nProtein\nSystem\n" for kwargs in pbc_calls)


def test_cycle_directories_are_sorted_numerically_and_noise_is_ignored(tmp_path):
    (tmp_path / "cycle10").mkdir()
    (tmp_path / "cycle2").mkdir()
    (tmp_path / "cycle-backup").mkdir()

    directories = pacs_trjcat._list_subdirs(tmp_path, "cycle")

    assert [path.name for path in directories] == ["cycle2", "cycle10"]


def test_missing_trajectory_is_rejected_before_external_commands(tmp_path, monkeypatch):
    trial_dir = tmp_path / "trial"
    replica_dir = trial_dir / "cycle000" / "replica001"
    replica_dir.mkdir(parents=True)
    (replica_dir / "rmmol_top.tpr").write_text("topology")
    monkeypatch.setattr(
        pacs_trjcat,
        "run_cmd",
        lambda *args, **kwargs: pytest.fail("external command should not run"),
    )
    args = Namespace(
        trial_dir=str(trial_dir),
        ref_structure=None,
        fit_selection="Protein",
        skip=1,
        trjname="missing.xtc",
        keep_selection="System",
        centering_selection="Protein",
        index=None,
        pbc="mol",
        keep_cycle_trj=False,
    )

    with pytest.raises(FileNotFoundError, match="missing.xtc"):
        pacs_trjcat.run(args)
