"""
trajectory/fit unit tests

Tests trajectory fitting with MDtraj (the gmx=False path).
"""

import types
from pathlib import Path

import pytest


class TestFitRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            file=traj_files["xtc"],
            topology=traj_files["pdb"],
            output=str(output),
            selection=selection,
            gmx=False,
            pbc="mol",
            index="index.ndx",
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """run() generates a fitted trajectory file"""

        from mdtbx.trajectory.fit import run

        out = tmp_path / "fitted.xtc"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_has_same_frames(self, trajectory_files, tmp_path):
        """The fitted trajectory keeps the original frame count"""
        import mdtraj as md

        from mdtbx.trajectory.fit import run

        out = tmp_path / "fitted.xtc"
        run(self._make_args(trajectory_files, out))

        fitted = md.load(str(out), top=trajectory_files["pdb"])
        original_n_frames = trajectory_files["traj"].n_frames
        assert fitted.n_frames == original_n_frames

    def test_output_has_same_atoms(self, trajectory_files, tmp_path):
        """The fitted trajectory keeps the original atom count"""
        import mdtraj as md

        from mdtbx.trajectory.fit import run

        out = tmp_path / "fitted.xtc"
        run(self._make_args(trajectory_files, out))

        fitted = md.load(str(out), top=trajectory_files["pdb"])
        original_n_atoms = trajectory_files["traj"].n_atoms
        assert fitted.n_atoms == original_n_atoms

    def test_empty_selection_raises(self, trajectory_files, tmp_path):
        from mdtbx.trajectory.fit import run

        args = self._make_args(
            trajectory_files, tmp_path / "fitted.xtc", selection="name XXXX"
        )

        with pytest.raises(ValueError, match="No atoms selected"):
            run(args)

    def test_nested_output_directory_is_created(self, trajectory_files, tmp_path):
        from mdtbx.trajectory.fit import run

        out = tmp_path / "nested" / "fitted.xtc"
        run(self._make_args(trajectory_files, out))

        assert out.exists()

    def test_gmx_uses_unique_temporary_path(self, tmp_path, monkeypatch):
        from mdtbx.trajectory import fit

        calls = []

        def fake_run_cmd(command, **kwargs):
            calls.append((command, kwargs["input"]))

        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(fit, "run_cmd", fake_run_cmd)
        args = types.SimpleNamespace(
            file="trajectory with spaces.xtc",
            topology="topology with spaces.tpr",
            output=str(tmp_path / "fitted.xtc"),
            selection="Protein",
            gmx=True,
            pbc="mol",
            index="index with spaces.ndx",
        )

        fit.run(args)

        assert len(calls) == 2
        assert all(isinstance(command, list) for command, _input in calls)
        assert all(input_text == "Protein\nSystem\n" for _command, input_text in calls)
        centered_path = Path(calls[0][0][calls[0][0].index("-o") + 1])
        assert not centered_path.exists()
