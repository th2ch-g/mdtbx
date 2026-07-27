"""
analysis/extract_ave_str unit tests

Tests average structure extraction with MDtraj (the gmx=False path).
"""

import types

import numpy as np
import pytest


class TestExtractAveStrRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            output=str(output),
            gmx=False,
            index=None,
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """run() writes a PDB file"""
        from mdtbx.analysis.extract_ave_str import run

        out = tmp_path / "ave.pdb"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_has_single_frame(self, trajectory_files, tmp_path):
        """The average structure has exactly one frame"""
        import mdtraj as md

        from mdtbx.analysis.extract_ave_str import run

        out = tmp_path / "ave_single.pdb"
        run(self._make_args(trajectory_files, out))

        loaded = md.load_pdb(str(out))
        assert loaded.n_frames == 1

    def test_output_has_correct_atom_count_all(self, trajectory_files, tmp_path):
        """Selecting all atoms keeps the atom count of the source trajectory"""
        import mdtraj as md

        from mdtbx.analysis.extract_ave_str import run

        out = tmp_path / "ave_all.pdb"
        run(self._make_args(trajectory_files, out, selection="all"))

        loaded = md.load_pdb(str(out))
        assert loaded.n_atoms == trajectory_files["traj"].n_atoms

    def test_output_has_correct_atom_count_subset(self, trajectory_files, tmp_path):
        """A partial selection yields exactly the selected atom count"""
        import mdtraj as md

        from mdtbx.analysis.extract_ave_str import run

        selection = "resid 0"
        out = tmp_path / "ave_subset.pdb"
        run(self._make_args(trajectory_files, out, selection=selection))

        loaded = md.load_pdb(str(out))
        n_selected = len(trajectory_files["traj"].topology.select(selection))
        assert loaded.n_atoms == n_selected

    def test_average_xyz_is_within_trajectory_range(self, trajectory_files, tmp_path):
        """Average coordinates stay within the per-axis min/max of the source"""
        import mdtraj as md

        from mdtbx.analysis.extract_ave_str import run

        out = tmp_path / "ave_range.pdb"
        run(self._make_args(trajectory_files, out))

        ave = md.load_pdb(str(out))
        traj = trajectory_files["traj"]

        assert np.all(ave.xyz >= traj.xyz.min(axis=0) - 1e-6)
        assert np.all(ave.xyz <= traj.xyz.max(axis=0) + 1e-6)

    def test_nested_output_directory_is_created(self, trajectory_files, tmp_path):
        from mdtbx.analysis.extract_ave_str import run

        out = tmp_path / "nested" / "average.pdb"
        run(self._make_args(trajectory_files, out))

        assert out.exists()

    def test_empty_selection_raises(self, trajectory_files, tmp_path):
        from mdtbx.analysis.extract_ave_str import run

        args = self._make_args(
            trajectory_files,
            tmp_path / "average.pdb",
            selection="name XXXX",
        )
        with pytest.raises(ValueError, match="No atoms selected"):
            run(args)

    def test_gmx_supplies_both_selection_prompts(self, tmp_path, monkeypatch):
        from mdtbx.analysis import extract_ave_str

        captured = {}

        def fake_run_cmd(command, **kwargs):
            captured["command"] = command
            captured["input"] = kwargs["input"]

        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(extract_ave_str, "run_cmd", fake_run_cmd)
        args = types.SimpleNamespace(
            topology="topology.tpr",
            trajectory="trajectory.xtc",
            selection="Backbone",
            output=str(tmp_path / "average.pdb"),
            gmx=True,
            index=None,
        )

        extract_ave_str.run(args)

        assert captured["command"][:2] == ["gmx", "covar"]
        assert captured["input"] == "Backbone\nBackbone\n"
        assert not list(tmp_path.glob("tmp*.xvg"))
        assert not list(tmp_path.glob("tmp*.trr"))
