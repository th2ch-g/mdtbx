"""
cv/rmsf unit tests

Tests the RMSF computation with MDtraj (the gmx=False path).
"""

import types
from pathlib import Path

import numpy as np


class TestRmsfRun:
    def _make_args(self, traj_files, output, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            output=str(output),
            gmx=False,
            resolution="atom",
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy file is generated"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape_all_atoms(self, trajectory_files, tmp_path):
        """Selecting all atoms gives an output length equal to the atom count"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_all.npy"
        run(self._make_args(trajectory_files, out))
        rmsf = np.load(str(out))

        n_atoms = trajectory_files["traj"].n_atoms
        assert rmsf.shape == (n_atoms,)

    def test_output_shape_subset(self, trajectory_files, tmp_path):
        """A partial selection gives an output length equal to the selected count"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_subset.npy"
        run(self._make_args(trajectory_files, out, selection="resid 0"))
        rmsf = np.load(str(out))

        n_selected = len(trajectory_files["traj"].topology.select("resid 0"))
        assert rmsf.shape == (n_selected,)

    def test_output_nonnegative(self, trajectory_files, tmp_path):
        """RMSF is always non-negative"""
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_nn.npy"
        run(self._make_args(trajectory_files, out))
        rmsf = np.load(str(out))
        assert np.all(rmsf >= 0)

    def test_static_trajectory_gives_zero_rmsf(self, tmp_path_factory):
        """RMSF is 0 for a trajectory whose frames share identical coordinates"""
        import mdtraj as md

        from src.cv.rmsf import run

        tmp = tmp_path_factory.mktemp("static")

        # Build a static trajectory
        top = md.Topology()
        chain = top.add_chain()
        res = top.add_residue("ALA", chain)
        top.add_atom("CA", md.element.carbon, res)

        xyz = np.zeros((5, 1, 3))
        traj = md.Trajectory(xyz, top)

        pdb_path = str(tmp / "static.pdb")
        xtc_path = str(tmp / "static.xtc")
        traj[0].save_pdb(pdb_path)
        traj.save_xtc(xtc_path)

        out = tmp / "rmsf_static.npy"
        run(
            types.SimpleNamespace(
                topology=pdb_path,
                trajectory=xtc_path,
                selection="all",
                output=str(out),
                gmx=False,
                resolution="atom",
            )
        )
        rmsf = np.load(str(out))
        assert np.allclose(rmsf, 0.0, atol=1e-6)

    def test_residue_resolution_returns_one_value_per_selected_residue(
        self, trajectory_files, tmp_path
    ):
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_residue.npy"
        args = self._make_args(trajectory_files, out)
        args.resolution = "residue"

        run(args)

        rmsf = np.load(out)
        selected_residues = {
            atom.residue.index for atom in trajectory_files["traj"].topology.atoms
        }
        assert rmsf.shape == (len(selected_residues),)

    def test_gromacs_style_protein_selection_works_with_mdtraj(
        self, trajectory_files, tmp_path
    ):
        from src.cv.rmsf import run

        out = tmp_path / "rmsf_protein.npy"
        args = self._make_args(trajectory_files, out, selection="Protein")

        run(args)

        assert out.exists()

    def test_empty_selection_raises(self, trajectory_files, tmp_path):
        import pytest

        from src.cv.rmsf import run

        args = self._make_args(
            trajectory_files, tmp_path / "rmsf_empty.npy", selection="name XXXX"
        )
        with pytest.raises(ValueError, match="No atoms selected"):
            run(args)

    def test_gmx_single_value_output(self, tmp_path, monkeypatch):
        from src.cv import rmsf

        monkeypatch.chdir(tmp_path)

        def fake_run_cmd(command, **_kwargs):
            output = command[command.index("-o") + 1]
            Path(output).write_text("1 0.25\n")

        monkeypatch.setattr(rmsf, "run_cmd", fake_run_cmd)
        out = tmp_path / "rmsf_gmx.npy"
        args = types.SimpleNamespace(
            topology="topol.tpr",
            trajectory="traj.xtc",
            selection="Protein",
            output=str(out),
            gmx=True,
            resolution="atom",
        )

        rmsf.run(args)

        assert np.load(out).tolist() == [0.25]
