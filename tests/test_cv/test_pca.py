"""
cv/pca unit tests

Tests PCA with MDtraj + scikit-learn (the gmx=False path).
"""

import types
import shutil

import numpy as np
import pytest


def _mock_gmx_outputs(
    monkeypatch, eigenvector_source, eigenvalue_source, projection_source
):
    def fake_run(command, **_kwargs):
        if command[1] == "covar":
            shutil.copy2(eigenvalue_source, command[command.index("-o") + 1])
            shutil.copy2(eigenvector_source, command[command.index("-v") + 1])
        elif command[1] == "anaeig":
            shutil.copy2(projection_source, command[command.index("-proj") + 1])

    monkeypatch.setattr("subprocess.run", fake_run)


class TestPcaRun:
    def _make_args(
        self,
        traj_files,
        output,
        n_components=2,
        output_npz=None,
        output_average=None,
        selection_cal_trj="all",
        selection_cal_ref="all",
        selection_fit_trj="all",
        selection_fit_ref="all",
    ):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            reference=traj_files["pdb"],
            selection_cal_trj=selection_cal_trj,
            selection_cal_ref=selection_cal_ref,
            selection_fit_trj=selection_fit_trj,
            selection_fit_ref=selection_fit_ref,
            output=str(output),
            gmx=False,
            n_components=n_components,
            index=None,
            output_npz=None if output_npz is None else str(output_npz),
            output_average=None if output_average is None else str(output_average),
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """.npy file is generated"""
        from src.cv.pca import run

        out = tmp_path / "pca.npy"
        run(self._make_args(trajectory_files, out))
        assert out.exists()

    def test_output_shape(self, trajectory_files, tmp_path):
        """The output shape is (n_frames, n_components)"""
        from src.cv.pca import run

        n_components = 3
        out = tmp_path / "pca_shape.npy"
        run(self._make_args(trajectory_files, out, n_components=n_components))
        pc = np.load(str(out))

        n_frames = trajectory_files["traj"].n_frames
        assert pc.shape == (n_frames, n_components)

    def test_first_pc_has_max_variance(self, trajectory_files, tmp_path):
        """The variance of PC1 is at least that of PC2 (by definition of PCA)"""
        from src.cv.pca import run

        out = tmp_path / "pca_var.npy"
        run(self._make_args(trajectory_files, out, n_components=2))
        pc = np.load(str(out))

        var1 = np.var(pc[:, 0])
        var2 = np.var(pc[:, 1])
        assert var1 >= var2

    def test_subset_selection(self, trajectory_files, tmp_path):
        """A partial selection still computes"""
        from src.cv.pca import run

        out = tmp_path / "pca_subset.npy"
        args = self._make_args(
            trajectory_files,
            out,
            n_components=2,
            selection_cal_trj="resid 0",
            selection_cal_ref="resid 0",
            selection_fit_trj="resid 0",
            selection_fit_ref="resid 0",
        )
        run(args)
        pc = np.load(str(out))
        assert pc.shape[0] == trajectory_files["traj"].n_frames

    def test_output_npz_contains_metadata(self, trajectory_files, tmp_path):
        """The .npz for PyMOL visualisation is generated"""
        from src.cv.pca import run

        out = tmp_path / "pca.npy"
        out_npz = tmp_path / "pca.npz"
        run(
            self._make_args(
                trajectory_files,
                out,
                n_components=2,
                output_npz=out_npz,
            )
        )

        data = np.load(str(out_npz))
        pc = np.load(str(out))

        assert out_npz.exists()
        assert data["scores"].shape == pc.shape
        assert np.allclose(data["scores"], pc)
        assert data["components"].shape == (2, trajectory_files["traj"].n_atoms, 3)
        assert data["mean_xyz"].shape == (trajectory_files["traj"].n_atoms, 3)
        assert data["atom_indices"].shape == (trajectory_files["traj"].n_atoms,)
        assert data["coordinate_unit"].item() == "nm"
        assert data["unit_scale_to_angstrom"].item() == 10.0

    def test_output_average_is_created(self, trajectory_files, tmp_path):
        """The average structure file is generated"""
        from src.cv.pca import run

        out = tmp_path / "pca.npy"
        out_average = tmp_path / "pca_average.pdb"
        run(
            self._make_args(
                trajectory_files,
                out,
                output_average=out_average,
            )
        )

        assert out_average.exists()

    def test_nested_output_directories_are_created(self, trajectory_files, tmp_path):
        from src.cv.pca import run

        output_dir = tmp_path / "nested"
        out = output_dir / "pca"
        out_npz = output_dir / "metadata"
        out_average = output_dir / "average.pdb"

        run(
            self._make_args(
                trajectory_files,
                out,
                output_npz=out_npz,
                output_average=out_average,
            )
        )

        assert (output_dir / "pca.npy").exists()
        assert (output_dir / "metadata.npz").exists()
        assert out_average.exists()

    def test_output_npz_with_gmx(self, trajectory_files, tmp_path, monkeypatch):
        """The gmx path also generates the .npz for PyMOL visualisation"""
        import mdtraj as md

        from src.cv.pca import run

        monkeypatch.chdir(tmp_path)

        avg_path = tmp_path / "average.pdb"
        eigenvec_path = tmp_path / "eigenvec.trr"
        eigenval_path = tmp_path / "eigenval.xvg"
        proj_path = tmp_path / "proj.xvg"
        out = tmp_path / "pca_gmx.npy"
        out_npz = tmp_path / "pca_gmx.npz"

        traj = trajectory_files["traj"]
        average_xyz = traj.xyz.mean(axis=0)
        md.Trajectory(average_xyz[None, :, :], traj.topology).save_pdb(str(avg_path))

        eigenvec_xyz = np.stack(
            [
                traj.xyz[0],
                average_xyz,
                np.zeros_like(average_xyz),
                np.zeros_like(average_xyz),
            ]
        )
        eigenvec_xyz[2, 1, 0] = 0.1
        eigenvec_xyz[3, 6, 1] = 0.2
        eigenvec_trj = md.Trajectory(eigenvec_xyz, traj.topology)
        eigenvec_trj.time = np.array([-1.0, 0.0, 4.0, 1.0])
        eigenvec_trj.save_trr(str(eigenvec_path))

        eigenval_path.write_text("1 4.0\n2 1.0\n")
        proj_path.write_text("0.1 0.2\n0.3 0.4\n")

        _mock_gmx_outputs(
            monkeypatch,
            eigenvec_path,
            eigenval_path,
            proj_path,
        )

        args = types.SimpleNamespace(
            topology=trajectory_files["pdb"],
            trajectory=trajectory_files["xtc"],
            reference=trajectory_files["pdb"],
            selection_cal_trj="Backbone",
            selection_cal_ref="Backbone",
            selection_fit_trj="Backbone",
            selection_fit_ref="Backbone",
            output=str(out),
            gmx=True,
            n_components=2,
            index=None,
            output_npz=str(out_npz),
            output_average=str(avg_path),
        )
        run(args)

        data = np.load(str(out_npz))
        saved_pc = np.load(str(out))

        assert np.allclose(saved_pc, np.loadtxt(str(proj_path)))
        assert data["components"].shape == (2, traj.n_atoms, 3)
        assert data["mean_xyz"].shape == (traj.n_atoms, 3)
        assert data["atom_names"].shape == (traj.n_atoms,)
        assert data["explained_variance"].tolist() == pytest.approx([4.0, 1.0])
        assert data["explained_variance_ratio"].tolist() == pytest.approx([0.8, 0.2])
        assert not list(tmp_path.glob("mdtbx-pca-*"))

    def test_output_npz_with_gmx_single_component(
        self, trajectory_files, tmp_path, monkeypatch
    ):
        """The gmx path reads the eigenvalue column correctly for n_components=1"""
        import mdtraj as md

        from src.cv.pca import run

        monkeypatch.chdir(tmp_path)

        avg_path = tmp_path / "average_single.pdb"
        eigenvec_path = tmp_path / "eigenvec.trr"
        eigenval_path = tmp_path / "eigenval.xvg"
        proj_path = tmp_path / "proj.xvg"
        out = tmp_path / "pca_gmx_single.npy"
        out_npz = tmp_path / "pca_gmx_single.npz"

        traj = trajectory_files["traj"]
        average_xyz = traj.xyz.mean(axis=0)
        md.Trajectory(average_xyz[None, :, :], traj.topology).save_pdb(str(avg_path))

        eigenvec_xyz = np.stack(
            [
                traj.xyz[0],
                average_xyz,
                np.zeros_like(average_xyz),
            ]
        )
        eigenvec_xyz[2, 1, 0] = 0.1
        eigenvec_trj = md.Trajectory(eigenvec_xyz, traj.topology)
        eigenvec_trj.time = np.array([-1.0, 0.0, 4.0])
        eigenvec_trj.save_trr(str(eigenvec_path))

        eigenval_path.write_text("1 4.0\n")
        proj_path.write_text("0.1\n0.3\n")

        _mock_gmx_outputs(
            monkeypatch,
            eigenvec_path,
            eigenval_path,
            proj_path,
        )

        args = types.SimpleNamespace(
            topology=trajectory_files["pdb"],
            trajectory=trajectory_files["xtc"],
            reference=trajectory_files["pdb"],
            selection_cal_trj="Backbone",
            selection_cal_ref="Backbone",
            selection_fit_trj="Backbone",
            selection_fit_ref="Backbone",
            output=str(out),
            gmx=True,
            n_components=1,
            index=None,
            output_npz=str(out_npz),
            output_average=str(avg_path),
        )
        run(args)

        data = np.load(str(out_npz))

        assert data["explained_variance"].tolist() == pytest.approx([4.0])
        assert data["explained_variance_ratio"].tolist() == pytest.approx([1.0])

    def test_rejects_non_positive_component_count(self, trajectory_files, tmp_path):
        from src.cv.pca import run

        args = self._make_args(trajectory_files, tmp_path / "pca.npy", n_components=0)

        with pytest.raises(ValueError, match="n_components"):
            run(args)

    def test_rejects_too_many_components(self, trajectory_files, tmp_path):
        from src.cv.pca import run

        n_frames = trajectory_files["traj"].n_frames
        args = self._make_args(
            trajectory_files,
            tmp_path / "pca.npy",
            n_components=n_frames + 1,
        )

        with pytest.raises(ValueError, match="must not exceed"):
            run(args)


def test_load_gmx_projections_removes_time_column(tmp_path):
    from src.cv.pca import _load_gmx_projections

    projections = tmp_path / "proj.xvg"
    projections.write_text("0.0 1.0 2.0\n1.0 3.0 4.0\n")

    values = _load_gmx_projections(projections, n_components=2)

    assert values.tolist() == [[1.0, 2.0], [3.0, 4.0]]
