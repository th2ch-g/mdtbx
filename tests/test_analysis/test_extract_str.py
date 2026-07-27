"""
analysis/extract_str unit tests

Tests single-frame structure extraction with MDtraj (the gmx=False path).
"""

import types

import pytest


class TestExtractStrRun:
    def _make_args(self, traj_files, output, time=1, selection="all"):
        return types.SimpleNamespace(
            topology=traj_files["pdb"],
            trajectory=traj_files["xtc"],
            selection=selection,
            time=time,
            output=str(output),
            gmx=False,
            index=None,
        )

    def test_output_file_created(self, trajectory_files, tmp_path):
        """run() writes a PDB file"""
        from src.analysis.extract_str import run

        out = tmp_path / "frame1.pdb"
        run(self._make_args(trajectory_files, out, time=1))
        assert out.exists()

    def test_output_is_valid_pdb(self, trajectory_files, tmp_path):
        """The generated PDB can be read back by MDtraj"""
        import mdtraj as md

        from src.analysis.extract_str import run

        out = tmp_path / "frame1.pdb"
        run(self._make_args(trajectory_files, out, time=1))

        loaded = md.load_pdb(str(out))
        assert loaded.n_frames == 1

    def test_output_atom_count(self, trajectory_files, tmp_path):
        """The extracted structure keeps the atom count of the source trajectory"""
        import mdtraj as md

        from src.analysis.extract_str import run

        out = tmp_path / "frame1.pdb"
        run(self._make_args(trajectory_files, out, time=1))

        loaded = md.load_pdb(str(out))
        assert loaded.n_atoms == trajectory_files["traj"].n_atoms

    def test_extract_different_frames(self, trajectory_files, tmp_path):
        """Different frames give different coordinates"""
        import mdtraj as md

        from src.analysis.extract_str import run

        out1 = tmp_path / "frame1.pdb"
        out2 = tmp_path / "frame2.pdb"
        run(self._make_args(trajectory_files, out1, time=1))
        run(self._make_args(trajectory_files, out2, time=2))

        t1 = md.load_pdb(str(out1))
        t2 = md.load_pdb(str(out2))
        # The coordinates are random, so the two frames must differ
        import numpy as np

        assert not np.allclose(t1.xyz, t2.xyz)

    def test_nested_output_directory_is_created(self, trajectory_files, tmp_path):
        from src.analysis.extract_str import run

        out = tmp_path / "nested" / "frame.pdb"
        run(self._make_args(trajectory_files, out))

        assert out.exists()

    def test_time_is_ps_not_frame_index(self, trajectory_files, tmp_path):
        """--time is interpreted in ps, not as a frame index"""
        import numpy as np

        from src.analysis.extract_str import run

        traj = trajectory_files["traj"]
        # The fixture trajectory carries time = 0, 1, ..., 9 ps
        out = tmp_path / "at_3ps.pdb"
        run(self._make_args(trajectory_files, out, time=3))

        import mdtraj as md

        loaded = md.load_pdb(str(out))
        # ps semantics selects frame 3; a frame index would have given frame 2
        assert np.allclose(loaded.xyz[0], traj.xyz[3], atol=1e-3)
        assert not np.allclose(loaded.xyz[0], traj.xyz[2], atol=1e-3)

    def test_time_zero_is_the_first_frame(self, trajectory_files, tmp_path):
        """--time 0 selects the first frame (the old code rejected it)"""
        import mdtraj as md
        import numpy as np

        from src.analysis.extract_str import run

        out = tmp_path / "at_0ps.pdb"
        run(self._make_args(trajectory_files, out, time=0))

        loaded = md.load_pdb(str(out))
        assert np.allclose(loaded.xyz[0], trajectory_files["traj"].xyz[0], atol=1e-3)

    def test_time_out_of_range_exits(self, trajectory_files, tmp_path):
        """A time outside the trajectory range exits"""
        from src.analysis.extract_str import run

        with pytest.raises(SystemExit):
            run(self._make_args(trajectory_files, tmp_path / "x.pdb", time=999))

    def test_time_between_frames_uses_nearest(self, trajectory_files, tmp_path):
        """A time between frames snaps to the nearest frame"""
        import mdtraj as md
        import numpy as np

        from src.analysis.extract_str import run

        out = tmp_path / "at_2p6ps.pdb"
        run(self._make_args(trajectory_files, out, time=2.6))

        loaded = md.load_pdb(str(out))
        assert np.allclose(loaded.xyz[0], trajectory_files["traj"].xyz[3], atol=1e-3)

    def test_empty_selection_raises(self, trajectory_files, tmp_path):
        from src.analysis.extract_str import run

        args = self._make_args(
            trajectory_files,
            tmp_path / "frame.pdb",
            selection="name XXXX",
        )
        with pytest.raises(ValueError, match="No atoms selected"):
            run(args)

    def test_gmx_uses_argv_and_text_input(self, tmp_path, monkeypatch):
        from src.analysis import extract_str

        captured = {}

        def fake_run_cmd(command, **kwargs):
            captured["command"] = command
            captured["input"] = kwargs["input"]

        monkeypatch.setattr(extract_str, "run_cmd", fake_run_cmd)
        args = types.SimpleNamespace(
            topology="topology with spaces.tpr",
            trajectory="trajectory with spaces.xtc",
            selection="Protein",
            time=10,
            output=str(tmp_path / "frame.pdb"),
            gmx=True,
            index="index with spaces.ndx",
        )

        extract_str.run(args)

        assert captured["command"][:2] == ["gmx", "trjconv"]
        assert "topology with spaces.tpr" in captured["command"]
        assert "index with spaces.ndx" in captured["command"]
        assert captured["input"] == "Protein\n"
