"""
Shared pytest fixtures and PyMOL mock setup

NOTE: the mock injection into sys.modules has to happen before modules that
      import PyMOL, so it is done at the top of this file.
"""

import sys
from unittest.mock import MagicMock

# Some commands import PyMOL at module load time, so it is mocked to let the
# test environment work without a PyMOL GUI.
sys.modules["pymol_plugins"] = MagicMock()
sys.modules["pymol"] = MagicMock()
sys.modules["pymol.cmd"] = MagicMock()

import pathlib  # noqa: E402

import numpy as np  # noqa: E402
import pytest  # noqa: E402

FIXTURES_DIR = pathlib.Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def fixtures_dir() -> pathlib.Path:
    return FIXTURES_DIR


@pytest.fixture(scope="session")
def sample_mdp_path() -> pathlib.Path:
    return FIXTURES_DIR / "sample.mdp"


@pytest.fixture(scope="session")
def sample_top_path() -> pathlib.Path:
    return FIXTURES_DIR / "sample.top"


@pytest.fixture(scope="session")
def sample_pdb_path() -> pathlib.Path:
    return FIXTURES_DIR / "sample.pdb"


@pytest.fixture(scope="session")
def trajectory_files(tmp_path_factory):
    """
    Prepare a trajectory synthesised with mdtraj as temporary files.
    A minimal system: ALA + GLY, 2 residues, 9 atoms, 10 frames.
    scope="session" builds it once for the whole test session.
    """
    import mdtraj as md

    tmp = tmp_path_factory.mktemp("traj")

    top = md.Topology()
    chain = top.add_chain()

    res1 = top.add_residue("ALA", chain)
    top.add_atom("N", md.element.nitrogen, res1)
    top.add_atom("CA", md.element.carbon, res1)
    top.add_atom("CB", md.element.carbon, res1)
    top.add_atom("C", md.element.carbon, res1)
    top.add_atom("O", md.element.oxygen, res1)

    res2 = top.add_residue("GLY", chain)
    top.add_atom("N", md.element.nitrogen, res2)
    top.add_atom("CA", md.element.carbon, res2)
    top.add_atom("C", md.element.carbon, res2)
    top.add_atom("O", md.element.oxygen, res2)

    n_frames = 10
    n_atoms = top.n_atoms  # 9

    np.random.seed(42)
    # Random coordinates in the 0-2 nm range (a different position per frame)
    xyz = np.random.rand(n_frames, n_atoms, 3) * 2.0

    traj = md.Trajectory(xyz, top)

    pdb_path = str(tmp / "sample.pdb")
    xtc_path = str(tmp / "sample.xtc")

    traj[0].save_pdb(pdb_path)
    traj.save_xtc(xtc_path)

    return {"pdb": pdb_path, "xtc": xtc_path, "traj": traj}
