"""
pytest 共有 fixture と PyMOL モック設定

NOTE: sys.modules へのモック注入は他のどの src インポートよりも前に
      実行される必要があるため、このファイルの先頭で行う。
"""

import sys
from unittest.mock import MagicMock

# config.py が `import pymol_plugins` を実行し、
# find_bond.py / convert.py が `from pymol import cmd` を実行するため、
# テスト環境では PyMOL GUI なしで動作できるようモックする。
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
    mdtraj で合成した軌跡を一時ファイルとして用意する。
    ALA + GLY の 2 残基 9 原子、10 フレームの最小系。
    scope="session" で全テストセッション中に 1 回だけ生成する。
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
    # 0〜2 nm の範囲でランダムな座標（各フレームで異なる位置）
    xyz = np.random.rand(n_frames, n_atoms, 3) * 2.0

    traj = md.Trajectory(xyz, top)

    pdb_path = str(tmp / "sample.pdb")
    xtc_path = str(tmp / "sample.xtc")

    traj[0].save_pdb(pdb_path)
    traj.save_xtc(xtc_path)

    return {"pdb": pdb_path, "xtc": xtc_path, "traj": traj}
