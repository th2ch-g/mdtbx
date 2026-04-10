"""
PyMOL plugin PCA visualization tests.
"""

import importlib
import sys
from pathlib import Path

import numpy as np

PLUGIN_ROOT = Path(__file__).resolve().parents[2] / "pymol-plugins"
if str(PLUGIN_ROOT) not in sys.path:
    sys.path.insert(0, str(PLUGIN_ROOT))


class FakeCmd:
    def __init__(self):
        self.group_calls = []

    def get_unused_name(self, prefix):
        return f"{prefix}_001"

    def group(self, group_name, member_name):
        self.group_calls.append((group_name, member_name))

    def extend(self, *_args, **_kwargs):
        return None


def test_show_pca_mode_from_npz(monkeypatch, tmp_path):
    sys.modules.pop("pymol_plugins", None)
    sys.modules.pop("pymol_plugins.pca", None)
    pymol_pca = importlib.import_module("pymol_plugins.pca")

    npz_path = tmp_path / "pca_backbone.npz"
    np.savez(
        npz_path,
        components=np.array(
            [
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ]
        ),
        explained_variance=np.array([0.25]),
        explained_variance_ratio=np.array([0.6]),
        mean_xyz=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 2.0, 3.0],
                [0.0, 0.0, 0.0],
                [4.0, 5.0, 6.0],
            ]
        ),
        atom_names=np.array(["N", "CA", "C", "CA"]),
        residue_ids=np.array([1, 1, 1, 2]),
        residue_names=np.array(["ALA", "ALA", "ALA", "GLY"]),
        unit_scale_to_angstrom=10.0,
    )

    fake_cmd = FakeCmd()
    arrow_calls = []

    def fake_cgo_arrow(start, end, radius, color, name):
        arrow_calls.append(
            {
                "start": start,
                "end": end,
                "radius": radius,
                "color": color,
                "name": name,
            }
        )

    monkeypatch.setattr(pymol_pca, "cmd", fake_cmd)
    monkeypatch.setattr(pymol_pca, "cgo_arrow", fake_cgo_arrow)

    pymol_pca.show_pca_mode_from_npz(
        str(npz_path),
        pc="1",
        scale="2.0",
        stride="1",
        atom_name="CA",
        radius="0.2",
        color="yellow red",
        name_prefix="test_mode",
    )

    assert len(arrow_calls) == 2
    assert arrow_calls[0]["start"] == [10.0, 20.0, 30.0]
    assert arrow_calls[0]["end"] == [20.0, 20.0, 30.0]
    assert arrow_calls[1]["start"] == [40.0, 50.0, 60.0]
    assert arrow_calls[1]["end"] == [40.0, 60.0, 60.0]
    assert fake_cmd.group_calls == [
        ("test_mode_pc1_001", "test_mode_pc1_001_1"),
        ("test_mode_pc1_001", "test_mode_pc1_001_2"),
    ]
