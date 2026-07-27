import json
import math
from types import SimpleNamespace

from mdtbx.analysis import analyze_abfe


def _manifest(tmp_path):
    legs = {
        "solvent_charge": "solvent_charge",
        "solvent_vdw": "solvent_vdw",
        "complex_charge": "complex_charge",
        "complex_vdw": "complex_vdw",
        "complex_restraint": "complex_restraint",
    }
    manifest = {
        "schema_version": 1,
        "temperature": 300.0,
        "geometry": {
            "distance_nm": 0.4,
            "angles_rad": [math.pi / 2, math.pi / 2],
            "dihedrals_rad": [0.0, 0.0, 0.0],
        },
        "springs": {
            "distance": 4184.0,
            "angle": 41.84,
            "dihedral": 41.84,
        },
        "long_range_method": "LJ-PME",
        "legs": legs,
    }
    (tmp_path / "abfe_manifest.json").write_text(json.dumps(manifest))
    return manifest


def test_run_combines_signed_legs_standard_state_and_correction(
    tmp_path,
    monkeypatch,
):
    _manifest(tmp_path)
    values = {
        "solvent_charge": 1.0,
        "solvent_vdw": 2.0,
        "complex_charge": 3.0,
        "complex_vdw": 4.0,
        "complex_restraint": 5.0,
    }

    def fake_analyze(_args, _base, relative_path, _temperature):
        return {
            "delta_g_kj_mol": values[relative_path],
            "uncertainty_kj_mol": 0.1,
        }

    monkeypatch.setattr(analyze_abfe, "_analyze_leg", fake_analyze)
    args = SimpleNamespace(
        path=str(tmp_path),
        begin=0.0,
        end=None,
        temperature=None,
        symmetry_number=1,
        correction=[["charge_correction", "1.5", "0.2"]],
        output=None,
        gmx="gmx",
    )

    analyze_abfe.run(args)

    result = json.loads((tmp_path / "abfe_result.json").read_text())
    statistical_cycle = 1.0 + 2.0 - 3.0 - 4.0 + 5.0
    standard = next(
        item["value_kj_mol"]
        for item in result["contributions"]
        if item["name"] == "boresch_standard_state"
    )
    assert result["delta_g_binding_kj_mol"] == statistical_cycle + standard + 1.5
    assert result["uncertainty_kj_mol"] == math.sqrt(5 * 0.1**2 + 0.2**2)
    assert (tmp_path / "abfe_result.txt").is_file()
