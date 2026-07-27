from pathlib import Path

import pytest

from mdtbx.analysis import analyze_fep_rest


def _energy(path, values):
    path.write_text(
        '@ title "Potential"\n'
        + "".join(f"{index}.0 {value}\n" for index, value in enumerate(values))
    )
    return path


def test_calculate_bar_uses_cross_hamiltonian_energies(tmp_path):
    paths = {
        (0, 0): _energy(tmp_path / "00.xvg", [0.0, 1.0, 2.0]),
        (0, 1): _energy(tmp_path / "01.xvg", [2.0, 3.0, 4.0]),
        (1, 1): _energy(tmp_path / "11.xvg", [5.0, 6.0, 7.0]),
        (1, 0): _energy(tmp_path / "10.xvg", [3.0, 4.0, 5.0]),
    }

    pairs, total, uncertainty = analyze_fep_rest._calculate_bar(
        paths,
        2,
        300.0,
        0.0,
        None,
        1,
    )

    assert total == pytest.approx(2.0)
    assert uncertainty == pytest.approx(0.0, abs=1e-7)
    assert pairs[0]["forward_samples"] == 3


def test_rerun_energy_uses_cache(tmp_path):
    output_dir = tmp_path / "rerun"
    output_dir.mkdir()
    cached = output_dir / "potential.xvg"
    cached.write_text("0 1\n")

    result = analyze_fep_rest._rerun_energy(
        type("Args", (), {"force": False, "gmx": "gmx"})(),
        Path("missing.trr"),
        Path("missing.tpr"),
        output_dir,
    )

    assert result == cached
