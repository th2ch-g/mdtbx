import pytest

from mdtbx.analysis.analyze_rbfe import calculate


def test_calculate_uses_complex_minus_solvent_and_ki_ratio():
    result = calculate(
        complex_delta_g=5.0,
        complex_uncertainty=0.6,
        solvent_delta_g=2.0,
        solvent_uncertainty=0.8,
        ki_a=0.32,
        ki_a_error=0.01,
        ki_b=1.3,
        ki_b_error=0.52,
        temperature=300.0,
    )

    assert result["calculated_delta_delta_g_kj_mol"] == pytest.approx(3.0)
    assert result["calculated_uncertainty_kj_mol"] == pytest.approx(1.0)
    assert result["experimental_delta_delta_g_kj_mol"] == pytest.approx(
        3.4966, abs=1e-4
    )
    assert result["experimental_uncertainty_kj_mol"] == pytest.approx(1.0008, abs=1e-4)
