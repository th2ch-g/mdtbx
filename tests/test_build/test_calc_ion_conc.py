"""
calc_ion_conc unit tests

Tests the pure computation helpers and the PDB parser.
"""

import types

import pytest

from src.build.calc_ion_conc import (
    calc_ion_conc_from_volume,
    get_boxsize_from_pdb,
    get_water_number_from_pdb,
)

# AVOGADRO_CONST = 6.022 (from config.py)
AVOGADRO = 6.022


class TestCalcIonConcFromVolume:
    def test_known_value(self):
        """
        For volume=1e6 A^3 and concentration=0.15 M,
        ionnum = 1e6 * 0.15 * 6.022 // 10000 = 90
        """
        result = calc_ion_conc_from_volume(1e6, 0.15)
        assert result == 90

    def test_zero_concentration(self):
        result = calc_ion_conc_from_volume(1e6, 0.0)
        assert result == 0

    def test_returns_int(self):
        result = calc_ion_conc_from_volume(1e6, 0.15)
        assert isinstance(result, int)

    def test_proportional_to_volume(self):
        """Doubling the volume roughly doubles the ion count"""
        n1 = calc_ion_conc_from_volume(1e6, 0.15)
        n2 = calc_ion_conc_from_volume(2e6, 0.15)
        assert abs(n2 - 2 * n1) <= 1  # Allow for integer truncation

    def test_proportional_to_concentration(self):
        """Doubling the concentration roughly doubles the ion count"""
        n1 = calc_ion_conc_from_volume(1e6, 0.10)
        n2 = calc_ion_conc_from_volume(1e6, 0.20)
        assert abs(n2 - 2 * n1) <= 1

    @pytest.mark.parametrize(
        ("volume", "concentration"),
        [(-1.0, 0.15), (1e6, -0.15), (float("nan"), 0.15)],
    )
    def test_invalid_inputs_raise(self, volume, concentration):
        with pytest.raises(ValueError):
            calc_ion_conc_from_volume(volume, concentration)


class TestGetBoxsizeFromPdb:
    def test_reads_cryst_line(self, sample_pdb_path):
        """The box size is read from the CRYST1 record"""
        args = types.SimpleNamespace(pdb=str(sample_pdb_path))
        x, y, z = get_boxsize_from_pdb(args)
        assert x == pytest.approx(50.0)
        assert y == pytest.approx(50.0)
        assert z == pytest.approx(50.0)

    def test_raises_if_no_cryst(self, tmp_path):
        """A PDB without a CRYST1 record raises"""
        pdb = tmp_path / "no_cryst.pdb"
        pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
        args = types.SimpleNamespace(pdb=str(pdb))
        with pytest.raises(Exception, match="CRYST"):
            get_boxsize_from_pdb(args)

    def test_ignores_cryst_text_outside_cryst1_record(self, tmp_path):
        pdb = tmp_path / "remark.pdb"
        pdb.write_text("REMARK CRYST values are unavailable\n")
        args = types.SimpleNamespace(pdb=str(pdb))

        with pytest.raises(ValueError, match="CRYST1"):
            get_boxsize_from_pdb(args)


class TestGetWaterNumberFromPdb:
    def test_counts_wat_oxygens(self, sample_pdb_path):
        """WAT oxygen atoms are counted correctly (the fixture holds 2 molecules)"""
        args = types.SimpleNamespace(pdb=str(sample_pdb_path), water_name="WAT")
        count = get_water_number_from_pdb(args)
        assert count == 2

    def test_counts_with_custom_water_name(self, tmp_path):
        """Counting is still correct after changing water_name"""
        pdb = tmp_path / "wat.pdb"
        # "WAT" does not contain "O", so only the O atom lines match
        pdb.write_text(
            "HETATM    1  O   WAT A   1       0.000   0.000   0.000\n"
            "HETATM    2  H1  WAT A   1       0.800   0.600   0.000\n"
            "HETATM    3  H2  WAT A   1      -0.800   0.600   0.000\n"
        )
        args = types.SimpleNamespace(pdb=str(pdb), water_name="WAT")
        count = get_water_number_from_pdb(args)
        assert count == 1

    def test_no_water_returns_zero(self, tmp_path):
        """Returns 0 when there is no water"""
        pdb = tmp_path / "protein_only.pdb"
        pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
        args = types.SimpleNamespace(pdb=str(pdb), water_name="WAT")
        count = get_water_number_from_pdb(args)
        assert count == 0
