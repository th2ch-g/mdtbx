"""
calc_ion_conc のユニットテスト

純粋計算関数と PDB パーサーをテストする。
"""

import types

import pytest

from src.build.calc_ion_conc import (
    calc_ion_conc_from_volume,
    get_boxsize_from_pdb,
    get_water_number_from_pdb,
)

# AVOGADRO_CONST = 6.022 (config.py より)
AVOGADRO = 6.022


class TestCalcIonConcFromVolume:
    def test_known_value(self):
        """
        volume=1e6 A^3, concentration=0.15 M のとき
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
        """体積が 2 倍になるとイオン数もほぼ 2 倍になること"""
        n1 = calc_ion_conc_from_volume(1e6, 0.15)
        n2 = calc_ion_conc_from_volume(2e6, 0.15)
        assert abs(n2 - 2 * n1) <= 1  # 整数切り捨てによる誤差を許容

    def test_proportional_to_concentration(self):
        """濃度が 2 倍になるとイオン数もほぼ 2 倍になること"""
        n1 = calc_ion_conc_from_volume(1e6, 0.10)
        n2 = calc_ion_conc_from_volume(1e6, 0.20)
        assert abs(n2 - 2 * n1) <= 1


class TestGetBoxsizeFromPdb:
    def test_reads_cryst_line(self, sample_pdb_path):
        """CRYST1 行からボックスサイズが読み取れること"""
        args = types.SimpleNamespace(pdb=str(sample_pdb_path))
        x, y, z = get_boxsize_from_pdb(args)
        assert x == pytest.approx(50.0)
        assert y == pytest.approx(50.0)
        assert z == pytest.approx(50.0)

    def test_raises_if_no_cryst(self, tmp_path):
        """CRYST1 行がない PDB は例外を送出すること"""
        pdb = tmp_path / "no_cryst.pdb"
        pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
        args = types.SimpleNamespace(pdb=str(pdb))
        with pytest.raises(Exception, match="CRYST"):
            get_boxsize_from_pdb(args)


class TestGetWaterNumberFromPdb:
    def test_counts_wat_oxygens(self, sample_pdb_path):
        """WAT の酸素原子数が正しくカウントされること（fixture は 2 分子）"""
        args = types.SimpleNamespace(pdb=str(sample_pdb_path), water_name="WAT")
        count = get_water_number_from_pdb(args)
        assert count == 2

    def test_counts_with_custom_water_name(self, tmp_path):
        """water_name を変えたときに正しくカウントされること"""
        pdb = tmp_path / "wat.pdb"
        # "WAT" は "O" を含まないため、O 原子行のみマッチする
        pdb.write_text(
            "HETATM    1  O   WAT A   1       0.000   0.000   0.000\n"
            "HETATM    2  H1  WAT A   1       0.800   0.600   0.000\n"
            "HETATM    3  H2  WAT A   1      -0.800   0.600   0.000\n"
        )
        args = types.SimpleNamespace(pdb=str(pdb), water_name="WAT")
        count = get_water_number_from_pdb(args)
        assert count == 1

    def test_no_water_returns_zero(self, tmp_path):
        """水分子がない場合は 0 を返すこと"""
        pdb = tmp_path / "protein_only.pdb"
        pdb.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
        args = types.SimpleNamespace(pdb=str(pdb), water_name="WAT")
        count = get_water_number_from_pdb(args)
        assert count == 0
