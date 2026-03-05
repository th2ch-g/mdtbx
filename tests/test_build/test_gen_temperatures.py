"""
gen_temperatures のユニットテスト

純粋計算関数（calc_mu）とバリデーションロジック（run）をテストする。
"""

import types

import pytest

from src.build.gen_temperatures import A0, A1, B0, B1, calc_mu, run


class TestCalcMu:
    def test_known_value(self):
        """calc_mu の数値が手計算と一致すること"""
        nw, np_val, temp, fener = 100, 50, 300.0, 0.0
        expected = (A0 + A1 * temp) * nw + (B0 + B1 * temp) * np_val - temp * fener
        assert calc_mu(nw, np_val, temp, fener) == pytest.approx(expected)

    def test_zero_system(self):
        """nw=0, np=0 のとき 0 になること"""
        result = calc_mu(0, 0, 300.0, 0.0)
        assert result == 0.0

    def test_temperature_effect(self):
        """温度が高いほど calc_mu の値が変化すること（単調ではないが変わること）"""
        mu1 = calc_mu(1000, 500, 300.0, 0.0)
        mu2 = calc_mu(1000, 500, 400.0, 0.0)
        assert mu1 != mu2


class TestGenTemperaturesRun:
    def _base_args(self, **kwargs):
        defaults = dict(
            pdes=0.2,
            tlow=300.0,
            thigh=350.0,
            nw=1000,
            np=500,
            tol=1e-3,
            pc=0,
            wc=3,
            hff=0,
            vs=0,
            alg=0,
        )
        defaults.update(kwargs)
        return types.SimpleNamespace(**defaults)

    def test_generates_temperature_ladder(self, capsys):
        """正常な入力で温度ラダーが出力されること"""
        args = self._base_args()
        run(args)
        captured = capsys.readouterr()
        assert "Temperature" in captured.out
        assert "300.00" in captured.out

    def test_first_temp_equals_tlow(self, capsys):
        """最初の温度が tlow であること"""
        args = self._base_args(tlow=310.0, thigh=360.0)
        run(args)
        captured = capsys.readouterr()
        assert "310.00" in captured.out

    def test_invalid_pdes_raises(self):
        """0〜1 の範囲外の pdes は ValueError を送出すること"""
        args = self._base_args(pdes=1.5)
        with pytest.raises(ValueError, match="Pdes"):
            run(args)

    def test_thigh_must_be_greater_than_tlow(self):
        """thigh <= tlow のとき ValueError を送出すること"""
        args = self._base_args(tlow=350.0, thigh=300.0)
        with pytest.raises(ValueError):
            run(args)

    def test_tlow_must_be_positive(self):
        """tlow が 0 以下のとき ValueError を送出すること"""
        args = self._base_args(tlow=0.0)
        with pytest.raises(ValueError):
            run(args)

    def test_zero_protein_atoms_raises(self):
        """np=0 のとき ValueError を送出すること"""
        args = self._base_args(np=0)
        with pytest.raises(ValueError, match="protein atoms"):
            run(args)

    def test_nvt_not_supported(self):
        """alg=1（NVT）は未対応なので ValueError を送出すること"""
        args = self._base_args(alg=1)
        with pytest.raises(ValueError, match="constant volume"):
            run(args)
