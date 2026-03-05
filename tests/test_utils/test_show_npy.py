"""
utils/show_npy のユニットテスト

run() が .npy ファイルの内容と shape を標準出力に表示することをテストする。
"""

import types

import numpy as np


class TestShowNpyRun:
    def _make_args(self, npy_path):
        return types.SimpleNamespace(npy=str(npy_path))

    def test_prints_array_content(self, tmp_path, capsys):
        """配列の内容が標準出力に表示されること"""
        from src.utils.show_npy import run

        arr = np.array([1.0, 2.0, 3.0])
        npy = tmp_path / "test.npy"
        np.save(str(npy), arr)

        run(self._make_args(npy))
        captured = capsys.readouterr()
        assert "1." in captured.out or "1.0" in captured.out

    def test_prints_shape(self, tmp_path, capsys):
        """配列の shape が標準出力に表示されること"""
        from src.utils.show_npy import run

        arr = np.zeros((3, 4))
        npy = tmp_path / "shape_test.npy"
        np.save(str(npy), arr)

        run(self._make_args(npy))
        captured = capsys.readouterr()
        assert "(3, 4)" in captured.out

    def test_2d_array(self, tmp_path, capsys):
        """2D 配列でも正常に動作すること"""
        from src.utils.show_npy import run

        arr = np.arange(6).reshape(2, 3)
        npy = tmp_path / "2d.npy"
        np.save(str(npy), arr)

        run(self._make_args(npy))
        captured = capsys.readouterr()
        assert "(2, 3)" in captured.out

    def test_scalar_array(self, tmp_path, capsys):
        """スカラー配列でも正常に動作すること"""
        from src.utils.show_npy import run

        arr = np.float64(42.0)
        npy = tmp_path / "scalar.npy"
        np.save(str(npy), arr)

        run(self._make_args(npy))
        captured = capsys.readouterr()
        assert len(captured.out.strip()) > 0
