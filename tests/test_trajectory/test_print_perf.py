"""
trajectory/print_perf のユニットテスト

parse_log_file() の純粋なログパースロジックをテストする。
外部ツール・gmx に依存しない。
"""

import pytest

from src.trajectory.print_perf import parse_log_file

# GROMACS ログの最小限サンプル
_SAMPLE_LOG = """\
                      GROMACS version:    2023.3
                         Executable:   /usr/local/bin/gmx
Hardware detected on host myhost01:
GPU info:
    Number of GPUs detected: 2
    GPU 0: NVIDIA A100
CPU info:
    Vendor: Intel
    Model name:  Xeon Gold 6338
Command line:
  gmx mdrun -deffnm prd -ntmpi 4 -ntomp 8

               Performance:       123.45          0.194
"""


class TestParseLogFile:
    def test_parses_performance(self, tmp_path):
        """`Performance:` 行から ns/day 値を正しく取得すること"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert data["performance"] == pytest.approx(123.45)

    def test_parses_version(self, tmp_path):
        """`GROMACS version:` から版数文字列を取得すること"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert data["version"] == "2023.3"

    def test_parses_executable(self, tmp_path):
        """`Executable:` からパスを取得すること"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert data["executable"] == "/usr/local/bin/gmx"

    def test_parses_hostname(self, tmp_path):
        """`Hardware detected on host` からホスト名を取得すること"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert data["hostname"] == "myhost01"

    def test_cmd_strips_deffnm(self, tmp_path):
        """`-deffnm` オプションが除去されたコマンド文字列を返すこと"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert "-deffnm" not in data["cmd"]
        assert "gmx mdrun" in data["cmd"]

    def test_nonexistent_file_returns_none(self, tmp_path):
        """存在しないファイルは None を返すこと"""
        data = parse_log_file(tmp_path / "nonexistent.log")
        assert data is None

    def test_empty_log_has_no_performance(self, tmp_path):
        """Performance 行がないログでは performance が None になること"""
        log = tmp_path / "empty.log"
        log.write_text("Some GROMACS output without performance line\n")
        data = parse_log_file(log)
        assert data is not None
        assert data["performance"] is None

    def test_default_values_for_missing_fields(self, tmp_path):
        """フィールドが見つからない場合は 'N/A' がデフォルトであること"""
        log = tmp_path / "minimal.log"
        log.write_text("               Performance:       50.0   0.5\n")
        data = parse_log_file(log)
        assert data["version"] == "N/A"
        assert data["hostname"] == "N/A"
