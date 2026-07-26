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

    def test_collects_all_gpu_descriptions(self, tmp_path):
        log = tmp_path / "multi-gpu.log"
        log.write_text(
            "GPU info:\n"
            "    Number of GPUs detected: 2\n"
            "    GPU 0: NVIDIA A100\n"
            "    GPU 1: NVIDIA A100\n"
            "CPU info:\n"
        )

        data = parse_log_file(log)

        assert data["n_GPU"] == "2"
        assert data["GPU_info"] == "GPU 0: NVIDIA A100; GPU 1: NVIDIA A100"

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

    def test_malformed_performance_preserves_other_log_data(self, tmp_path):
        log = tmp_path / "malformed.log"
        log.write_text("GROMACS version: 2024.1\nPerformance: not-a-number\n")

        data = parse_log_file(log)
        assert data["version"] == "2024.1"
        assert data["performance"] is None

    @pytest.mark.parametrize("value", ["nan", "inf", "-1.0"])
    def test_invalid_performance_is_ignored(self, tmp_path, value):
        log = tmp_path / "invalid.log"
        log.write_text(f"Performance: {value}\n")

        assert parse_log_file(log)["performance"] is None

    def test_deffnm_equals_syntax_is_removed(self, tmp_path):
        log = tmp_path / "command.log"
        log.write_text("Command line:\n  gmx mdrun -deffnm=prd -nt 8\n")

        assert parse_log_file(log)["cmd"] == "gmx mdrun -nt 8"
