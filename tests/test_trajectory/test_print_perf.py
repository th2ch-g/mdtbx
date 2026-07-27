"""
trajectory/print_perf unit tests

Tests the pure log parsing logic of parse_log_file(). It depends on no
external tool and not on gmx.
"""

import pytest

from src.trajectory.print_perf import parse_log_file

# A minimal sample of a GROMACS log
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
        """The ns/day value is read from the `Performance:` line"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert data["performance"] == pytest.approx(123.45)

    def test_parses_version(self, tmp_path):
        """The version string is read from `GROMACS version:`"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert data["version"] == "2023.3"

    def test_parses_executable(self, tmp_path):
        """The path is read from `Executable:`"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert data["executable"] == "/usr/local/bin/gmx"

    def test_parses_hostname(self, tmp_path):
        """The host name is read from `Hardware detected on host`"""
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
        """The returned command string has the `-deffnm` option stripped"""
        log = tmp_path / "prd.log"
        log.write_text(_SAMPLE_LOG)
        data = parse_log_file(log)
        assert "-deffnm" not in data["cmd"]
        assert "gmx mdrun" in data["cmd"]

    def test_nonexistent_file_returns_none(self, tmp_path):
        """A missing file returns None"""
        data = parse_log_file(tmp_path / "nonexistent.log")
        assert data is None

    def test_empty_log_has_no_performance(self, tmp_path):
        """performance is None for a log without a Performance line"""
        log = tmp_path / "empty.log"
        log.write_text("Some GROMACS output without performance line\n")
        data = parse_log_file(log)
        assert data is not None
        assert data["performance"] is None

    def test_default_values_for_missing_fields(self, tmp_path):
        """Fields that are not found default to 'N/A'"""
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
