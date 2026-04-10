"""
trajectory/opt_perf のユニットテスト

小さなダミーコマンドを grid sampler で最適化し、ベスト条件を検証する。
"""

import json
import types


class TestOptPerfRun:
    def _make_args(self, script_path, tmp_path):
        return types.SimpleNamespace(
            command_template=None,
            mdrun_command=(
                f"sh {script_path} {{log_path}} {{n_gpu}} {{n_core}} "
                "{ntomp} {ntmpi} {gpu_ids}"
            ),
            mpi_launcher=None,
            mpi_mode="auto",
            log_name="md.log",
            workdir=str(tmp_path / "trials"),
            sampler="grid",
            seed=0,
            n_trials=None,
            n_gpu=[1, 2],
            n_core=[4, 8],
            ntomp=[1, 2],
            ntmpi=[1, 4],
            gpu_env="cuda_visible_devices",
            output=str(tmp_path / "best.json"),
            history_output=str(tmp_path / "history.csv"),
        )

    def test_run_finds_best_trial(self, tmp_path):
        from src.trajectory.opt_perf import run

        script = tmp_path / "fake_mdrun.sh"
        script.write_text(
            "#!/bin/sh\n"
            "log_path=$1\n"
            "n_gpu=$2\n"
            "n_core=$3\n"
            "ntomp=$4\n"
            "ntmpi=$5\n"
            "gpu_ids=$6\n"
            "perf=$((n_gpu * 1000 + n_core * 100 + ntomp * 10 + ntmpi + ${#gpu_ids}))\n"
            'cat > "$log_path" <<EOF\n'
            "               Performance:       $perf          0.194\n"
            "EOF\n"
        )
        script.chmod(0o755)

        args = self._make_args(script, tmp_path)
        run(args)

        best = json.loads((tmp_path / "best.json").read_text())
        assert best["params"] == {
            "n_core": 8,
            "n_gpu": 2,
            "ntmpi": 4,
            "ntomp": 2,
        }
        assert (tmp_path / "history.csv").exists()

    def test_build_command_template_appends_missing_mdrun_args(self):
        from src.trajectory.opt_perf import _build_command_template

        args = types.SimpleNamespace(
            command_template=None,
            mdrun_command="gmx mdrun -deffnm prd",
            mpi_launcher=None,
            mpi_mode="auto",
        )

        command = _build_command_template(args)
        assert "{log_path}" in command
        assert "{ntomp}" in command
        assert "{ntmpi}" in command

    def test_build_command_template_supports_external_mpi_launcher(self):
        from src.trajectory.opt_perf import _build_command_template

        args = types.SimpleNamespace(
            command_template=None,
            mdrun_command="gmx_mpi mdrun -deffnm prd -gpu_id {gpu_ids}",
            mpi_launcher="mpirun -np {ntmpi}",
            mpi_mode="auto",
        )

        command = _build_command_template(args)
        assert command.startswith("mpirun -np {ntmpi} gmx_mpi mdrun")
        assert command.count("{ntmpi}") == 1
        assert "{log_path}" in command
        assert "{ntomp}" in command
