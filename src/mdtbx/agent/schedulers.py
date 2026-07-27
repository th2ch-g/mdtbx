"""Batch scheduler adapters with normalized probe and status output."""

from __future__ import annotations

import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from .model import SCHEMA_VERSION, utc_now

NORMAL_STATES = {
    "queued",
    "running",
    "succeeded",
    "failed",
    "cancelled",
    "unknown",
}


def _run(argv: list[str], *, check: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(
        argv,
        check=check,
        capture_output=True,
        text=True,
        timeout=30,
    )


def _available(*commands: str) -> bool:
    return all(shutil.which(command) is not None for command in commands)


def detect_scheduler() -> tuple[str | None, list[str]]:
    matches = []
    if _available("sbatch", "squeue"):
        matches.append("slurm")
    if _available("qsub", "qstat"):
        matches.append("age")
    if _available("pjsub", "pjstat"):
        matches.append("pjm")
    return (matches[0] if len(matches) == 1 else None), matches


def _walltime(seconds: int) -> str:
    hours, remainder = divmod(int(seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def _directive_lines(resource: dict[str, Any]) -> list[str]:
    return list(resource.get("scheduler_options", {}).get("extra_directives", []))


class Scheduler(ABC):
    name: str

    @abstractmethod
    def probe(self, profile: dict[str, Any] | None = None) -> dict[str, Any]:
        raise NotImplementedError

    @abstractmethod
    def directives(
        self,
        step: dict[str, Any],
        resource: dict[str, Any],
        profile: dict[str, Any],
    ) -> list[str]:
        raise NotImplementedError

    @abstractmethod
    def submit(
        self,
        script: Path,
        *,
        dependencies: list[str],
        resource: dict[str, Any],
    ) -> str:
        raise NotImplementedError

    @abstractmethod
    def status(self, job_id: str) -> dict[str, Any]:
        raise NotImplementedError


class SlurmScheduler(Scheduler):
    name = "slurm"

    def probe(self, profile=None):
        commands = [
            ["sinfo", "-h", "-o", "%P|%a|%l|%D|%C|%G"],
            ["scontrol", "show", "node", "-o"],
        ]
        results = []
        for argv in commands:
            result = _run(argv)
            results.append(
                {
                    "argv": argv,
                    "returncode": result.returncode,
                    "stdout": result.stdout[:200000],
                    "stderr": result.stderr[:20000],
                }
            )
        return {
            "schema_version": SCHEMA_VERSION,
            "scheduler": self.name,
            "probed_at": utc_now(),
            "commands": results,
            "resources": _parse_slurm_nodes(results[1]["stdout"]),
        }

    def directives(self, step, resource, profile):
        request = step["resources"]
        options = resource.get("scheduler_options", {})
        lines = [
            f"#SBATCH --job-name={step['id'][:64]}",
            f"#SBATCH --nodes={request['nodes']}",
            f"#SBATCH --time={_walltime(request['walltime_seconds'])}",
            f"#SBATCH --cpus-per-task={request['cpus_per_node']}",
            f"#SBATCH --mem={request['memory_mb']}M",
        ]
        if request["gpus_per_node"]:
            gres = options.get("gres_name", "gpu")
            lines.append(f"#SBATCH --gres={gres}:{request['gpus_per_node']}")
        for key, option in (
            ("partition", "--partition"),
            ("account", "--account"),
            ("qos", "--qos"),
            ("constraint", "--constraint"),
        ):
            if options.get(key):
                lines.append(f"#SBATCH {option}={options[key]}")
        lines.extend(_directive_lines(resource))
        return lines

    def submit(self, script, *, dependencies, resource):
        argv = ["sbatch", "--parsable"]
        if dependencies:
            argv.append(f"--dependency=afterok:{':'.join(dependencies)}")
        argv.append(str(script))
        result = _run(argv, check=True)
        job_id = result.stdout.strip().split(";", 1)[0]
        if not job_id:
            raise ValueError("sbatch did not return a job ID")
        return job_id

    def status(self, job_id):
        result = _run(["squeue", "-h", "-j", job_id, "-o", "%T|%R|%M"])
        line = result.stdout.strip().splitlines()
        if line:
            state, reason, elapsed = (line[0].split("|") + ["", ""])[:3]
            return _status(job_id, _slurm_state(state), state, reason, elapsed)
        result = _run(
            [
                "sacct",
                "-n",
                "-P",
                "-j",
                job_id,
                "--format=State,ExitCode,Elapsed,NodeList",
            ]
        )
        line = result.stdout.strip().splitlines()
        if not line:
            return _status(job_id, "unknown", "", "not found", "")
        fields = line[0].split("|")
        state = fields[0].split()[0]
        return _status(
            job_id,
            _slurm_state(state),
            state,
            fields[1] if len(fields) > 1 else "",
            fields[2] if len(fields) > 2 else "",
        )


class AgeScheduler(Scheduler):
    name = "age"

    def probe(self, profile=None):
        commands = [["qconf", "-sql"], ["qstat", "-xml"]]
        results = []
        for argv in commands:
            result = _run(argv)
            results.append(
                {
                    "argv": argv,
                    "returncode": result.returncode,
                    "stdout": result.stdout[:200000],
                    "stderr": result.stderr[:20000],
                }
            )
        return {
            "schema_version": SCHEMA_VERSION,
            "scheduler": self.name,
            "probed_at": utc_now(),
            "commands": results,
            "queues": [
                line for line in results[0]["stdout"].splitlines() if line.strip()
            ],
            "resources": profile.get("resources", []) if profile else [],
        }

    def directives(self, step, resource, profile):
        request = step["resources"]
        options = resource.get("scheduler_options", {})
        lines = [
            "#$ -cwd",
            f"#$ -N {step['id'][:64]}",
            f"#$ -l h_rt={_walltime(request['walltime_seconds'])}",
        ]
        resource_name = options.get("resource_name")
        if resource_name:
            lines.append(f"#$ -l {resource_name}={request['nodes']}")
        if options.get("queue"):
            lines.append(f"#$ -q {options['queue']}")
        if options.get("parallel_environment") and request["cpus_per_node"] > 1:
            slots = request["nodes"] * request["cpus_per_node"]
            lines.append(f"#$ -pe {options['parallel_environment']} {slots}")
        lines.extend(_directive_lines(resource))
        return lines

    def submit(self, script, *, dependencies, resource):
        options = resource.get("scheduler_options", {})
        argv = ["qsub", "-terse"]
        if options.get("group"):
            argv.extend(["-g", options["group"]])
        if dependencies:
            argv.extend(["-hold_jid", ",".join(dependencies)])
        argv.append(str(script))
        result = _run(argv, check=True)
        match = re.search(r"\d+", result.stdout)
        if not match:
            raise ValueError("qsub did not return a job ID")
        return match.group(0)

    def status(self, job_id):
        result = _run(["qstat", "-xml"])
        if result.returncode == 0 and result.stdout.strip():
            try:
                root = ET.fromstring(result.stdout)
                for job in root.findall(".//job_list"):
                    number = job.findtext("JB_job_number")
                    if number != str(job_id):
                        continue
                    state = job.findtext("state", default="")
                    return _status(job_id, _age_state(state), state, "", "")
            except ET.ParseError:
                pass
        result = _run(["qacct", "-j", str(job_id)])
        if result.returncode != 0 or not result.stdout.strip():
            return _status(job_id, "unknown", "", "not found", "")
        values = {}
        for line in result.stdout.splitlines():
            fields = line.split(None, 1)
            if len(fields) == 2:
                values[fields[0]] = fields[1]
        failed = values.get("failed", "0") != "0"
        exit_status = values.get("exit_status", "0")
        state = "failed" if failed or exit_status != "0" else "succeeded"
        return _status(
            job_id, state, state, exit_status, values.get("ru_wallclock", "")
        )


class PjmScheduler(Scheduler):
    name = "pjm"

    def probe(self, profile=None):
        commands = [["pjstat", "--rsc"], ["pjstat", "--limit"]]
        results = []
        for argv in commands:
            result = _run(argv)
            results.append(
                {
                    "argv": argv,
                    "returncode": result.returncode,
                    "stdout": result.stdout[:200000],
                    "stderr": result.stderr[:20000],
                }
            )
        return {
            "schema_version": SCHEMA_VERSION,
            "scheduler": self.name,
            "probed_at": utc_now(),
            "commands": results,
            "resources": profile.get("resources", []) if profile else [],
        }

    def directives(self, step, resource, profile):
        request = step["resources"]
        options = resource.get("scheduler_options", {})
        lines = [
            f"#PJM -N {step['id'][:64]}",
            f'#PJM -L "node={request["nodes"]}"',
            f'#PJM -L "elapse={_walltime(request["walltime_seconds"])}"',
        ]
        if options.get("resource_group"):
            lines.append(f'#PJM -L "rscgrp={options["resource_group"]}"')
        if options.get("group"):
            lines.append(f"#PJM -g {options['group']}")
        if options.get("max_proc_per_node"):
            lines.append(
                f'#PJM --mpi "max-proc-per-node={options["max_proc_per_node"]}"'
            )
        lines.extend(_directive_lines(resource))
        return lines

    def submit(self, script, *, dependencies, resource):
        options = resource.get("scheduler_options", {})
        argv = ["pjsub"]
        if dependencies:
            template = options.get("dependency_args")
            if not isinstance(template, list):
                raise ValueError(
                    "PJM dependencies require scheduler_options.dependency_args"
                )
            for item in template:
                argv.append(item.format(job_ids=",".join(dependencies)))
        argv.append(str(script))
        result = _run(argv, check=True)
        match = re.search(r"\bJob\s+(\d+)\s+submitted\b", result.stdout + result.stderr)
        if not match:
            raise ValueError("pjsub did not return a job ID")
        return match.group(1)

    def status(self, job_id):
        result = _run(["pjstat", "-v", str(job_id)])
        if result.returncode != 0 or not result.stdout.strip():
            return _status(job_id, "unknown", "", "not found", "")
        for line in result.stdout.splitlines():
            if not line.strip().startswith(str(job_id)):
                continue
            fields = line.split()
            raw = fields[3] if len(fields) > 3 else ""
            return _status(job_id, _pjm_state(raw), raw, "", "")
        return _status(job_id, "unknown", "", "not found", "")


def scheduler(name: str) -> Scheduler:
    mapping = {
        "slurm": SlurmScheduler,
        "age": AgeScheduler,
        "pjm": PjmScheduler,
    }
    try:
        return mapping[name]()
    except KeyError as error:
        raise ValueError(f"Unsupported batch scheduler: {name}") from error


def _status(job_id, state, raw_state, reason, elapsed):
    if state not in NORMAL_STATES:
        state = "unknown"
    return {
        "schema_version": SCHEMA_VERSION,
        "job_id": str(job_id),
        "state": state,
        "raw_state": raw_state,
        "reason": reason,
        "elapsed": elapsed,
        "checked_at": utc_now(),
    }


def _slurm_state(value: str) -> str:
    value = value.upper().split("+", 1)[0]
    if value in {"PENDING", "CONFIGURING", "SUSPENDED"}:
        return "queued"
    if value in {"RUNNING", "COMPLETING"}:
        return "running"
    if value == "COMPLETED":
        return "succeeded"
    if value.startswith("CANCELLED"):
        return "cancelled"
    if value in {
        "FAILED",
        "TIMEOUT",
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "PREEMPTED",
        "BOOT_FAIL",
        "DEADLINE",
    }:
        return "failed"
    return "unknown"


def _age_state(value: str) -> str:
    value = value.lower()
    if "r" in value and "qw" not in value:
        return "running"
    if "qw" in value or "h" in value or "t" in value:
        return "queued"
    if "e" in value:
        return "failed"
    return "unknown"


def _pjm_state(value: str) -> str:
    value = value.upper()
    if value in {"ACC", "HLD", "QUE", "RNA"}:
        return "queued"
    if value in {"RUN"}:
        return "running"
    if value in {"EXT"}:
        return "succeeded"
    if value in {"CCL"}:
        return "cancelled"
    if value in {"ERR", "RJT"}:
        return "failed"
    return "unknown"


def _parse_slurm_nodes(text: str) -> list[dict[str, Any]]:
    groups: dict[tuple[int, int, int, str], int] = {}
    for line in text.splitlines():
        values = {}
        for field in line.split():
            key, separator, value = field.partition("=")
            if separator:
                values[key] = value
        try:
            cpus = int(values.get("CPUTot", 0))
            memory = int(values.get("RealMemory", 0))
        except ValueError:
            continue
        gres = values.get("Gres", "")
        gpu_match = re.search(r"gpu(?::([^:,()]+))?:(\d+)", gres)
        gpu_type = gpu_match.group(1) if gpu_match else ""
        gpus = int(gpu_match.group(2)) if gpu_match else 0
        state = values.get("State", "UNKNOWN").split("+", 1)[0]
        key = (cpus, gpus, memory, gpu_type or "")
        groups[key] = groups.get(key, 0) + (0 if state == "DOWN" else 1)
    return [
        {
            "name": f"{gpu_type or 'cpu'}-{cpus}c-{memory}m",
            "count": count,
            "capabilities": {
                "cpus_per_node": cpus,
                "gpus_per_node": gpus,
                "memory_mb_per_node": memory,
                "gpu_type": gpu_type or None,
            },
        }
        for (cpus, gpus, memory, gpu_type), count in sorted(groups.items())
    ]
