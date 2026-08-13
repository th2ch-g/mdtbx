import argparse
import json
import math
from pathlib import Path

import numpy as np

from ..config import GAS_CONSTANT_KJ
from ..logger import generate_logger
from ..utils.convergence import convergence_ranges
from ..utils.fep import KJ_TO_KCAL, load_fep_manifest
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "analyze_fep_rest",
        help="Analyze PLUMED FEP/REST replicas by cross-Hamiltonian reruns",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--path",
        default="fep_rest",
        help="FEP/REST setup directory or manifest",
    )
    parser.add_argument(
        "-b", "--begin", type=float, default=0.0, help="Begin time [ps]"
    )
    parser.add_argument("-e", "--end", type=float, help="End time [ps]")
    parser.add_argument(
        "--subsample",
        type=int,
        default=1,
        help="Use every Nth energy sample",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output JSON file",
    )
    parser.add_argument(
        "--gmx",
        help="GROMACS executable; defaults to the setup executable",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate cached cross-Hamiltonian energies",
    )
    parser.add_argument(
        "--gpu-offload",
        action="store_true",
        help="Enable GPU nonbonded and PME offload for reruns",
    )
    parser.add_argument(
        "--ntomp",
        type=int,
        default=1,
        help="OpenMP thread count for each rerun",
    )
    parser.add_argument(
        "--convergence-blocks",
        type=int,
        default=0,
        help="Number of equal-duration convergence blocks; 0 disables analysis",
    )
    parser.set_defaults(func=run)


def _parse_energy_xvg(path):
    times = []
    values = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line[0] in {"#", "@"}:
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            times.append(float(fields[0]))
            values.append(float(fields[1]))
    if not times:
        raise ValueError(f"No energy samples found in {path}")
    return np.asarray(times), np.asarray(values)


def _filtered_energy(path, begin, end, subsample):
    times, values = _parse_energy_xvg(path)
    mask = times >= begin
    if end is not None:
        mask &= times <= end
    return times[mask][::subsample], values[mask][::subsample]


def _aligned_difference(reference_path, target_path, begin, end, subsample):
    reference_times, reference = _filtered_energy(
        reference_path,
        begin,
        end,
        subsample,
    )
    target_times, target = _filtered_energy(
        target_path,
        begin,
        end,
        subsample,
    )
    target_by_time = {
        round(float(time), 6): float(value) for time, value in zip(target_times, target)
    }
    differences = [
        target_by_time[round(float(time), 6)] - float(value)
        for time, value in zip(reference_times, reference)
        if round(float(time), 6) in target_by_time
    ]
    if len(differences) < 2:
        raise ValueError("FEP/REST analysis requires at least 2 aligned samples")
    return np.asarray(differences)


def _rerun_energy(args, trajectory, tpr, output_dir):
    potential = output_dir / "potential.xvg"
    if potential.is_file() and not args.force:
        cache_time = potential.stat().st_mtime_ns
        if all(
            source.is_file() and source.stat().st_mtime_ns <= cache_time
            for source in (trajectory, tpr)
        ):
            return potential
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        args.gmx,
        "mdrun",
        "-s",
        str(tpr),
        "-rerun",
        str(trajectory),
        "-e",
        "rerun.edr",
        "-g",
        "rerun.log",
        "-deffnm",
        "rerun",
        "-ntomp",
        str(args.ntomp),
    ]
    if getattr(args, "gpu_offload", False):
        command.extend(
            [
                "-pin",
                "on",
                "-nb",
                "gpu",
                "-pme",
                "gpu",
                "-pmefft",
                "gpu",
                "-bonded",
                "gpu",
                "-update",
                "cpu",
            ]
        )
    run_cmd(command, cwd=output_dir)
    run_cmd(
        [
            args.gmx,
            "energy",
            "-f",
            "rerun.edr",
            "-o",
            "potential.xvg",
        ],
        cwd=output_dir,
        input="Potential\n",
    )
    if not potential.is_file():
        raise FileNotFoundError(f"Potential energy output not found: {potential}")
    for trajectory_output in (output_dir / "rerun.trr", output_dir / "rerun.xtc"):
        if trajectory_output.is_file():
            trajectory_output.unlink()
    return potential


def _calculate_bar(
    energy_paths,
    state_count,
    temperature,
    begin,
    end,
    subsample,
):
    from pymbar.other_estimators import bar

    beta = 1.0 / (GAS_CONSTANT_KJ * temperature)
    pair_results = []
    total = 0.0
    total_variance = 0.0
    for state in range(state_count - 1):
        forward = (
            _aligned_difference(
                energy_paths[(state, state)],
                energy_paths[(state, state + 1)],
                begin,
                end,
                subsample,
            )
            * beta
        )
        reverse = (
            _aligned_difference(
                energy_paths[(state + 1, state + 1)],
                energy_paths[(state + 1, state)],
                begin,
                end,
                subsample,
            )
            * beta
        )
        result = bar(forward, reverse)
        delta_g = float(result["Delta_f"]) / beta
        uncertainty = float(result["dDelta_f"]) / beta
        if not math.isfinite(uncertainty):
            if np.var(forward) < 1e-20 and np.var(reverse) < 1e-20:
                uncertainty = 0.0
            else:
                raise ValueError(
                    f"BAR uncertainty is not finite for states {state} and {state + 1}"
                )
        pair_results.append(
            {
                "state_a": state,
                "state_b": state + 1,
                "delta_g_kj_mol": delta_g,
                "uncertainty_kj_mol": uncertainty,
                "forward_samples": len(forward),
                "reverse_samples": len(reverse),
            }
        )
        total += delta_g
        total_variance += uncertainty**2
    return pair_results, total, math.sqrt(total_variance)


def _convergence_entry(index, begin, end, delta_g, uncertainty):
    return {
        "index": index,
        "begin_ps": begin,
        "end_ps": end,
        "delta_g_kj_mol": delta_g,
        "uncertainty_kj_mol": uncertainty,
        "delta_g_kcal_mol": delta_g * KJ_TO_KCAL,
        "uncertainty_kcal_mol": uncertainty * KJ_TO_KCAL,
    }


def run(args):
    if args.begin < 0:
        raise ValueError("--begin must be non-negative")
    if args.end is not None and args.end <= args.begin:
        raise ValueError("--end must be greater than --begin")
    if args.subsample <= 0:
        raise ValueError("--subsample must be positive")
    if args.ntomp <= 0:
        raise ValueError("--ntomp must be positive")
    convergence_blocks = getattr(args, "convergence_blocks", 0)
    if convergence_blocks == 1 or convergence_blocks < 0:
        raise ValueError("--convergence-blocks must be 0 or at least 2")

    base, manifest = load_fep_manifest(args.path)
    if manifest.get("workflow") != "fep-rest":
        raise ValueError("Manifest is not a FEP/REST calculation")
    temperature = float(manifest["physical_temperature"])
    args.gmx = args.gmx or manifest.get("gmx_executable", "gmx_mpi")
    deffnm = manifest["deffnm"]
    windows = manifest["windows"]
    trajectories = [
        (base / window["directory"] / f"{deffnm}.trr").resolve() for window in windows
    ]
    tprs = [
        (base / window["directory"] / f"{deffnm}.tpr").resolve() for window in windows
    ]
    for path in [*trajectories, *tprs]:
        if not path.is_file():
            raise FileNotFoundError(f"FEP/REST input not found: {path}")

    rerun_root = base / "rerun_energy"
    energy_paths = {}
    for simulation in range(len(windows)):
        evaluation_states = {
            simulation,
            max(0, simulation - 1),
            min(len(windows) - 1, simulation + 1),
        }
        for evaluation in sorted(evaluation_states):
            energy_paths[(simulation, evaluation)] = _rerun_energy(
                args,
                trajectories[simulation],
                tprs[evaluation],
                rerun_root / f"sim_{simulation:03d}_eval_{evaluation:03d}",
            )
    ranges = (
        convergence_ranges(
            energy_paths.values(),
            args.begin,
            args.end,
            convergence_blocks,
        )
        if convergence_blocks
        else None
    )
    analysis_begin = ranges["effective_begin_ps"] if ranges else args.begin
    analysis_end = ranges["effective_end_ps"] if ranges else args.end

    pair_results, delta_g, uncertainty = _calculate_bar(
        energy_paths,
        len(windows),
        temperature,
        analysis_begin,
        analysis_end,
        args.subsample,
    )
    output_path = (
        Path(args.output).expanduser() if args.output else base / "fep_rest_result.json"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output = {
        "method": "neighbor BAR from cross-Hamiltonian reruns",
        "workflow": "fep-rest",
        "temperature_k": temperature,
        "delta_g_kj_mol": delta_g,
        "uncertainty_kj_mol": uncertainty,
        "delta_g_kcal_mol": delta_g * KJ_TO_KCAL,
        "uncertainty_kcal_mol": uncertainty * KJ_TO_KCAL,
        "begin_ps": analysis_begin,
        "end_ps": analysis_end,
        "subsample": args.subsample,
        "pairs": pair_results,
    }
    if convergence_blocks:
        block_estimates = []
        cumulative_estimates = []
        for index, (range_begin, range_end) in enumerate(ranges["blocks"]):
            _pairs, value, error = _calculate_bar(
                energy_paths,
                len(windows),
                temperature,
                range_begin,
                range_end,
                args.subsample,
            )
            block_estimates.append(
                _convergence_entry(index, range_begin, range_end, value, error)
            )
        for index, (range_begin, range_end) in enumerate(ranges["cumulative"][:-1]):
            _pairs, value, error = _calculate_bar(
                energy_paths,
                len(windows),
                temperature,
                range_begin,
                range_end,
                args.subsample,
            )
            cumulative_estimates.append(
                _convergence_entry(index, range_begin, range_end, value, error)
            )
        final_begin, final_end = ranges["cumulative"][-1]
        cumulative_estimates.append(
            _convergence_entry(
                convergence_blocks - 1,
                final_begin,
                final_end,
                delta_g,
                uncertainty,
            )
        )
        output["convergence"] = {
            "effective_begin_ps": ranges["effective_begin_ps"],
            "effective_end_ps": ranges["effective_end_ps"],
            "block_estimates": block_estimates,
            "cumulative_estimates": cumulative_estimates,
        }
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    LOGGER.info(
        "FEP/REST Delta G = %.3f +/- %.3f kJ/mol",
        delta_g,
        uncertainty,
    )
