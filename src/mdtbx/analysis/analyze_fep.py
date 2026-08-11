import argparse
import json
import tempfile
from pathlib import Path

from ..logger import generate_logger
from ..utils.fep import KJ_TO_KCAL, load_fep_manifest
from ..utils.convergence import convergence_ranges
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "analyze_fep",
        help="Analyze GROMACS FEP windows with BAR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--path",
        default="fep",
        type=str,
        help="FEP setup directory or manifest",
    )
    parser.add_argument(
        "-b",
        "--begin",
        default=0.0,
        type=float,
        help="Beginning time [ps]",
    )
    parser.add_argument(
        "-e",
        "--end",
        type=float,
        help="Ending time [ps]",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        help="Override temperature [K]",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        help="Output prefix",
    )
    parser.add_argument(
        "--gmx",
        type=str,
        help="GROMACS executable; defaults to the setup executable",
    )
    parser.add_argument(
        "--convergence-blocks",
        type=int,
        default=0,
        help="Number of equal-duration convergence blocks; 0 disables analysis",
    )
    parser.set_defaults(func=run)


def _parse_total(stdout):
    for line in reversed(stdout.splitlines()):
        fields = line.split()
        if not fields or fields[0].lower() != "total":
            continue
        for separator in ("+/-", "+-", "±"):
            if separator not in fields:
                continue
            index = fields.index(separator)
            if index == 0 or index + 1 >= len(fields):
                break
            return float(fields[index - 1]), float(fields[index + 1])
    raise ValueError("Could not parse the total BAR estimate")


def _prefixed_path(prefix, suffix):
    return prefix.parent / f"{prefix.name}{suffix}"


def _bar_estimate(gmx, dhdl_files, prefix, begin, end, temperature):
    command = [
        gmx,
        "bar",
        "-f",
        *[str(path) for path in dhdl_files],
        "-o",
        str(_prefixed_path(prefix, ".xvg")),
        "-oi",
        str(_prefixed_path(prefix, "_integral.xvg")),
        "-oh",
        str(_prefixed_path(prefix, "_histogram.xvg")),
        "-b",
        str(begin),
    ]
    if end is not None:
        command.extend(["-e", str(end)])
    if temperature is not None:
        command.extend(["-temp", str(temperature)])
    result = run_cmd(command, capture_output=True, text=True)
    log_text = result.stdout
    if result.stderr:
        log_text += f"\n{result.stderr}"
    delta_g, uncertainty = _parse_total(log_text)
    return delta_g, uncertainty, log_text


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
    if args.temperature is not None and args.temperature <= 0:
        raise ValueError("--temperature must be positive")
    convergence_blocks = getattr(args, "convergence_blocks", 0)
    if convergence_blocks == 1 or convergence_blocks < 0:
        raise ValueError("--convergence-blocks must be 0 or at least 2")

    base, manifest = load_fep_manifest(args.path)
    deffnm = manifest["deffnm"]
    gmx = args.gmx or manifest.get("gmx_executable", "gmx")
    dhdl_files = [
        (base / window["directory"] / f"{deffnm}.xvg").resolve()
        for window in manifest["windows"]
    ]
    missing = [path for path in dhdl_files if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Delta H file: {missing[0]}")
    ranges = (
        convergence_ranges(
            dhdl_files,
            args.begin,
            args.end,
            convergence_blocks,
        )
        if convergence_blocks
        else None
    )
    analysis_begin = ranges["effective_begin_ps"] if ranges else args.begin
    analysis_end = ranges["effective_end_ps"] if ranges else args.end

    prefix = (
        Path(args.output_prefix).expanduser() if args.output_prefix else base / "bar"
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    log_path = _prefixed_path(prefix, ".log")
    json_path = _prefixed_path(prefix, ".json")

    delta_g, uncertainty, log_text = _bar_estimate(
        gmx,
        dhdl_files,
        prefix,
        analysis_begin,
        analysis_end,
        args.temperature,
    )
    log_path.write_text(log_text)
    output = {
        "method": "BAR",
        "mode": manifest["mode"],
        "delta_g_kj_mol": delta_g,
        "uncertainty_kj_mol": uncertainty,
        "delta_g_kcal_mol": delta_g * KJ_TO_KCAL,
        "uncertainty_kcal_mol": uncertainty * KJ_TO_KCAL,
        "begin_ps": analysis_begin,
        "end_ps": analysis_end,
        "temperature_k": args.temperature,
        "windows": len(dhdl_files),
    }
    if convergence_blocks:
        block_estimates = []
        cumulative_estimates = []
        with tempfile.TemporaryDirectory(prefix="mdtbx-bar-") as temporary:
            temporary_root = Path(temporary)
            for index, (range_begin, range_end) in enumerate(ranges["blocks"]):
                value, error, _log = _bar_estimate(
                    gmx,
                    dhdl_files,
                    temporary_root / f"block_{index:03d}",
                    range_begin,
                    range_end,
                    args.temperature,
                )
                block_estimates.append(
                    _convergence_entry(index, range_begin, range_end, value, error)
                )
            for index, (range_begin, range_end) in enumerate(ranges["cumulative"][:-1]):
                value, error, _log = _bar_estimate(
                    gmx,
                    dhdl_files,
                    temporary_root / f"cumulative_{index:03d}",
                    range_begin,
                    range_end,
                    args.temperature,
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
    json_path.write_text(json.dumps(output, indent=2) + "\n")
    LOGGER.info(
        "Delta G = %.3f +/- %.3f kJ/mol (%.3f +/- %.3f kcal/mol)",
        delta_g,
        uncertainty,
        output["delta_g_kcal_mol"],
        output["uncertainty_kcal_mol"],
    )
