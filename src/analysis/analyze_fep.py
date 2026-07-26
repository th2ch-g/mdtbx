import argparse
import json
from pathlib import Path

from ..logger import generate_logger
from ..utils.fep import KJ_TO_KCAL, load_fep_manifest
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


def run(args):
    if args.begin < 0:
        raise ValueError("--begin must be non-negative")
    if args.end is not None and args.end <= args.begin:
        raise ValueError("--end must be greater than --begin")
    if args.temperature is not None and args.temperature <= 0:
        raise ValueError("--temperature must be positive")

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

    prefix = (
        Path(args.output_prefix).expanduser() if args.output_prefix else base / "bar"
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    bar_path = _prefixed_path(prefix, ".xvg")
    integral_path = _prefixed_path(prefix, "_integral.xvg")
    histogram_path = _prefixed_path(prefix, "_histogram.xvg")
    log_path = _prefixed_path(prefix, ".log")
    json_path = _prefixed_path(prefix, ".json")

    command = [
        gmx,
        "bar",
        "-f",
        *[str(path) for path in dhdl_files],
        "-o",
        str(bar_path),
        "-oi",
        str(integral_path),
        "-oh",
        str(histogram_path),
        "-b",
        str(args.begin),
    ]
    if args.end is not None:
        command.extend(["-e", str(args.end)])
    if args.temperature is not None:
        command.extend(["-temp", str(args.temperature)])

    result = run_cmd(command, capture_output=True, text=True)
    log_text = result.stdout
    if result.stderr:
        log_text += f"\n{result.stderr}"
    log_path.write_text(log_text)
    delta_g, uncertainty = _parse_total(log_text)
    output = {
        "method": "BAR",
        "mode": manifest["mode"],
        "delta_g_kj_mol": delta_g,
        "uncertainty_kj_mol": uncertainty,
        "delta_g_kcal_mol": delta_g * KJ_TO_KCAL,
        "uncertainty_kcal_mol": uncertainty * KJ_TO_KCAL,
        "begin_ps": args.begin,
        "end_ps": args.end,
        "temperature_k": args.temperature,
        "windows": len(dhdl_files),
    }
    json_path.write_text(json.dumps(output, indent=2) + "\n")
    LOGGER.info(
        "Delta G = %.3f +/- %.3f kJ/mol (%.3f +/- %.3f kcal/mol)",
        delta_g,
        uncertainty,
        output["delta_g_kcal_mol"],
        output["uncertainty_kcal_mol"],
    )
