"""Combine two alchemical legs into a relative binding free energy."""

import argparse
import json
import math
from pathlib import Path

from ..logger import generate_logger

LOGGER = generate_logger(__name__)
GAS_CONSTANT_KJ = 0.00831446261815324


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "analyze_rbfe",
        help="Combine complex and solvent free energies and compare with Ki data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--complex", required=True, help="Complex-leg result JSON")
    parser.add_argument("--solvent", required=True, help="Solvent-leg result JSON")
    parser.add_argument("--ki-a", required=True, type=float)
    parser.add_argument("--ki-a-error", required=True, type=float)
    parser.add_argument("--ki-b", required=True, type=float)
    parser.add_argument("--ki-b-error", required=True, type=float)
    parser.add_argument("--ki-unit", default="uM")
    parser.add_argument("--temperature", default=300.0, type=float)
    parser.add_argument("--ligand-a", required=True)
    parser.add_argument("--ligand-b", required=True)
    parser.add_argument("-o", "--output", default="rbfe_result.json")
    parser.set_defaults(func=run)


def _result(path):
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    data = json.loads(resolved.read_text())
    try:
        delta_g = float(data["delta_g_kj_mol"])
        uncertainty = float(data["uncertainty_kj_mol"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid free-energy result: {resolved}") from error
    if not all(math.isfinite(value) for value in (delta_g, uncertainty)):
        raise ValueError(f"Non-finite free-energy result: {resolved}")
    if uncertainty < 0:
        raise ValueError(f"Negative uncertainty in {resolved}")
    return resolved, data, delta_g, uncertainty


def calculate(
    complex_delta_g,
    complex_uncertainty,
    solvent_delta_g,
    solvent_uncertainty,
    ki_a,
    ki_a_error,
    ki_b,
    ki_b_error,
    temperature,
):
    if temperature <= 0 or not math.isfinite(temperature):
        raise ValueError("Temperature must be positive")
    if ki_a <= 0 or ki_b <= 0:
        raise ValueError("Ki values must be positive")
    if ki_a_error < 0 or ki_b_error < 0:
        raise ValueError("Ki errors must be non-negative")
    delta_delta_g = complex_delta_g - solvent_delta_g
    calculated_uncertainty = math.hypot(complex_uncertainty, solvent_uncertainty)
    rt = GAS_CONSTANT_KJ * temperature
    experimental = rt * math.log(ki_b / ki_a)
    experimental_uncertainty = rt * math.hypot(ki_a_error / ki_a, ki_b_error / ki_b)
    residual = delta_delta_g - experimental
    combined_uncertainty = math.hypot(calculated_uncertainty, experimental_uncertainty)
    z_score = residual / combined_uncertainty if combined_uncertainty else None
    return {
        "calculated_delta_delta_g_kj_mol": delta_delta_g,
        "calculated_uncertainty_kj_mol": calculated_uncertainty,
        "experimental_delta_delta_g_kj_mol": experimental,
        "experimental_uncertainty_kj_mol": experimental_uncertainty,
        "residual_kj_mol": residual,
        "combined_uncertainty_kj_mol": combined_uncertainty,
        "z_score": z_score,
    }


def run(args):
    complex_path, complex_data, complex_delta_g, complex_uncertainty = _result(
        args.complex
    )
    solvent_path, solvent_data, solvent_delta_g, solvent_uncertainty = _result(
        args.solvent
    )
    values = calculate(
        complex_delta_g,
        complex_uncertainty,
        solvent_delta_g,
        solvent_uncertainty,
        args.ki_a,
        args.ki_a_error,
        args.ki_b,
        args.ki_b_error,
        args.temperature,
    )
    output = {
        "schema_version": 1,
        "workflow": "rbfe",
        "direction": f"{args.ligand_a}->{args.ligand_b}",
        "definition": "DeltaG_complex(A->B) - DeltaG_solvent(A->B)",
        "temperature_k": args.temperature,
        "complex": {
            "path": str(complex_path),
            "delta_g_kj_mol": complex_delta_g,
            "uncertainty_kj_mol": complex_uncertainty,
            "begin_ps": complex_data.get("begin_ps"),
            "end_ps": complex_data.get("end_ps"),
        },
        "solvent": {
            "path": str(solvent_path),
            "delta_g_kj_mol": solvent_delta_g,
            "uncertainty_kj_mol": solvent_uncertainty,
            "begin_ps": solvent_data.get("begin_ps"),
            "end_ps": solvent_data.get("end_ps"),
        },
        "experiment": {
            "ki_a": args.ki_a,
            "ki_a_error": args.ki_a_error,
            "ki_b": args.ki_b,
            "ki_b_error": args.ki_b_error,
            "ki_unit": args.ki_unit,
        },
        **values,
    }
    output_path = Path(args.output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    tsv_path = output_path.with_suffix(".tsv")
    columns = [
        "direction",
        "calculated_delta_delta_g_kj_mol",
        "calculated_uncertainty_kj_mol",
        "experimental_delta_delta_g_kj_mol",
        "experimental_uncertainty_kj_mol",
        "residual_kj_mol",
        "z_score",
    ]
    tsv_path.write_text(
        "\t".join(columns)
        + "\n"
        + "\t".join(str(output[column]) for column in columns)
        + "\n"
    )
    LOGGER.info(
        "RBFE Delta Delta G = %.3f +/- %.3f kJ/mol",
        values["calculated_delta_delta_g_kj_mol"],
        values["calculated_uncertainty_kj_mol"],
    )
