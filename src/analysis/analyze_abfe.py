import argparse
import json
import math
from pathlib import Path
from types import SimpleNamespace

from ..logger import generate_logger
from ..utils.abfe import (
    GAS_CONSTANT_KJ,
    boresch_standard_state_correction,
    load_abfe_manifest,
)
from ..utils.fep import KJ_TO_KCAL
from . import analyze_fep

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "analyze_abfe",
        help="Analyze and combine an ABFE thermodynamic cycle",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--path", default="abfe", help="ABFE setup directory")
    parser.add_argument(
        "-b", "--begin", type=float, default=0.0, help="Begin time [ps]"
    )
    parser.add_argument("-e", "--end", type=float, help="End time [ps]")
    parser.add_argument("--temperature", type=float, help="Override temperature [K]")
    parser.add_argument(
        "--symmetry-number",
        type=int,
        default=1,
        help="Ligand rotational symmetry number",
    )
    parser.add_argument(
        "--correction",
        nargs=3,
        action="append",
        metavar=("NAME", "VALUE_KJ_MOL", "ERROR_KJ_MOL"),
        help="Additional signed ABFE correction; repeat as needed",
    )
    parser.add_argument("-o", "--output", help="Output JSON file")
    parser.add_argument(
        "--gmx",
        help="GROMACS executable; defaults to the setup executable",
    )
    parser.set_defaults(func=run)


def _analyze_leg(args, base, relative_path, temperature):
    leg_path = base / relative_path
    analyze_fep.run(
        SimpleNamespace(
            path=str(leg_path),
            begin=args.begin,
            end=args.end,
            temperature=temperature,
            output_prefix=str(leg_path / "bar"),
            gmx=args.gmx,
        )
    )
    return json.loads((leg_path / "bar.json").read_text())


def run(args):
    if args.symmetry_number <= 0:
        raise ValueError("--symmetry-number must be positive")
    base, manifest = load_abfe_manifest(args.path)
    temperature = (
        args.temperature if args.temperature is not None else manifest["temperature"]
    )
    if temperature <= 0:
        raise ValueError("Temperature must be positive")

    signs = {
        "solvent_charge": 1.0,
        "solvent_vdw": 1.0,
        "complex_charge": -1.0,
        "complex_vdw": -1.0,
        "complex_restraint": 1.0,
    }
    contributions = []
    total = 0.0
    variance = 0.0
    for leg, sign in signs.items():
        result = _analyze_leg(
            args,
            base,
            manifest["legs"][leg],
            temperature,
        )
        signed_value = sign * result["delta_g_kj_mol"]
        uncertainty = result["uncertainty_kj_mol"]
        contributions.append(
            {
                "name": leg,
                "coefficient": sign,
                "value_kj_mol": signed_value,
                "uncertainty_kj_mol": uncertainty,
            }
        )
        total += signed_value
        variance += uncertainty**2

    springs = manifest["springs"]
    standard_state = boresch_standard_state_correction(
        manifest["geometry"],
        temperature=temperature,
        distance_spring=springs["distance"],
        angle_spring=springs["angle"],
        dihedral_spring=springs["dihedral"],
        symmetry_number=args.symmetry_number,
    )
    contributions.append(
        {
            "name": "boresch_standard_state",
            "coefficient": 1.0,
            "value_kj_mol": standard_state,
            "uncertainty_kj_mol": 0.0,
        }
    )
    total += standard_state

    for correction in args.correction or []:
        name, value_text, error_text = correction
        value = float(value_text)
        uncertainty = float(error_text)
        if not math.isfinite(value) or not math.isfinite(uncertainty):
            raise ValueError("Correction values must be finite")
        if uncertainty < 0:
            raise ValueError("Correction uncertainty must be non-negative")
        contributions.append(
            {
                "name": name,
                "coefficient": 1.0,
                "value_kj_mol": value,
                "uncertainty_kj_mol": uncertainty,
            }
        )
        total += value
        variance += uncertainty**2

    uncertainty = math.sqrt(variance)
    exponent = total / (GAS_CONSTANT_KJ * temperature)
    if exponent > 709:
        kd_molar = None
    elif exponent < -745:
        kd_molar = 0.0
    else:
        kd_molar = math.exp(exponent)
    output = {
        "method": "Boresch-restraint double decoupling",
        "temperature_k": temperature,
        "delta_g_binding_kj_mol": total,
        "uncertainty_kj_mol": uncertainty,
        "delta_g_binding_kcal_mol": total * KJ_TO_KCAL,
        "uncertainty_kcal_mol": uncertainty * KJ_TO_KCAL,
        "kd_molar": kd_molar,
        "log10_kd_molar": exponent / math.log(10.0),
        "symmetry_number": args.symmetry_number,
        "long_range_method": manifest["long_range_method"],
        "contributions": contributions,
    }
    output_path = (
        Path(args.output).expanduser() if args.output else base / "abfe_result.json"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    text_path = output_path.with_suffix(".txt")
    lines = [
        f"{item['name']} {item['value_kj_mol']:.6f} {item['uncertainty_kj_mol']:.6f}"
        for item in contributions
    ]
    lines.append("----")
    lines.append(f"total {total:.6f} {uncertainty:.6f}")
    text_path.write_text("\n".join(lines) + "\n")
    LOGGER.info(
        "ABFE Delta G = %.3f +/- %.3f kJ/mol",
        total,
        uncertainty,
    )
