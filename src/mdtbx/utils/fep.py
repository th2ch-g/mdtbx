import json
import math
from pathlib import Path

FEP_MANIFEST = "fep_manifest.json"
FEP_SCHEMA_VERSION = 1
KJ_TO_KCAL = 1.0 / 4.184

FEP_MDP_KEYS = {
    "free-energy",
    "init-lambda-state",
    "fep-lambdas",
    "coul-lambdas",
    "vdw-lambdas",
    "bonded-lambdas",
    "restraint-lambdas",
    "mass-lambdas",
    "calc-lambda-neighbors",
    "dhdl-derivatives",
    "dhdl-print-energy",
    "separate-dhdl-file",
    "nstdhdl",
    "couple-moltype",
    "couple-lambda0",
    "couple-lambda1",
    "couple-intramol",
    "sc-alpha",
    "sc-power",
    "sc-r-power",
    "sc-sigma",
    "sc-coul",
}


def _normalize_mdp_key(key):
    return key.strip().lower().replace("_", "-")


def render_fep_mdp(base_text, settings, *, remove_prefixes=()):
    normalized_settings = {
        _normalize_mdp_key(key): (key, str(value)) for key, value in settings.items()
    }
    normalized_prefixes = {_normalize_mdp_key(prefix) for prefix in remove_prefixes}
    controlled_keys = {
        *(_normalize_mdp_key(key) for key in FEP_MDP_KEYS),
        *normalized_settings,
    }
    emitted = set()
    lines = []

    for line in base_text.splitlines(keepends=True):
        setting = line.partition(";")[0]
        key, separator, _value = setting.partition("=")
        if not separator:
            lines.append(line)
            continue

        normalized_key = _normalize_mdp_key(key)
        remove_by_prefix = any(
            normalized_key == prefix or normalized_key.startswith(f"{prefix}-")
            for prefix in normalized_prefixes
        )
        if normalized_key not in controlled_keys and not remove_by_prefix:
            lines.append(line)
            continue
        if normalized_key not in normalized_settings or normalized_key in emitted:
            continue

        output_key, output_value = normalized_settings[normalized_key]
        lines.append(f"{output_key:<24} = {output_value}\n")
        emitted.add(normalized_key)

    if lines and not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    if lines and lines[-1].strip():
        lines.append("\n")
    lines.append("; mdtbx FEP settings\n")
    for normalized_key, (output_key, output_value) in normalized_settings.items():
        if normalized_key not in emitted:
            lines.append(f"{output_key:<24} = {output_value}\n")
    return "".join(lines)


def temperature_mdp_overrides(base_text, temperature):
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("Temperature must be positive")
    overrides = {}
    for line in base_text.splitlines():
        setting = line.partition(";")[0]
        key, separator, value = setting.partition("=")
        if not separator:
            continue
        normalized_key = _normalize_mdp_key(key)
        if normalized_key == "ref-t":
            count = len(value.split())
            if count:
                overrides["ref-t"] = " ".join(f"{temperature:g}" for _ in range(count))
        elif normalized_key == "gen-temp":
            overrides["gen-temp"] = f"{temperature:g}"
    return overrides


def _linspace(count):
    if count < 2:
        raise ValueError("Lambda window count must be at least 2")
    return [index / (count - 1) for index in range(count)]


def _validate_lambdas(values, label):
    if len(values) < 2:
        raise ValueError(f"{label} must contain at least 2 values")
    if not all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in values):
        raise ValueError(f"{label} values must be finite and between 0 and 1")
    if any(current > following for current, following in zip(values, values[1:])):
        raise ValueError(f"{label} values must be non-decreasing")
    if values[0] != 0.0 or values[-1] != 1.0:
        raise ValueError(f"{label} must start at 0 and end at 1")
    return [float(value) for value in values]


def build_lambda_schedule(
    mode,
    *,
    windows=12,
    coul_windows=12,
    vdw_windows=12,
    fep_lambdas=None,
    coul_lambdas=None,
    vdw_lambdas=None,
):
    if mode == "transform":
        if coul_lambdas is not None or vdw_lambdas is not None:
            raise ValueError("transform mode only accepts --fep-lambdas")
        values = _validate_lambdas(
            fep_lambdas if fep_lambdas is not None else _linspace(windows),
            "fep lambdas",
        )
        return {"fep": values}

    if fep_lambdas is not None:
        raise ValueError(f"{mode} mode does not accept --fep-lambdas")

    if mode == "charge":
        if vdw_lambdas is not None:
            raise ValueError("charge mode does not accept --vdw-lambdas")
        coul = _validate_lambdas(
            coul_lambdas if coul_lambdas is not None else _linspace(windows),
            "coul lambdas",
        )
        return {"coul": coul, "vdw": [0.0] * len(coul)}

    if mode == "vdw":
        if coul_lambdas is not None:
            raise ValueError("vdw mode does not accept --coul-lambdas")
        vdw = _validate_lambdas(
            vdw_lambdas if vdw_lambdas is not None else _linspace(windows),
            "vdw lambdas",
        )
        return {"coul": [0.0] * len(vdw), "vdw": vdw}

    if mode != "decouple":
        raise ValueError(f"Unsupported FEP mode: {mode}")

    if coul_lambdas is not None or vdw_lambdas is not None:
        if coul_lambdas is None or vdw_lambdas is None:
            raise ValueError(
                "decouple mode requires both --coul-lambdas and --vdw-lambdas"
            )
        coul = _validate_lambdas(coul_lambdas, "coul lambdas")
        vdw = _validate_lambdas(vdw_lambdas, "vdw lambdas")
        if len(coul) != len(vdw):
            raise ValueError("coul and vdw lambda lists must have the same length")
        return {"coul": coul, "vdw": vdw}

    coul_stage = _linspace(coul_windows)
    vdw_stage = _linspace(vdw_windows)
    return {
        "coul": coul_stage + [1.0] * (vdw_windows - 1),
        "vdw": [0.0] * (coul_windows - 1) + vdw_stage,
    }


def format_lambdas(values):
    return " ".join(f"{value:.6f}" for value in values)


def load_fep_manifest(path):
    source = Path(path).expanduser()
    manifest_path = source if source.is_file() else source / FEP_MANIFEST
    if not manifest_path.is_file():
        raise FileNotFoundError(f"FEP manifest not found: {manifest_path}")

    data = json.loads(manifest_path.read_text())
    if data.get("schema_version") != FEP_SCHEMA_VERSION:
        raise ValueError("Unsupported FEP manifest schema")
    windows = data.get("windows")
    if not isinstance(windows, list) or len(windows) < 2:
        raise ValueError("FEP manifest must contain at least 2 windows")
    deffnm = data.get("deffnm")
    if not isinstance(deffnm, str) or not deffnm or Path(deffnm).name != deffnm:
        raise ValueError("FEP manifest contains an invalid deffnm")
    for index, window in enumerate(windows):
        if (
            not isinstance(window, dict)
            or window.get("index") != index
            or not isinstance(window.get("directory"), str)
            or not window["directory"]
            or Path(window["directory"]).name != window["directory"]
        ):
            raise ValueError(
                f"FEP manifest contains an invalid window at index {index}"
            )
    return manifest_path.parent, data
