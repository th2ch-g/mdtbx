import bisect
import json
import math
import re
from pathlib import Path


SCHEDULE_SCHEMA_VERSION = 1
SCHEDULE_WORKFLOW = "optimized-fep-schedule"
_FLOAT_TOKEN = re.compile(r"(?<!\S)[+-]?(?:\d+\.\d*|\.\d+|\d+(?:[eE][+-]?\d+))(?!\S)")


def parse_exchange_probabilities(text, state_count):
    if state_count < 2:
        raise ValueError("At least 2 states are required")
    lines = text.splitlines()
    markers = [
        index
        for index, line in enumerate(lines)
        if "average probabilities" in line.lower()
    ]
    expected = state_count - 1
    for marker_index, marker in reversed(list(enumerate(markers))):
        stop = (
            markers[marker_index + 1] if marker_index + 1 < len(markers) else len(lines)
        )
        candidates = []
        for line in lines[marker + 1 : stop]:
            if not line.lstrip().lower().startswith("repl"):
                continue
            values = [float(value) for value in _FLOAT_TOKEN.findall(line)]
            if len(values) == expected and all(0.0 <= value <= 1.0 for value in values):
                candidates.append(values)
        if candidates:
            return candidates[-1]
    raise ValueError(
        f"Could not find a complete exchange-probability table for {state_count} states"
    )


def _validate_coordinates(coordinates, *, maximum):
    values = [float(value) for value in coordinates]
    if len(values) < 2:
        raise ValueError("Lambda schedule must contain at least 2 states")
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Lambda coordinates must be finite")
    if not math.isclose(values[0], 0.0) or not math.isclose(values[-1], maximum):
        raise ValueError(f"Lambda coordinates must start at 0 and end at {maximum:g}")
    if any(current >= following for current, following in zip(values, values[1:])):
        raise ValueError("Lambda coordinates must be strictly increasing")
    return values


def progress_coordinates(manifest):
    mode = manifest.get("mode")
    components = manifest.get("lambda_components")
    if not isinstance(components, dict):
        raise ValueError("FEP manifest does not contain lambda components")
    source_workflow = manifest.get("workflow", "fep")
    if source_workflow == "fep-rest":
        coordinates = _validate_coordinates(
            components.get("coordinates", []), maximum=3.0
        )
        for boundary in (1.0, 2.0):
            if not any(math.isclose(value, boundary) for value in coordinates):
                raise ValueError(
                    "FEP/REST schedule must contain the phase boundaries 1 and 2"
                )
        return coordinates, [1.0, 2.0]

    if mode == "transform":
        return _validate_coordinates(components.get("fep", []), maximum=1.0), []
    if mode == "charge":
        return _validate_coordinates(components.get("coul", []), maximum=1.0), []
    if mode == "vdw":
        return _validate_coordinates(components.get("vdw", []), maximum=1.0), []
    if mode != "decouple":
        raise ValueError(f"Unsupported FEP mode: {mode}")

    coul = [float(value) for value in components.get("coul", [])]
    vdw = [float(value) for value in components.get("vdw", [])]
    if len(coul) != len(vdw) or len(coul) < 2:
        raise ValueError(
            "Decouple schedule must contain matching coul and vdw components"
        )
    coordinates = []
    for coul_value, vdw_value in zip(coul, vdw):
        if not all(
            math.isfinite(value) and 0.0 <= value <= 1.0
            for value in (coul_value, vdw_value)
        ):
            raise ValueError(
                "Decouple lambda values must be finite and between 0 and 1"
            )
        if not math.isclose(vdw_value, 0.0) and not math.isclose(coul_value, 1.0):
            raise ValueError(
                "Only staged charge-then-vdw decouple schedules can be optimized"
            )
        coordinates.append(
            coul_value if math.isclose(vdw_value, 0.0) else 1.0 + vdw_value
        )
    coordinates = _validate_coordinates(coordinates, maximum=2.0)
    if not any(math.isclose(value, 1.0) for value in coordinates):
        raise ValueError("Decouple schedule must contain the charge/vdw phase boundary")
    return coordinates, [1.0]


def components_from_coordinates(mode, source_workflow, coordinates):
    values = [float(value) for value in coordinates]
    if source_workflow == "fep-rest":
        return {
            "coordinates": values,
            "fep": [value / 3.0 for value in values],
            "vdw": [min(max(value - 1.0, 0.0), 1.0) for value in values],
            "charge_a": [min(value, 1.0) for value in values],
            "charge_b": [max(value - 2.0, 0.0) for value in values],
        }
    if mode == "transform":
        return {"fep": values}
    if mode == "charge":
        return {"coul": values, "vdw": [0.0] * len(values)}
    if mode == "vdw":
        return {"coul": [0.0] * len(values), "vdw": values}
    if mode == "decouple":
        return {
            "coul": [min(value, 1.0) for value in values],
            "vdw": [max(value - 1.0, 0.0) for value in values],
        }
    raise ValueError(f"Unsupported FEP mode: {mode}")


def _redistribute_segment(coordinates, probabilities, start, stop):
    edge_count = stop - start
    difficulties = [
        max(-math.log(probabilities[index]), 1.0e-12) for index in range(start, stop)
    ]
    cumulative = [0.0]
    for difficulty in difficulties:
        cumulative.append(cumulative[-1] + difficulty)
    proposal = [coordinates[start]]
    for offset in range(1, edge_count):
        target = cumulative[-1] * offset / edge_count
        edge = min(bisect.bisect_right(cumulative, target) - 1, edge_count - 1)
        low = cumulative[edge]
        high = cumulative[edge + 1]
        fraction = (target - low) / (high - low)
        left = coordinates[start + edge]
        right = coordinates[start + edge + 1]
        proposal.append(left + fraction * (right - left))
    proposal.append(coordinates[stop])
    return proposal


def optimize_coordinates(
    coordinates,
    probabilities,
    *,
    boundaries=(),
    iteration=1,
    minimum_probability=0.01,
):
    if iteration < 1:
        raise ValueError("Iteration must be at least 1")
    if not math.isfinite(minimum_probability) or not 0.0 < minimum_probability <= 1.0:
        raise ValueError("Minimum probability must be finite and in (0, 1]")
    if len(probabilities) != len(coordinates) - 1:
        raise ValueError("Exchange probability count does not match the state count")
    probabilities = [float(value) for value in probabilities]
    if not all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in probabilities):
        raise ValueError("Exchange probabilities must be finite and between 0 and 1")
    clipped = [max(value, minimum_probability) for value in probabilities]
    boundary_indices = [0]
    for boundary in boundaries:
        matches = [
            index
            for index, value in enumerate(coordinates)
            if math.isclose(value, boundary)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Lambda schedule must contain boundary {boundary:g} exactly once"
            )
        boundary_indices.append(matches[0])
    boundary_indices.append(len(coordinates) - 1)

    proposal = list(coordinates)
    for start, stop in zip(boundary_indices, boundary_indices[1:]):
        proposal[start : stop + 1] = _redistribute_segment(
            coordinates,
            clipped,
            start,
            stop,
        )
    weight = 1.0 / (iteration + 1.0)
    optimized = [old + weight * (new - old) for old, new in zip(coordinates, proposal)]
    for index in boundary_indices:
        optimized[index] = coordinates[index]
    return optimized


def make_optimized_schedule(
    manifest,
    probabilities,
    *,
    source_manifest,
    source_log,
    iteration=1,
    minimum_probability=0.01,
):
    coordinates, boundaries = progress_coordinates(manifest)
    optimized = optimize_coordinates(
        coordinates,
        probabilities,
        boundaries=boundaries,
        iteration=iteration,
        minimum_probability=minimum_probability,
    )
    source_workflow = manifest.get("workflow", "fep")
    mode = manifest.get("mode")
    return {
        "schema_version": SCHEDULE_SCHEMA_VERSION,
        "workflow": SCHEDULE_WORKFLOW,
        "source_manifest": str(Path(source_manifest).resolve()),
        "source_log": str(Path(source_log).resolve()),
        "source_workflow": source_workflow,
        "mode": mode,
        "state_count": len(coordinates),
        "iteration": iteration,
        "minimum_probability": minimum_probability,
        "exchange_probabilities": [float(value) for value in probabilities],
        "coordinates": optimized,
        "lambda_components": components_from_coordinates(
            mode,
            source_workflow,
            optimized,
        ),
    }


def load_optimized_schedule(path, *, expected_mode=None, expected_workflow=None):
    schedule_path = Path(path).expanduser().resolve()
    if not schedule_path.is_file():
        raise FileNotFoundError(f"FEP schedule not found: {schedule_path}")
    data = json.loads(schedule_path.read_text())
    if data.get("schema_version") != SCHEDULE_SCHEMA_VERSION:
        raise ValueError("Unsupported FEP schedule schema")
    if data.get("workflow") != SCHEDULE_WORKFLOW:
        raise ValueError("File is not an optimized FEP schedule")
    if expected_mode is not None and data.get("mode") != expected_mode:
        raise ValueError(
            f"FEP schedule mode {data.get('mode')!r} does not match {expected_mode!r}"
        )
    if (
        expected_workflow is not None
        and data.get("source_workflow") != expected_workflow
    ):
        raise ValueError(
            f"FEP schedule source workflow does not match {expected_workflow!r}"
        )
    state_count = data.get("state_count")
    components = data.get("lambda_components")
    if (
        not isinstance(state_count, int)
        or state_count < 2
        or not isinstance(components, dict)
    ):
        raise ValueError("FEP schedule contains invalid state metadata")
    for name, values in components.items():
        if not isinstance(values, list) or len(values) != state_count:
            raise ValueError(f"FEP schedule component {name!r} has an invalid length")
        if not all(
            isinstance(value, (int, float)) and math.isfinite(value) for value in values
        ):
            raise ValueError(f"FEP schedule component {name!r} contains invalid values")
    coordinates, _boundaries = progress_coordinates(
        {
            "workflow": data.get("source_workflow"),
            "mode": data.get("mode"),
            "lambda_components": components,
        }
    )
    stored_coordinates = data.get("coordinates")
    if (
        not isinstance(stored_coordinates, list)
        or len(stored_coordinates) != state_count
        or not all(isinstance(value, (int, float)) for value in stored_coordinates)
        or any(
            not math.isclose(float(stored), calculated)
            for stored, calculated in zip(stored_coordinates, coordinates)
        )
    ):
        raise ValueError("FEP schedule coordinates do not match its lambda components")
    return schedule_path, data
