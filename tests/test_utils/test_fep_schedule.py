import json

import pytest

from mdtbx.utils.fep_schedule import (
    components_from_coordinates,
    load_optimized_schedule,
    make_optimized_schedule,
    optimize_coordinates,
    parse_exchange_probabilities,
    progress_coordinates,
)


def test_parse_exchange_probabilities_uses_last_complete_table():
    log = """\
Repl  average probabilities:
Repl     0    1    2
Repl    .10  .20
Repl  average probabilities:
Repl     0    1    2
Repl    .30
"""

    assert parse_exchange_probabilities(log, 3) == [0.1, 0.2]


def test_optimize_coordinates_preserves_phase_boundary_and_damps_update():
    optimized = optimize_coordinates(
        [0.0, 0.5, 1.0, 1.5, 2.0],
        [0.01, 0.8, 0.8, 0.8],
        boundaries=[1.0],
        iteration=1,
    )

    assert optimized[0] == 0.0
    assert optimized[2] == 1.0
    assert optimized[-1] == 2.0
    assert optimized[1] < 0.5
    assert optimized[3] == pytest.approx(1.5)


def test_progress_coordinates_rejects_overlapping_decouple_stages():
    manifest = {
        "mode": "decouple",
        "lambda_components": {
            "coul": [0.0, 0.5, 1.0],
            "vdw": [0.0, 0.2, 1.0],
        },
    }

    with pytest.raises(ValueError, match="staged"):
        progress_coordinates(manifest)


def test_optimized_schedule_round_trip(tmp_path):
    manifest = {
        "mode": "transform",
        "lambda_components": {"fep": [0.0, 0.5, 1.0]},
    }
    data = make_optimized_schedule(
        manifest,
        [0.2, 0.8],
        source_manifest=tmp_path / "fep_manifest.json",
        source_log=tmp_path / "fep.log",
    )
    path = tmp_path / "schedule.json"
    path.write_text(json.dumps(data))

    resolved, loaded = load_optimized_schedule(
        path,
        expected_mode="transform",
        expected_workflow="fep",
    )

    assert resolved == path.resolve()
    assert loaded["state_count"] == 3
    with pytest.raises(ValueError, match="does not match"):
        load_optimized_schedule(path, expected_mode="vdw")


def test_fep_rest_components_keep_three_phases():
    components = components_from_coordinates(
        "transform",
        "fep-rest",
        [0.0, 1.0, 1.5, 2.0, 3.0],
    )

    assert components["charge_a"] == [0.0, 1.0, 1.0, 1.0, 1.0]
    assert components["vdw"] == [0.0, 0.0, 0.5, 1.0, 1.0]
    assert components["charge_b"] == [0.0, 0.0, 0.0, 0.0, 1.0]
