import json

import pytest

from src.utils.fep import (
    FEP_MANIFEST,
    build_lambda_schedule,
    load_fep_manifest,
    render_fep_mdp,
    temperature_mdp_overrides,
)


def test_decouple_schedule_stages_charge_before_vdw():
    schedule = build_lambda_schedule("decouple", coul_windows=3, vdw_windows=4)

    assert schedule["coul"] == [0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
    assert schedule["vdw"] == [0.0, 0.0, 0.0, 1 / 3, 2 / 3, 1.0]


@pytest.mark.parametrize(
    ("mode", "kwargs", "expected"),
    [
        ("charge", {"windows": 3}, {"coul": [0.0, 0.5, 1.0], "vdw": [0.0] * 3}),
        ("vdw", {"windows": 2}, {"coul": [0.0, 0.0], "vdw": [0.0, 1.0]}),
        ("transform", {"windows": 3}, {"fep": [0.0, 0.5, 1.0]}),
    ],
)
def test_component_schedules(mode, kwargs, expected):
    assert build_lambda_schedule(mode, **kwargs) == expected


def test_custom_decouple_schedule_requires_matching_components():
    with pytest.raises(ValueError, match="requires both"):
        build_lambda_schedule(
            "decouple",
            coul_lambdas=[0.0, 1.0],
        )


@pytest.mark.parametrize(
    "values",
    [
        [0.1, 1.0],
        [0.0, 0.8],
        [0.0, 0.7, 0.6, 1.0],
        [0.0, float("nan"), 1.0],
    ],
)
def test_lambda_validation_rejects_invalid_schedules(values):
    with pytest.raises(ValueError):
        build_lambda_schedule("transform", fep_lambdas=values)


def test_render_fep_mdp_replaces_and_removes_controlled_settings():
    base = (
        "integrator = md\n"
        "free_energy = no ; stale\n"
        "fep-lambdas = 0 1\n"
        "couple_moltype = OLD\n"
        "nsteps = 100\n"
    )

    rendered = render_fep_mdp(
        base,
        {
            "free-energy": "yes",
            "init-lambda-state": 1,
            "fep-lambdas": "0.000000 0.500000 1.000000",
        },
    )

    assert "integrator = md" in rendered
    assert "nsteps = 100" in rendered
    assert "free-energy" in rendered
    assert "free_energy" not in rendered
    assert "couple_moltype" not in rendered
    assert rendered.count("fep-lambdas") == 1
    assert rendered.count("init-lambda-state") == 1


def test_render_fep_mdp_removes_matching_prefixes():
    rendered = render_fep_mdp(
        "pull = yes\npull-coord7-k = 10\nnsteps = 100\n",
        {"pull": "no"},
        remove_prefixes=("pull",),
    )

    assert "pull-coord7-k" not in rendered
    assert "pull" in rendered
    assert "nsteps = 100" in rendered


def test_temperature_overrides_preserve_coupling_group_count():
    overrides = temperature_mdp_overrides(
        "ref_t = 298 298\ngen-temp = 298\n",
        310.0,
    )

    assert overrides == {"ref-t": "310 310", "gen-temp": "310"}


def test_load_fep_manifest_accepts_directory_and_file(tmp_path):
    manifest = {
        "schema_version": 1,
        "deffnm": "fep",
        "windows": [
            {"index": 0, "directory": "lambda_000"},
            {"index": 1, "directory": "lambda_001"},
        ],
    }
    path = tmp_path / FEP_MANIFEST
    path.write_text(json.dumps(manifest))

    assert load_fep_manifest(tmp_path) == (tmp_path, manifest)
    assert load_fep_manifest(path) == (tmp_path, manifest)


def test_load_fep_manifest_rejects_unsafe_window_path(tmp_path):
    manifest = {
        "schema_version": 1,
        "deffnm": "fep",
        "windows": [
            {"index": 0, "directory": "../lambda_000"},
            {"index": 1, "directory": "lambda_001"},
        ],
    }
    (tmp_path / FEP_MANIFEST).write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="invalid window"):
        load_fep_manifest(tmp_path)
