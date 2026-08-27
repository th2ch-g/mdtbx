import types

import pytest


def _args(**overrides):
    defaults = {
        "keep_selection": "non-Water",
        "centering_selection": "Protein",
        "num_of_step": 2,
        "index": "index with spaces.ndx",
        "pbc": "mol",
        "prefix": "run with spaces/prd",
        "skip": 1,
        "no_resnr": False,
        "output": None,
        "keep_concatenated": False,
        "preserve_time": False,
    }
    defaults.update(overrides)
    return types.SimpleNamespace(**defaults)


def test_run_uses_argv_and_stdin(monkeypatch):
    from mdtbx.trajectory import trjcat

    calls = []
    monkeypatch.setattr(
        trjcat,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs.get("input"))),
    )

    trjcat.run(_args())

    assert all(isinstance(command, list) for command, _input in calls)
    assert calls[0][1] == "non-Water\n"
    trjcat_call = next(call for call in calls if call[0][1] == "trjcat")
    assert trjcat_call[1] == "c\nc\n"
    assert "run with spaces/prd1.xtc" in trjcat_call[0]


def test_preserve_time_uses_existing_trajectory_timestamps(monkeypatch):
    from mdtbx.trajectory import trjcat

    calls = []
    monkeypatch.setattr(
        trjcat,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs.get("input"))),
    )

    trjcat.run(_args(preserve_time=True))

    trjcat_call = next(call for call in calls if call[0][1] == "trjcat")
    assert "-settime" not in trjcat_call[0]
    assert trjcat_call[1] is None


def test_cluster_supplies_three_selection_prompts(monkeypatch):
    from mdtbx.trajectory import trjcat

    calls = []
    monkeypatch.setattr(
        trjcat,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs.get("input"))),
    )

    trjcat.run(_args(pbc="cluster"))

    trjconv_call = next(call for call in calls if call[0][1] == "trjconv")
    assert trjconv_call[1] == "Protein\nProtein\nnon-Water\n"


def test_custom_output_and_intermediate_retention(tmp_path, monkeypatch):
    from mdtbx.trajectory import trjcat

    calls = []
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        trjcat,
        "run_cmd",
        lambda command, **kwargs: calls.append((command, kwargs.get("input"))),
    )
    intermediate = tmp_path / "prd_all.xtc"
    intermediate.touch()
    output = tmp_path / "custom output.xtc"

    trjcat.run(
        _args(
            prefix="prd",
            output=str(output),
            keep_concatenated=True,
        )
    )

    trjconv_call = next(call for call in calls if call[0][1] == "trjconv")
    assert str(output) in trjconv_call[0]
    assert intermediate.exists()


@pytest.mark.parametrize(
    "overrides",
    [{"num_of_step": 0}, {"skip": 0}],
)
def test_rejects_non_positive_counts(overrides):
    from mdtbx.trajectory.trjcat import run

    with pytest.raises(ValueError):
        run(_args(**overrides))
