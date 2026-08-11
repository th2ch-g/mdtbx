import json
from types import SimpleNamespace

from mdtbx.analysis import optimize_fep_schedule


def test_run_writes_reusable_schedule(tmp_path):
    base = tmp_path / "fep"
    first = base / "lambda_000"
    first.mkdir(parents=True)
    manifest = {
        "schema_version": 1,
        "mode": "transform",
        "deffnm": "fep",
        "lambda_components": {"fep": [0.0, 0.5, 1.0]},
        "windows": [
            {"index": 0, "directory": "lambda_000"},
            {"index": 1, "directory": "lambda_001"},
            {"index": 2, "directory": "lambda_002"},
        ],
    }
    (base / "fep_manifest.json").write_text(json.dumps(manifest))
    (first / "fep.log").write_text(
        "Repl  average probabilities:\nRepl  0 1 2\nRepl  .10 .80\n"
    )

    optimize_fep_schedule.run(
        SimpleNamespace(
            path=str(base),
            log=None,
            iteration=1,
            min_probability=0.01,
            output=None,
        )
    )

    result = json.loads((base / "optimized_schedule.json").read_text())
    assert result["workflow"] == "optimized-fep-schedule"
    assert result["exchange_probabilities"] == [0.1, 0.8]
    assert result["coordinates"][1] < 0.5
