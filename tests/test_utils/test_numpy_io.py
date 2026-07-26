from pathlib import Path

import numpy as np

from src.utils.numpy_io import save_npy, save_npz


def test_save_npy_creates_parent_and_adds_suffix(tmp_path):
    requested_path = tmp_path / "nested" / "values"

    output_path = save_npy(requested_path, np.array([1.0, 2.0]))

    assert output_path == Path(f"{requested_path}.npy")
    assert np.load(output_path).tolist() == [1.0, 2.0]


def test_save_npz_creates_parent_and_adds_suffix(tmp_path):
    requested_path = tmp_path / "nested" / "values"

    output_path = save_npz(requested_path, values=np.array([1.0, 2.0]))

    assert output_path == Path(f"{requested_path}.npz")
    with np.load(output_path) as archive:
        assert archive["values"].tolist() == [1.0, 2.0]
