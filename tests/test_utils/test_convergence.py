import pytest

from mdtbx.utils.convergence import convergence_ranges


def _xvg(path, times):
    path.write_text('@ title "test"\n' + "".join(f"{time} 0\n" for time in times))
    return path


def test_convergence_ranges_use_common_equal_duration_interval(tmp_path):
    first = _xvg(tmp_path / "first.xvg", [0, 1, 2, 3, 4])
    second = _xvg(tmp_path / "second.xvg", [1, 2, 3, 4, 5])

    ranges = convergence_ranges([first, second], 0.0, None, 3)

    assert ranges["effective_begin_ps"] == 1.0
    assert ranges["effective_end_ps"] == 4.0
    assert ranges["blocks"] == [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0)]


def test_convergence_ranges_require_two_samples_per_block(tmp_path):
    sparse = _xvg(tmp_path / "sparse.xvg", [0, 20])

    with pytest.raises(ValueError, match="at least 2 samples"):
        convergence_ranges([sparse], 0.0, None, 2)
