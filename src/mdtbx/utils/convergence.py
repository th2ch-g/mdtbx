import math
from pathlib import Path


def xvg_times(path):
    times = []
    with Path(path).open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped[0] in {"#", "@"}:
                continue
            try:
                value = float(stripped.split()[0])
            except (IndexError, ValueError):
                continue
            if not math.isfinite(value):
                raise ValueError(f"Non-finite time value in {path}")
            times.append(value)
    if not times:
        raise ValueError(f"No time samples found in {path}")
    if any(current > following for current, following in zip(times, times[1:])):
        raise ValueError(f"Time samples are not ordered in {path}")
    return times


def convergence_ranges(paths, begin, end, block_count):
    if block_count < 2:
        raise ValueError("--convergence-blocks must be at least 2")
    samples = [xvg_times(path) for path in dict.fromkeys(Path(path) for path in paths)]
    effective_begin = max(float(begin), *(values[0] for values in samples))
    effective_end = min(
        *(values[-1] for values in samples),
        float(end) if end is not None else math.inf,
    )
    if effective_end <= effective_begin:
        raise ValueError("FEP inputs have no common convergence-analysis time range")
    edges = [
        effective_begin + (effective_end - effective_begin) * index / block_count
        for index in range(block_count + 1)
    ]
    blocks = list(zip(edges, edges[1:]))
    for block_begin, block_end in blocks:
        for values in samples:
            count = sum(block_begin <= value <= block_end for value in values)
            if count < 2:
                raise ValueError(
                    "Each convergence block requires at least 2 samples in every input"
                )
    return {
        "effective_begin_ps": effective_begin,
        "effective_end_ps": effective_end,
        "blocks": blocks,
        "cumulative": [(effective_begin, block_end) for block_end in edges[1:]],
    }
