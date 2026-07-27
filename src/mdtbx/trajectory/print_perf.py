import argparse
import math
import re
from pathlib import Path

import polars as pl

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx print_perf
    """
    parser = subparsers.add_parser(
        "print_perf",
        help="Print performance (ns/day)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p", "--path", default=".", type=str, help="Path to the directory"
    )

    parser.add_argument(
        "--prefix", help="prefix of gromacs log file", type=str, default="prd"
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output file name",
        type=str,
    )

    parser.set_defaults(func=run)


def _value_after_colon(line: str) -> str | None:
    _, separator, value = line.partition(":")
    value = value.strip()
    return value if separator and value else None


def parse_log_file(log_path: str | Path) -> dict[str, str | float | None] | None:
    data: dict[str, str | float | None] = {
        "cmd": "N/A",
        "version": "N/A",
        "executable": "N/A",
        "hostname": "N/A",
        "n_GPU": "N/A",
        "GPU_info": "N/A",
        "CPU_info": "N/A",
        "performance": None,
    }

    try:
        with Path(log_path).open(errors="replace") as f:
            lines = f.readlines()
    except OSError as exc:
        LOGGER.error("Error reading %s: %s", log_path, exc)
        return None

    for i, line in enumerate(lines):
        if "Command line" in line and i + 1 < len(lines):
            cmd_full = lines[i + 1].strip()
            cmd_normalized = re.sub(r"(?:^|\s)-deffnm(?:\s+\S+|=\S+)", " ", cmd_full)
            data["cmd"] = " ".join(cmd_normalized.split()) or "N/A"
        elif "GROMACS version:" in line:
            data["version"] = _value_after_colon(line) or "N/A"
        elif "Executable:" in line:
            data["executable"] = _value_after_colon(line) or "N/A"
        elif "Hardware detected on host" in line:
            match = re.search(r"Hardware detected on host\s+([^:\s]+)", line)
            if match:
                data["hostname"] = match.group(1)
        elif "GPU info:" in line:
            gpu_block = lines[i + 1 : i + 12]
            gpu_descriptions = []
            for gpu_line in gpu_block:
                if "CPU info:" in gpu_line:
                    break
                count_match = re.search(r"Number of GPUs detected:\s*(\d+)", gpu_line)
                if count_match:
                    data["n_GPU"] = count_match.group(1)
                if re.match(r"\s*GPU\s+\d+\s*:", gpu_line):
                    gpu_descriptions.append(gpu_line.strip())
            if gpu_descriptions:
                data["GPU_info"] = "; ".join(gpu_descriptions)
        elif "CPU info:" in line:
            cpu_block = lines[i + 1 : i + 20]
            for cpu_line in cpu_block:
                if "Command line" in cpu_line:
                    break
                if "Model name:" in cpu_line:
                    data["CPU_info"] = _value_after_colon(cpu_line) or "N/A"
                    break
        elif "Performance:" in line:
            match = re.search(r"Performance:\s+(\S+)", line)
            if not match:
                continue
            try:
                performance = float(match.group(1))
            except ValueError:
                LOGGER.warning(
                    "Ignoring malformed performance value in %s: %s",
                    log_path,
                    match.group(1),
                )
                continue
            if math.isfinite(performance) and performance >= 0:
                data["performance"] = performance
            else:
                LOGGER.warning(
                    "Ignoring non-finite or negative performance in %s: %s",
                    log_path,
                    match.group(1),
                )

    return data


def run(args):
    log_dir = Path(args.path)
    log_files = sorted(
        path for path in log_dir.glob("*.log") if path.name.startswith(args.prefix)
    )
    if not log_files:
        LOGGER.warning(f"No log files with prefix '{args.prefix}' found in {args.path}")
        return

    all_data = [
        d
        for d in (parse_log_file(p) for p in log_files)
        if d and d["performance"] is not None
    ]

    if not all_data:
        LOGGER.warning("No performance data could be extracted from the log files.")
        return

    df = pl.DataFrame(all_data)

    agg_df = (
        df.group_by("cmd")
        .agg(
            pl.mean("performance").alias("mean_perf"),
            pl.std("performance").alias("std_perf"),
            pl.len().alias("count"),
            pl.first("version"),
            pl.first("hostname"),
            pl.first("n_GPU"),
            pl.first("GPU_info"),
            pl.first("CPU_info"),
        )
        .sort("cmd")
        .with_columns(pl.col("std_perf").fill_null(0.0))
    )

    print(f"{'Command':<80} {'Mean (ns/day)':>15} {'Std (ns/day)':>15} {'Count':>7}")
    print("=" * 120)

    for row in agg_df.iter_rows(named=True):
        print(
            f"{row['cmd']:<80} {row['mean_perf']:>15.2f} {row['std_perf']:>15.2f} {row['count']:>7}"
        )
        print(
            f"  - Version: {row['version']}, Host: {row['hostname']}, GPUs: {row['n_GPU']}"
        )
        print(f"  - GPU Info: {row['GPU_info']}")
        print(f"  - CPU Info: {row['CPU_info']}")
        print("-" * 120)

    if args.output:
        try:
            output_path = Path(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            agg_df.write_csv(output_path)
            LOGGER.info("Performance data successfully saved to %s", output_path)
        except OSError as exc:
            raise RuntimeError(
                f"Failed to write performance data to {args.output}: {exc}"
            ) from exc
