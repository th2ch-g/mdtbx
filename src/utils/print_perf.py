import argparse
import re
from pathlib import Path
import polars as pl

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx print_perf
    """
    parser = subparsers.add_parser(
        "print_perf",
        help="Print performance",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p", "--path", default=".", type=str, help="Path to the directory"
    )

    parser.add_argument(
        "--prefix", help="prefix of gromacs log file", type=str, default="prd"
    )

    parser.add_argument(
        "-o", "--output", help="Output file name", type=str,
    )


def parse_log_file(log_path):
    data = {
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
        with open(log_path, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if "Command line" in line:
                if i + 1 < len(lines):
                    cmd_full = lines[i + 1].strip()
                    # Remove -deffnm option and its argument, which often vary between runs
                    cmd_normalized = re.sub(r"-deffnm\s+\S+", "", cmd_full).strip()
                    # Replace multiple spaces with a single space
                    cmd_normalized = re.sub(r"\s\s+", " ", cmd_normalized)
                    data["cmd"] = cmd_normalized
            elif "GROMACS version:" in line:
                data["version"] = line.split()[2]
            elif "Executable:" in line:
                data["executable"] = line.split()[1]
            elif "Hardware detected on host" in line:
                data["hostname"] = line.split()[4].split(":")[0]
            elif "GPU info:" in line:
                if i + 2 < len(lines):
                    n_gpu_line = lines[i + 1].strip()
                    gpu_info_line = lines[i + 2].strip()
                    data["n_GPU"] = n_gpu_line.split()[4]
                    data["GPU_info"] = gpu_info_line
            elif "CPU info:" in line:
                if i + 2 < len(lines):
                    cpu_info_line = lines[i + 2].strip()
                    data["CPU_info"] = " ".join(cpu_info_line.split()[1:])
            elif "Performance:" in line:
                data["performance"] = float(line.split()[1])

    except (IOError, IndexError) as e:
        LOGGER.error(f"Error parsing {log_path}: {e}")
        return None

    return data


def run(args):
    log_files = list(Path(args.path).glob(f"{args.prefix}*.log"))
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
            pl.count("performance").alias("count"),
            pl.first("version"),
            pl.first("hostname"),
            pl.first("n_GPU"),
            pl.first("GPU_info"),
            pl.first("CPU_info"),
        )
        .sort("cmd")
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
            agg_df.write_csv(args.output)
            LOGGER.info(f"Performance data successfully saved to {args.output}")
        except IOError as e:
            LOGGER.error(f"Error writing to file {args.output}: {e}")
