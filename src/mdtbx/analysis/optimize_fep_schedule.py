import argparse
import json
from pathlib import Path

from ..logger import generate_logger
from ..utils.fep import FEP_MANIFEST, load_fep_manifest
from ..utils.fep_schedule import (
    make_optimized_schedule,
    parse_exchange_probabilities,
)


LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "optimize_fep_schedule",
        help="Redistribute an FEP lambda schedule from replica-exchange probabilities",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--path",
        default="fep",
        help="FEP setup directory or manifest",
    )
    parser.add_argument(
        "--log",
        help="GROMACS replica-exchange log; defaults to the first window log",
    )
    parser.add_argument(
        "--iteration",
        type=int,
        default=1,
        help="One-based optimization iteration used for damping",
    )
    parser.add_argument(
        "--min-probability",
        type=float,
        default=0.01,
        help="Probability floor used when redistributing states",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output schedule JSON",
    )
    parser.set_defaults(func=run)


def run(args):
    base, manifest = load_fep_manifest(args.path)
    source_path = Path(args.path).expanduser()
    manifest_path = (
        source_path.resolve() if source_path.is_file() else base / FEP_MANIFEST
    )
    log_path = (
        Path(args.log).expanduser().resolve()
        if args.log
        else (
            base / manifest["windows"][0]["directory"] / f"{manifest['deffnm']}.log"
        ).resolve()
    )
    if not log_path.is_file():
        raise FileNotFoundError(f"Replica-exchange log not found: {log_path}")
    probabilities = parse_exchange_probabilities(
        log_path.read_text(),
        len(manifest["windows"]),
    )
    schedule = make_optimized_schedule(
        manifest,
        probabilities,
        source_manifest=manifest_path,
        source_log=log_path,
        iteration=args.iteration,
        minimum_probability=args.min_probability,
    )
    output = (
        Path(args.output).expanduser().resolve()
        if args.output
        else (base / "optimized_schedule.json").resolve()
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(schedule, indent=2) + "\n")
    LOGGER.info("Wrote optimized FEP schedule to %s", output)
