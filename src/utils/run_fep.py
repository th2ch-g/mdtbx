import argparse
import shlex
from pathlib import Path

from ..logger import generate_logger
from .fep import load_fep_manifest
from .proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "run_fep",
        help="Run prepared GROMACS FEP windows",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--path",
        default="fep",
        type=str,
        help="FEP setup directory or manifest",
    )
    parser.add_argument(
        "--start-window",
        default=0,
        type=int,
        help="First window index",
    )
    parser.add_argument(
        "--end-window",
        type=int,
        help="Last window index, inclusive",
    )
    parser.add_argument(
        "--multidir",
        action="store_true",
        help="Run selected windows in one GROMACS multidir invocation",
    )
    parser.add_argument(
        "--replex",
        default=0,
        type=int,
        help="Replica-exchange interval for multidir runs",
    )
    parser.add_argument(
        "--ntmpi",
        type=int,
        help="MPI rank count for the launcher or thread-MPI mdrun",
    )
    parser.add_argument(
        "--ntomp",
        type=int,
        help="OpenMP thread count per rank",
    )
    parser.add_argument(
        "--maxh",
        type=float,
        help="Maximum run time [hours]",
    )
    parser.add_argument(
        "--nsteps",
        type=int,
        help="Override the number of integration steps",
    )
    parser.add_argument(
        "--no-continue",
        action="store_true",
        help="Do not continue from existing checkpoints",
    )
    parser.add_argument(
        "--gmx",
        type=str,
        help="GROMACS executable; FEP/REST defaults to the setup executable",
    )
    parser.add_argument(
        "--mpi-launcher",
        help=(
            "External MPI launcher template, for example "
            "'mpirun -np {ntmpi}' or 'srun --ntasks {ntmpi}'"
        ),
    )
    parser.set_defaults(func=run)


def _selected_windows(base, manifest, start, end):
    windows = manifest["windows"]
    last = len(windows) - 1 if end is None else end
    if start < 0 or last < start or last >= len(windows):
        raise ValueError(
            f"Window range must satisfy 0 <= start <= end < {len(windows)}"
        )
    return [
        (window["index"], (base / window["directory"]).resolve())
        for window in windows[start : last + 1]
    ]


def _runtime_options(args, *, external_mpi=False):
    options = []
    if args.ntmpi is not None:
        if args.ntmpi <= 0:
            raise ValueError("--ntmpi must be positive")
        if not external_mpi:
            options.extend(["-ntmpi", str(args.ntmpi)])
    if args.ntomp is not None:
        if args.ntomp <= 0:
            raise ValueError("--ntomp must be positive")
        options.extend(["-ntomp", str(args.ntomp)])
    if args.maxh is not None:
        if args.maxh <= 0:
            raise ValueError("--maxh must be positive")
        options.extend(["-maxh", str(args.maxh)])
    if args.nsteps is not None:
        if args.nsteps < -2:
            raise ValueError("--nsteps must be -2 or greater")
        options.extend(["-nsteps", str(args.nsteps)])
    return options


def _mpi_prefix(args, *, is_fep_rest, window_count, multidir):
    template = getattr(args, "mpi_launcher", None)
    if template is None and not is_fep_rest:
        return []
    if template is None:
        template = "mpirun -np {ntmpi}"
    ranks = args.ntmpi
    if ranks is None:
        ranks = window_count if multidir else 1
    if ranks <= 0:
        raise ValueError("--ntmpi must be positive")
    if multidir and (ranks < window_count or ranks % window_count):
        raise ValueError(
            "External MPI rank count must be a positive multiple of the "
            "selected window count"
        )
    try:
        command = template.format(
            ntmpi=ranks,
            ntomp=args.ntomp if args.ntomp is not None else 1,
        )
    except (KeyError, ValueError) as error:
        raise ValueError(f"Invalid --mpi-launcher template: {error}") from error
    prefix = shlex.split(command)
    if not prefix:
        raise ValueError("--mpi-launcher must not be empty")
    return prefix


def run(args):
    base, manifest = load_fep_manifest(args.path)
    selected = _selected_windows(
        base,
        manifest,
        args.start_window,
        args.end_window,
    )
    deffnm = manifest["deffnm"]
    is_fep_rest = manifest.get("workflow") == "fep-rest"
    mpi_prefix = _mpi_prefix(
        args,
        is_fep_rest=is_fep_rest,
        window_count=len(selected),
        multidir=args.multidir,
    )
    runtime_options = _runtime_options(args, external_mpi=bool(mpi_prefix))
    gmx = args.gmx or manifest.get("gmx_executable", "gmx")

    for index, directory in selected:
        tpr = directory / f"{deffnm}.tpr"
        if not tpr.is_file():
            raise FileNotFoundError(f"TPR for window {index} not found: {tpr}")

    if args.replex < 0:
        raise ValueError("--replex must be non-negative")
    if args.replex and not args.multidir:
        raise ValueError("--replex requires --multidir")
    if is_fep_rest and not args.multidir:
        raise ValueError("FEP/REST must be run with --multidir")
    if is_fep_rest and not args.replex:
        raise ValueError("FEP/REST requires a positive --replex interval")

    if args.multidir:
        checkpoints = [
            (directory / f"{deffnm}.cpt").is_file() for _index, directory in selected
        ]
        if not args.no_continue and any(checkpoints) and not all(checkpoints):
            raise ValueError(
                "All selected windows must have checkpoints for multidir continuation"
            )
        command = [
            *mpi_prefix,
            gmx,
            "mdrun",
            "-multidir",
            *[str(directory) for _index, directory in selected],
            "-s",
            f"{deffnm}.tpr",
            "-deffnm",
            deffnm,
            *runtime_options,
        ]
        if not args.no_continue and all(checkpoints):
            command.extend(["-cpi", f"{deffnm}.cpt"])
        if args.replex:
            command.extend(["-replex", str(args.replex)])
        if is_fep_rest:
            plumed_file = Path(manifest["plumed_file"]).expanduser().resolve()
            if not plumed_file.is_file():
                raise FileNotFoundError(f"PLUMED input not found: {plumed_file}")
            command.extend(["-plumed", str(plumed_file), "-hrex", "-dlb", "no"])
        run_cmd(command, log=f"Finished {len(selected)} FEP windows")
        return

    for index, directory in selected:
        command = [
            *mpi_prefix,
            gmx,
            "mdrun",
            "-s",
            f"{deffnm}.tpr",
            "-deffnm",
            deffnm,
            *runtime_options,
        ]
        checkpoint = directory / f"{deffnm}.cpt"
        if not args.no_continue and checkpoint.is_file():
            command.extend(["-cpi", checkpoint.name])
        run_cmd(command, cwd=directory, log=f"Finished FEP window {index}")
