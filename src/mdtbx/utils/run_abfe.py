import argparse
from types import SimpleNamespace

from ..logger import generate_logger
from .abfe import load_abfe_manifest
from . import run_fep

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "run_abfe",
        help="Run all prepared ABFE thermodynamic-cycle legs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--path", default="abfe", help="ABFE setup directory")
    parser.add_argument(
        "--legs",
        nargs="+",
        choices=[
            "solvent_charge",
            "solvent_vdw",
            "complex_charge",
            "complex_vdw",
            "complex_restraint",
        ],
        help="Selected legs; all legs are run by default",
    )
    parser.add_argument("--start-window", type=int, default=0)
    parser.add_argument("--end-window", type=int)
    parser.add_argument("--multidir", action="store_true")
    parser.add_argument("--replex", type=int, default=0)
    parser.add_argument(
        "--ntmpi",
        type=int,
        help="MPI rank count for the launcher or thread-MPI mdrun",
    )
    parser.add_argument("--ntomp", type=int)
    parser.add_argument("--maxh", type=float)
    parser.add_argument("--nsteps", type=int)
    parser.add_argument("--no-continue", action="store_true")
    parser.add_argument(
        "--gmx",
        help="GROMACS executable; defaults to the setup executable",
    )
    parser.add_argument(
        "--mpi-launcher",
        help=(
            "External MPI launcher template, for example "
            "'mpirun -np {ntmpi}' or 'srun --ntasks {ntmpi}'"
        ),
    )
    parser.add_argument(
        "--gpu-offload",
        action="store_true",
        help="Enable the validated single-GPU GROMACS offload preset",
    )
    parser.set_defaults(func=run)


def run(args):
    base, manifest = load_abfe_manifest(args.path)
    legs = args.legs if args.legs else list(manifest["legs"])
    for leg in legs:
        run_fep.run(
            SimpleNamespace(
                path=str(base / manifest["legs"][leg]),
                start_window=args.start_window,
                end_window=args.end_window,
                multidir=args.multidir,
                replex=args.replex,
                ntmpi=args.ntmpi,
                ntomp=args.ntomp,
                maxh=args.maxh,
                nsteps=args.nsteps,
                no_continue=args.no_continue,
                gmx=args.gmx,
                mpi_launcher=args.mpi_launcher,
                gpu_offload=getattr(args, "gpu_offload", False),
            )
        )
        LOGGER.info("Finished ABFE leg %s", leg)
