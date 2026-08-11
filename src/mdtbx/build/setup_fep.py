import argparse
import json
from pathlib import Path

from ..logger import generate_logger
from ..utils.fep import (
    FEP_MANIFEST,
    FEP_SCHEMA_VERSION,
    build_lambda_schedule,
    format_lambdas,
    render_fep_mdp,
    soft_core_mdp_settings,
)
from ..utils.fep_schedule import load_optimized_schedule
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "setup_fep",
        help="Set up GROMACS alchemical FEP windows",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-f",
        "--mdp",
        required=True,
        type=str,
        help="Base production MDP file",
    )
    parser.add_argument(
        "-p",
        "--topology",
        required=True,
        type=str,
        help="GROMACS topology file",
    )
    parser.add_argument(
        "-c",
        "--structure",
        required=True,
        type=str,
        help="Starting structure file",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="fep",
        type=str,
        help="Output directory",
    )
    parser.add_argument(
        "--mode",
        choices=["decouple", "charge", "vdw", "transform"],
        default="decouple",
        help="Alchemical transformation mode",
    )
    parser.add_argument(
        "--moltype",
        type=str,
        help="GROMACS molecule type to decouple",
    )
    parser.add_argument(
        "--windows",
        default=12,
        type=int,
        help="Window count for charge, vdw, or transform mode",
    )
    parser.add_argument(
        "--coul-windows",
        default=12,
        type=int,
        help="Electrostatic window count for decouple mode",
    )
    parser.add_argument(
        "--vdw-windows",
        default=12,
        type=int,
        help="Van der Waals window count for decouple mode",
    )
    parser.add_argument(
        "--fep-lambdas",
        nargs="+",
        type=float,
        help="Custom dual-state lambda schedule",
    )
    parser.add_argument(
        "--coul-lambdas",
        nargs="+",
        type=float,
        help="Custom electrostatic lambda schedule",
    )
    parser.add_argument(
        "--vdw-lambdas",
        nargs="+",
        type=float,
        help="Custom Van der Waals lambda schedule",
    )
    parser.add_argument(
        "--schedule",
        help="Optimized FEP schedule JSON",
    )
    parser.add_argument(
        "--calc-lambda-neighbors",
        default=1,
        type=int,
        help="Number of neighboring states evaluated for Delta H",
    )
    parser.add_argument(
        "--nstdhdl",
        default=100,
        type=int,
        help="Delta H output interval",
    )
    parser.add_argument(
        "--checkpoint",
        type=str,
        help="Checkpoint or trajectory supplied to grompp",
    )
    parser.add_argument(
        "--reference",
        type=str,
        help="Position-restraint reference structure",
    )
    parser.add_argument(
        "--b-reference",
        type=str,
        help="State-B restraint reference structure",
    )
    parser.add_argument(
        "--index",
        type=str,
        help="GROMACS index file",
    )
    parser.add_argument(
        "--deffnm",
        default="fep",
        type=str,
        help="Run file basename",
    )
    parser.add_argument(
        "--gmx",
        default="gmx",
        type=str,
        help="GROMACS executable",
    )
    parser.add_argument(
        "--maxwarn",
        default=0,
        type=int,
        help="Maximum grompp warnings",
    )
    parser.add_argument(
        "--no-grompp",
        action="store_true",
        help="Generate MDP files without creating TPR files",
    )
    parser.set_defaults(func=run)


def _existing_file(path, label):
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} not found: {resolved}")
    return resolved


def _validate_args(args):
    if args.mode != "transform" and not args.moltype:
        raise ValueError("--moltype is required for decouple, charge, and vdw modes")
    if args.mode == "transform" and args.moltype:
        raise ValueError("--moltype is not used in transform mode")
    if args.calc_lambda_neighbors == 0 or args.calc_lambda_neighbors < -1:
        raise ValueError("--calc-lambda-neighbors must be -1 or positive")
    if args.nstdhdl <= 0:
        raise ValueError("--nstdhdl must be positive")
    if args.maxwarn < 0:
        raise ValueError("--maxwarn must be non-negative")
    if not args.deffnm or Path(args.deffnm).name != args.deffnm:
        raise ValueError("--deffnm must be a filename without directories")
    if getattr(args, "schedule", None) and any(
        getattr(args, name, None) is not None
        for name in ("fep_lambdas", "coul_lambdas", "vdw_lambdas")
    ):
        raise ValueError("--schedule cannot be combined with explicit lambda lists")


def _fep_settings(args, schedule, state):
    settings = {
        "free-energy": "yes",
        "init-lambda-state": state,
        "calc-lambda-neighbors": args.calc_lambda_neighbors,
        "dhdl-derivatives": "no",
        "dhdl-print-energy": "potential",
        "separate-dhdl-file": "yes",
        "nstdhdl": args.nstdhdl,
    }
    if "fep" in schedule:
        settings["fep-lambdas"] = format_lambdas(schedule["fep"])
        settings.update(soft_core_mdp_settings(coulomb=True))
        return settings

    settings["coul-lambdas"] = format_lambdas(schedule["coul"])
    settings["vdw-lambdas"] = format_lambdas(schedule["vdw"])
    settings["couple-moltype"] = args.moltype
    settings["couple-lambda0"] = "vdw-q" if args.mode != "vdw" else "vdw"
    settings["couple-lambda1"] = "vdw" if args.mode == "charge" else "none"
    settings["couple-intramol"] = "no"
    settings.update(soft_core_mdp_settings(coulomb=False))
    return settings


def _grompp_command(
    args,
    *,
    mdp_path,
    topology,
    structure,
    reference,
    checkpoint,
    b_reference,
    index,
    window_dir,
):
    command = [
        args.gmx,
        "grompp",
        "-f",
        str(mdp_path),
        "-p",
        str(topology),
        "-c",
        str(structure),
        "-r",
        str(reference),
        "-o",
        str(window_dir / f"{args.deffnm}.tpr"),
        "-po",
        str(window_dir / "mdout.mdp"),
        "-maxwarn",
        str(args.maxwarn),
    ]
    if checkpoint is not None:
        command.extend(["-t", str(checkpoint)])
    if b_reference is not None:
        command.extend(["-rb", str(b_reference)])
    if index is not None:
        command.extend(["-n", str(index)])
    return command


def run(args):
    _validate_args(args)
    mdp_path = _existing_file(args.mdp, "MDP file")
    topology = _existing_file(args.topology, "Topology")
    structure = _existing_file(args.structure, "Structure")
    reference = (
        _existing_file(args.reference, "Reference structure")
        if args.reference
        else structure
    )
    checkpoint = (
        _existing_file(args.checkpoint, "Checkpoint") if args.checkpoint else None
    )
    b_reference = (
        _existing_file(args.b_reference, "State-B reference")
        if args.b_reference
        else None
    )
    index = _existing_file(args.index, "Index") if args.index else None

    schedule_source = None
    if getattr(args, "schedule", None):
        schedule_source, schedule_data = load_optimized_schedule(
            args.schedule,
            expected_mode=args.mode,
            expected_workflow="fep",
        )
        schedule = schedule_data["lambda_components"]
    else:
        schedule = build_lambda_schedule(
            args.mode,
            windows=args.windows,
            coul_windows=args.coul_windows,
            vdw_windows=args.vdw_windows,
            fep_lambdas=args.fep_lambdas,
            coul_lambdas=args.coul_lambdas,
            vdw_lambdas=args.vdw_lambdas,
        )
    window_count = len(next(iter(schedule.values())))
    outdir = Path(args.outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    base_mdp = mdp_path.read_text()
    windows = []
    width = max(3, len(str(window_count - 1)))
    for state in range(window_count):
        directory_name = f"lambda_{state:0{width}d}"
        window_dir = outdir / directory_name
        window_dir.mkdir()
        output_mdp = window_dir / f"{args.deffnm}.mdp"
        settings = _fep_settings(args, schedule, state)
        output_mdp.write_text(render_fep_mdp(base_mdp, settings))

        if not args.no_grompp:
            command = _grompp_command(
                args,
                mdp_path=output_mdp,
                topology=topology,
                structure=structure,
                reference=reference,
                checkpoint=checkpoint,
                b_reference=b_reference,
                index=index,
                window_dir=window_dir,
            )
            run_cmd(
                command,
                cwd=topology.parent,
                log=f"Prepared FEP window {state}",
            )

        window = {"index": state, "directory": directory_name}
        for component, values in schedule.items():
            window[f"{component}_lambda"] = values[state]
        windows.append(window)

    manifest = {
        "schema_version": FEP_SCHEMA_VERSION,
        "mode": args.mode,
        "molecule_type": args.moltype,
        "base_mdp": str(mdp_path),
        "topology": str(topology),
        "structure": str(structure),
        "checkpoint": str(checkpoint) if checkpoint else None,
        "reference": str(reference),
        "b_reference": str(b_reference) if b_reference else None,
        "index": str(index) if index else None,
        "deffnm": args.deffnm,
        "gmx_executable": args.gmx,
        "prepared": not args.no_grompp,
        "lambda_components": schedule,
        "schedule_source": str(schedule_source) if schedule_source else None,
        "windows": windows,
    }
    (outdir / FEP_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info(f"Generated {window_count} FEP windows in {outdir}")
