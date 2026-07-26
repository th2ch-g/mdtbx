import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from ..logger import generate_logger
from ..utils.fep import (
    FEP_MANIFEST,
    FEP_SCHEMA_VERSION,
    format_lambdas,
    render_fep_mdp,
    temperature_mdp_overrides,
)
from ..utils.fep_rest import (
    automatic_hot_region,
    build_fep_rest_schedule,
    hot_global_indices,
    prepare_plumed_dual_state,
    restore_plumed_dual_state,
    underline_global_atoms,
    unify_fep_rest_charges,
)
from ..utils.partial_tempering import run as mark_partial_tempering
from ..utils.proc import run_cmd
from .setup_fep import _existing_file

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "setup_fep_rest",
        help="Set up a PLUMED FEP/REST Hamiltonian replica-exchange calculation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--mdp", required=True, help="Base production MDP")
    parser.add_argument("-p", "--topology", required=True, help="Dual-state topology")
    parser.add_argument("-c", "--structure", required=True, help="Starting structure")
    parser.add_argument("-o", "--outdir", default="fep_rest", help="Output directory")
    parser.add_argument("--replicas", type=int, default=32, help="Replica count")
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Physical simulation temperature [K]",
    )
    parser.add_argument(
        "--max-temperature",
        type=float,
        default=1200.0,
        help="Maximum effective REST temperature [K]",
    )
    parser.add_argument(
        "--hot-selection",
        help="Explicit mdtbx atom selection for the REST hot region",
    )
    parser.add_argument(
        "--hot-distance",
        type=float,
        default=0.4,
        help="Automatic hot-region distance from perturbed atoms [nm]",
    )
    parser.add_argument(
        "--hot-region",
        default="not water",
        help="MDTraj selection searched by automatic hot-region detection",
    )
    parser.add_argument(
        "--nstdhdl",
        type=int,
        default=100,
        help="Delta H and full-precision trajectory interval",
    )
    parser.add_argument("--checkpoint", help="Checkpoint or trajectory for grompp")
    parser.add_argument("--reference", help="Position-restraint reference structure")
    parser.add_argument("--b-reference", help="State-B restraint reference structure")
    parser.add_argument("--index", help="GROMACS index file")
    parser.add_argument("--deffnm", default="fep_rest", help="Run file basename")
    parser.add_argument(
        "--gmx",
        default="gmx_mpi",
        help="PLUMED-patched GROMACS executable",
    )
    parser.add_argument("--plumed", default="plumed", help="PLUMED executable")
    parser.add_argument(
        "--maxwarn", type=int, default=0, help="Maximum grompp warnings"
    )
    parser.add_argument(
        "--no-tpr",
        action="store_true",
        help="Generate topologies and MDP files without final TPR files",
    )
    parser.set_defaults(func=run)


def _validate_args(args):
    if args.hot_distance <= 0:
        raise ValueError("--hot-distance must be positive")
    if args.nstdhdl <= 0:
        raise ValueError("--nstdhdl must be positive")
    if args.maxwarn < 0:
        raise ValueError("--maxwarn must be non-negative")
    if not args.deffnm or Path(args.deffnm).name != args.deffnm:
        raise ValueError("--deffnm must be a filename without directories")


def _grompp(
    args,
    *,
    mdp,
    topology,
    structure,
    reference,
    checkpoint,
    b_reference,
    index,
    output,
    mdout,
    cwd,
    processed=None,
):
    command = [
        args.gmx,
        "grompp",
        "-f",
        str(mdp),
        "-p",
        str(topology),
        "-c",
        str(structure),
        "-r",
        str(reference),
        "-o",
        str(output),
        "-po",
        str(mdout),
        "-maxwarn",
        str(args.maxwarn),
    ]
    if checkpoint is not None:
        command.extend(["-t", str(checkpoint)])
    if b_reference is not None:
        command.extend(["-rb", str(b_reference)])
    if index is not None:
        command.extend(["-n", str(index)])
    if processed is not None:
        command.extend(["-pp", str(processed)])
    run_cmd(command, cwd=cwd)


def _mark_hot_region(args, processed_topology, structure, underlined_topology):
    if args.hot_selection:
        mark_partial_tempering(
            SimpleNamespace(
                topology=str(processed_topology),
                selection=args.hot_selection,
                output=str(underlined_topology),
            )
        )
    else:
        text = processed_topology.read_text()
        selected = automatic_hot_region(
            text,
            structure,
            distance=args.hot_distance,
            region_selection=args.hot_region,
        )
        underlined_topology.write_text(underline_global_atoms(text, selected))
    hot = hot_global_indices(underlined_topology.read_text())
    if not hot:
        raise ValueError("The REST hot region is empty")
    return hot


def _run_partial_tempering(args, input_text, output_path, scale):
    compatible_text, dual_records = prepare_plumed_dual_state(input_text)
    temporary_input = output_path.with_suffix(".plumed-input.top")
    temporary_output = output_path.with_suffix(".plumed-output.top")
    temporary_input.write_text(compatible_text)
    with temporary_input.open() as stdin, temporary_output.open("w") as stdout:
        run_cmd(
            [args.plumed, "partial_tempering", f"{scale:.10g}"],
            stdin=stdin,
            stdout=stdout,
        )
    scaled_text = restore_plumed_dual_state(
        temporary_output.read_text(),
        dual_records,
        scale,
    )
    output_path.write_text(scaled_text)
    temporary_input.unlink()
    temporary_output.unlink()


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

    schedule = build_fep_rest_schedule(
        args.replicas,
        args.temperature,
        args.max_temperature,
    )
    outdir = Path(args.outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    base_text = mdp_path.read_text()
    preprocess_mdp = outdir / "preprocess.mdp"
    preprocess_mdp.write_text(
        render_fep_mdp(
            base_text,
            {
                **temperature_mdp_overrides(base_text, args.temperature),
                "free-energy": "yes",
                "init-lambda-state": 0,
                "fep-lambdas": "0.000000 1.000000",
                "calc-lambda-neighbors": 1,
            },
        )
    )
    processed_topology = outdir / "processed.top"
    _grompp(
        args,
        mdp=preprocess_mdp,
        topology=topology,
        structure=structure,
        reference=reference,
        checkpoint=checkpoint,
        b_reference=b_reference,
        index=index,
        output=outdir / "preprocess.tpr",
        mdout=outdir / "preprocess-mdout.mdp",
        processed=processed_topology,
        cwd=topology.parent,
    )

    underlined_topology = outdir / "hot.top"
    hot = _mark_hot_region(
        args,
        processed_topology,
        structure,
        underlined_topology,
    )
    plumed_file = outdir / "plumed.dat"
    plumed_file.write_text("# PLUMED is enabled for Hamiltonian replica exchange.\n")

    fep_values = format_lambdas(schedule["fep"])
    vdw_values = format_lambdas(schedule["vdw"])
    coul_values = format_lambdas([0.0] * args.replicas)
    windows = []
    width = max(3, len(str(args.replicas - 1)))
    underlined_text = underlined_topology.read_text()

    for state in range(args.replicas):
        directory_name = f"lambda_{state:0{width}d}"
        window_dir = outdir / directory_name
        window_dir.mkdir()
        topology_text = unify_fep_rest_charges(
            underlined_text,
            general_lambda=schedule["fep"][state],
            charge_a_lambda=schedule["charge_a"][state],
            charge_b_lambda=schedule["charge_b"][state],
        )
        output_topology = window_dir / "fep_rest.top"
        _run_partial_tempering(
            args,
            topology_text,
            output_topology,
            schedule["scales"][state],
        )

        output_mdp = window_dir / f"{args.deffnm}.mdp"
        output_mdp.write_text(
            render_fep_mdp(
                base_text,
                {
                    **temperature_mdp_overrides(base_text, args.temperature),
                    "free-energy": "yes",
                    "init-lambda-state": state,
                    "fep-lambdas": fep_values,
                    "vdw-lambdas": vdw_values,
                    "coul-lambdas": coul_values,
                    "calc-lambda-neighbors": 1,
                    "dhdl-derivatives": "no",
                    "dhdl-print-energy": "potential",
                    "separate-dhdl-file": "yes",
                    "nstdhdl": args.nstdhdl,
                    "nstxout": args.nstdhdl,
                    "nstxout-compressed": 0,
                    "nstvout": 0,
                },
            )
        )
        if not args.no_tpr:
            _grompp(
                args,
                mdp=output_mdp,
                topology=output_topology,
                structure=structure,
                reference=reference,
                checkpoint=checkpoint,
                b_reference=b_reference,
                index=index,
                output=window_dir / f"{args.deffnm}.tpr",
                mdout=window_dir / "mdout.mdp",
                cwd=window_dir,
            )
        windows.append(
            {
                "index": state,
                "directory": directory_name,
                "fep_lambda": schedule["fep"][state],
                "vdw_lambda": schedule["vdw"][state],
                "charge_a_lambda": schedule["charge_a"][state],
                "charge_b_lambda": schedule["charge_b"][state],
                "rest_temperature": schedule["rest_temperatures"][state],
                "rest_scale": schedule["scales"][state],
            }
        )

    manifest = {
        "schema_version": FEP_SCHEMA_VERSION,
        "workflow": "fep-rest",
        "mode": "transform",
        "topology": str(topology),
        "structure": str(structure),
        "deffnm": args.deffnm,
        "prepared": not args.no_tpr,
        "physical_temperature": args.temperature,
        "maximum_rest_temperature": args.max_temperature,
        "hot_atom_indices": [index + 1 for index in hot],
        "lambda_components": schedule,
        "plumed_file": str(plumed_file),
        "plumed_executable": args.plumed,
        "gmx_executable": args.gmx,
        "hrex": True,
        "windows": windows,
    }
    (outdir / FEP_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Generated %d PLUMED FEP/REST replicas in %s", args.replicas, outdir)
