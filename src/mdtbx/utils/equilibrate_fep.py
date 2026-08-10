"""Run staged state-A equilibration for a dual-state topology."""

import argparse
import json
from pathlib import Path

from ..logger import generate_logger
from .fep import render_fep_mdp, temperature_mdp_overrides
from .proc import run_cmd

LOGGER = generate_logger(__name__)
_GPU_DYNAMICAL_INTEGRATORS = {"bd", "md", "md-vv", "md-vv-avek", "sd"}


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "equilibrate_fep",
        help="Run staged equilibration at one endpoint of a dual-state topology",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-p", "--topology", required=True)
    parser.add_argument("-c", "--structure", required=True)
    parser.add_argument("-f", "--mdps", required=True, nargs="+")
    parser.add_argument("-o", "--outdir", default="equilibration")
    parser.add_argument("--reference")
    parser.add_argument("--index")
    parser.add_argument("--lambda-state", type=int, choices=[0, 1], default=0)
    parser.add_argument(
        "--temperature",
        type=float,
        help="Override ref_t and gen_temp in each stage [K]",
    )
    parser.add_argument("--gmx", default="gmx_mpi")
    parser.add_argument("--ntomp", type=int, default=1)
    parser.add_argument("--maxwarn", type=int, default=0)
    parser.add_argument("--gpu-offload", action="store_true")
    parser.set_defaults(func=run)


def _existing(path, label):
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} not found: {resolved}")
    return resolved


def _mdp_integrator(text):
    for raw in text.splitlines():
        code = raw.partition(";")[0]
        if "=" not in code:
            continue
        key, value = code.split("=", 1)
        if key.strip().lower().replace("_", "-") == "integrator":
            return value.strip().lower()
    return None


def _gpu_options(enabled, integrator):
    if not enabled or integrator not in _GPU_DYNAMICAL_INTEGRATORS:
        return []
    return [
        "-pin",
        "on",
        "-nb",
        "gpu",
        "-pme",
        "gpu",
        "-pmefft",
        "gpu",
        "-bonded",
        "gpu",
        "-update",
        "cpu",
    ]


def run(args):
    if args.ntomp <= 0:
        raise ValueError("--ntomp must be positive")
    if args.maxwarn < 0:
        raise ValueError("--maxwarn must be non-negative")
    topology = _existing(args.topology, "Topology")
    structure = _existing(args.structure, "Structure")
    reference = _existing(args.reference, "Reference") if args.reference else structure
    index = _existing(args.index, "Index") if args.index else None
    mdps = [_existing(path, "MDP") for path in args.mdps]
    outdir = Path(args.outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    current_structure = structure
    current_checkpoint = None
    stages = []
    lambda_values = "0.000000 1.000000"
    for number, source_mdp in enumerate(mdps):
        stage_name = f"{number:02d}_{source_mdp.stem}"
        stage_dir = outdir / stage_name
        stage_dir.mkdir()
        mdp = stage_dir / "stage.mdp"
        base_text = source_mdp.read_text()
        temperature_settings = (
            temperature_mdp_overrides(base_text, args.temperature)
            if args.temperature is not None
            else {}
        )
        mdp.write_text(
            render_fep_mdp(
                base_text,
                {
                    **temperature_settings,
                    "free-energy": "yes",
                    "init-lambda-state": args.lambda_state,
                    "fep-lambdas": lambda_values,
                    "calc-lambda-neighbors": 1,
                    "dhdl-derivatives": "no",
                    "dhdl-print-energy": "potential",
                    "separate-dhdl-file": "yes",
                    "nstdhdl": 1000,
                },
            )
        )
        integrator = _mdp_integrator(mdp.read_text())
        tpr = stage_dir / "stage.tpr"
        grompp = [
            args.gmx,
            "grompp",
            "-f",
            str(mdp),
            "-p",
            str(topology),
            "-c",
            str(current_structure),
            "-r",
            str(reference),
            "-o",
            str(tpr),
            "-po",
            str(stage_dir / "mdout.mdp"),
            "-maxwarn",
            str(args.maxwarn),
        ]
        if current_checkpoint is not None and current_checkpoint.is_file():
            grompp.extend(["-t", str(current_checkpoint)])
        if index is not None:
            grompp.extend(["-n", str(index)])
        run_cmd(grompp, cwd=topology.parent)

        mdrun = [
            args.gmx,
            "mdrun",
            "-s",
            "stage.tpr",
            "-deffnm",
            "stage",
            "-ntomp",
            str(args.ntomp),
            *_gpu_options(args.gpu_offload, integrator),
        ]
        run_cmd(mdrun, cwd=stage_dir)
        current_structure = stage_dir / "stage.gro"
        if not current_structure.is_file():
            raise FileNotFoundError(
                f"Equilibration output not found: {current_structure}"
            )
        current_checkpoint = stage_dir / "stage.cpt"
        stages.append(
            {
                "index": number,
                "name": stage_name,
                "source_mdp": str(source_mdp),
                "structure": str(current_structure),
                "checkpoint": str(current_checkpoint)
                if current_checkpoint.is_file()
                else None,
            }
        )

    manifest = {
        "schema_version": 1,
        "workflow": "equilibrate-fep",
        "lambda_state": args.lambda_state,
        "temperature_k": args.temperature,
        "topology": str(topology),
        "initial_structure": str(structure),
        "final_structure": str(current_structure),
        "final_checkpoint": str(current_checkpoint)
        if current_checkpoint and current_checkpoint.is_file()
        else None,
        "stages": stages,
    }
    (outdir / "equilibration_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    LOGGER.info("Completed %d equilibration stages in %s", len(stages), outdir)
