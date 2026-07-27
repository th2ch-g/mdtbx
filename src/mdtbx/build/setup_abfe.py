import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from ..logger import generate_logger
from ..utils.abfe import (
    ABFE_MANIFEST,
    ABFE_SCHEMA_VERSION,
    boresch_pull_settings,
    calculate_anchor_geometry,
    select_boresch_anchors,
    write_anchor_index,
)
from ..utils.fep import render_fep_mdp, temperature_mdp_overrides
from ..utils.proc import run_cmd
from . import setup_fep

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "setup_abfe",
        help="Set up a Boresch-restraint absolute binding free-energy cycle",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--mdp", required=True, help="Base production MDP")
    parser.add_argument("--complex-mdp", help="Complex-specific base MDP")
    parser.add_argument("--solvent-mdp", help="Solvent-specific base MDP")
    parser.add_argument(
        "--complex-topology",
        required=True,
        help="Complex GROMACS topology",
    )
    parser.add_argument(
        "--complex-structure",
        required=True,
        help="Equilibrated complex structure",
    )
    parser.add_argument(
        "--solvent-topology",
        required=True,
        help="Solvated ligand GROMACS topology",
    )
    parser.add_argument(
        "--solvent-structure",
        required=True,
        help="Equilibrated solvated-ligand structure",
    )
    parser.add_argument("--moltype", required=True, help="Ligand molecule type")
    parser.add_argument("-o", "--outdir", default="abfe", help="Output directory")
    parser.add_argument(
        "--anchors",
        nargs=6,
        type=int,
        metavar=("R1", "R2", "R3", "L1", "L2", "L3"),
        help="One-based Boresch anchor atom indices",
    )
    parser.add_argument(
        "--receptor-selection",
        default="protein",
        help="MDTraj receptor selection for automatic anchors",
    )
    parser.add_argument(
        "--ligand-selection",
        help="MDTraj ligand selection for automatic anchors",
    )
    parser.add_argument(
        "--anchor-search-distance",
        type=float,
        default=0.5,
        help="Maximum receptor-ligand anchor distance [nm]",
    )
    parser.add_argument(
        "--distance-spring",
        type=float,
        default=4184.0,
        help="Distance restraint spring [kJ/mol/nm^2]",
    )
    parser.add_argument(
        "--angle-spring",
        type=float,
        default=41.84,
        help="Angle restraint spring [kJ/mol/rad^2]",
    )
    parser.add_argument(
        "--dihedral-spring",
        type=float,
        default=41.84,
        help="Dihedral restraint spring [kJ/mol/rad^2]",
    )
    parser.add_argument(
        "--temperature", type=float, default=300.0, help="Temperature [K]"
    )
    parser.add_argument(
        "--lj-pme-comb-rule",
        choices=["Geometric", "Lorentz-Berthelot"],
        default="Lorentz-Berthelot",
        help="LJ-PME combination rule",
    )
    parser.add_argument("--charge-windows", type=int, default=12)
    parser.add_argument("--vdw-windows", type=int, default=12)
    parser.add_argument("--restraint-windows", type=int, default=12)
    parser.add_argument("--nstdhdl", type=int, default=100)
    parser.add_argument("--complex-checkpoint")
    parser.add_argument("--solvent-checkpoint")
    parser.add_argument("--complex-reference")
    parser.add_argument("--solvent-reference")
    parser.add_argument("--complex-index")
    parser.add_argument("--solvent-index")
    parser.add_argument("--deffnm", default="abfe")
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--maxwarn", type=int, default=0)
    parser.add_argument(
        "--no-tpr",
        action="store_true",
        help="Generate MDP files without final TPR files",
    )
    parser.set_defaults(func=run)


def _existing(path, label):
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} not found: {resolved}")
    return resolved


def _setup_args(
    *,
    mdp,
    topology,
    structure,
    outdir,
    mode,
    moltype,
    windows,
    lambdas,
    nstdhdl,
    checkpoint,
    reference,
    index,
    args,
):
    return SimpleNamespace(
        mdp=str(mdp),
        topology=str(topology),
        structure=str(structure),
        outdir=str(outdir),
        mode=mode,
        moltype=moltype,
        windows=windows,
        coul_windows=windows,
        vdw_windows=windows,
        fep_lambdas=lambdas if mode == "transform" else None,
        coul_lambdas=lambdas if mode == "charge" else None,
        vdw_lambdas=lambdas if mode == "vdw" else None,
        calc_lambda_neighbors=1,
        nstdhdl=nstdhdl,
        checkpoint=str(checkpoint) if checkpoint else None,
        reference=str(reference) if reference else None,
        b_reference=None,
        index=str(index) if index else None,
        deffnm=args.deffnm,
        gmx=args.gmx,
        maxwarn=args.maxwarn,
        no_grompp=args.no_tpr,
    )


def _linspace(count):
    if count < 2:
        raise ValueError("ABFE window counts must be at least 2")
    return [index / (count - 1) for index in range(count)]


def run(args):
    if args.temperature <= 0:
        raise ValueError("--temperature must be positive")
    if args.nstdhdl <= 0:
        raise ValueError("--nstdhdl must be positive")
    if args.maxwarn < 0:
        raise ValueError("--maxwarn must be non-negative")
    if args.anchors is None and not args.ligand_selection:
        raise ValueError("--anchors or --ligand-selection is required")

    common_mdp = _existing(args.mdp, "Base MDP")
    complex_mdp = (
        _existing(args.complex_mdp, "Complex MDP") if args.complex_mdp else common_mdp
    )
    solvent_mdp = (
        _existing(args.solvent_mdp, "Solvent MDP") if args.solvent_mdp else common_mdp
    )
    complex_topology = _existing(args.complex_topology, "Complex topology")
    complex_structure = _existing(args.complex_structure, "Complex structure")
    solvent_topology = _existing(args.solvent_topology, "Solvent topology")
    solvent_structure = _existing(args.solvent_structure, "Solvent structure")
    complex_checkpoint = (
        _existing(args.complex_checkpoint, "Complex checkpoint")
        if args.complex_checkpoint
        else None
    )
    solvent_checkpoint = (
        _existing(args.solvent_checkpoint, "Solvent checkpoint")
        if args.solvent_checkpoint
        else None
    )
    complex_reference = (
        _existing(args.complex_reference, "Complex reference")
        if args.complex_reference
        else complex_structure
    )
    solvent_reference = (
        _existing(args.solvent_reference, "Solvent reference")
        if args.solvent_reference
        else solvent_structure
    )
    complex_index = (
        _existing(args.complex_index, "Complex index") if args.complex_index else None
    )
    solvent_index = (
        _existing(args.solvent_index, "Solvent index") if args.solvent_index else None
    )

    outdir = Path(args.outdir).expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {outdir}")
    input_dir = outdir / "inputs"
    input_dir.mkdir(parents=True)

    if args.anchors is not None:
        anchors = list(args.anchors)
        geometry = calculate_anchor_geometry(complex_structure, anchors)
    else:
        anchors, geometry = select_boresch_anchors(
            complex_structure,
            receptor_selection=args.receptor_selection,
            ligand_selection=args.ligand_selection,
            search_distance=args.anchor_search_distance,
        )
    anchor_index = input_dir / "complex_abfe.ndx"
    if complex_index is None:
        run_cmd(
            [
                args.gmx,
                "make_ndx",
                "-f",
                str(complex_structure),
                "-o",
                str(anchor_index),
            ],
            input="q\n",
            cwd=complex_topology.parent,
            log="Generated default complex index groups",
        )
        complex_index = anchor_index
    write_anchor_index(anchor_index, anchors, complex_index)

    common_overrides = {
        "vdwtype": "PME",
        "DispCorr": "no",
        "lj-pme-comb-rule": args.lj_pme_comb_rule,
    }
    restrained_settings = boresch_pull_settings(
        geometry,
        distance_spring=args.distance_spring,
        angle_spring=args.angle_spring,
        dihedral_spring=args.dihedral_spring,
        release=False,
    )
    release_settings = boresch_pull_settings(
        geometry,
        distance_spring=args.distance_spring,
        angle_spring=args.angle_spring,
        dihedral_spring=args.dihedral_spring,
        release=True,
    )
    complex_restrained_mdp = input_dir / "complex_restrained.mdp"
    complex_restrained_mdp.write_text(
        render_fep_mdp(
            complex_mdp.read_text(),
            {
                **common_overrides,
                **temperature_mdp_overrides(
                    complex_mdp.read_text(),
                    args.temperature,
                ),
                **restrained_settings,
            },
            remove_prefixes=("pull",),
        )
    )
    complex_release_mdp = input_dir / "complex_restraint_release.mdp"
    complex_release_mdp.write_text(
        render_fep_mdp(
            complex_mdp.read_text(),
            {
                **common_overrides,
                **temperature_mdp_overrides(
                    complex_mdp.read_text(),
                    args.temperature,
                ),
                **release_settings,
            },
            remove_prefixes=("pull",),
        )
    )
    solvent_output_mdp = input_dir / "solvent.mdp"
    solvent_output_mdp.write_text(
        render_fep_mdp(
            solvent_mdp.read_text(),
            {
                **common_overrides,
                **temperature_mdp_overrides(
                    solvent_mdp.read_text(),
                    args.temperature,
                ),
            },
            remove_prefixes=("pull",),
        )
    )

    charge_lambdas = _linspace(args.charge_windows)
    x = _linspace(args.vdw_windows)
    vdw_lambdas = [1.0 - (value - 1.0) ** 2 for value in x]
    restraint_lambdas = _linspace(args.restraint_windows)
    leg_definitions = {
        "solvent_charge": (
            solvent_output_mdp,
            solvent_topology,
            solvent_structure,
            "charge",
            charge_lambdas,
            solvent_checkpoint,
            solvent_reference,
            solvent_index,
        ),
        "solvent_vdw": (
            solvent_output_mdp,
            solvent_topology,
            solvent_structure,
            "vdw",
            vdw_lambdas,
            solvent_checkpoint,
            solvent_reference,
            solvent_index,
        ),
        "complex_charge": (
            complex_restrained_mdp,
            complex_topology,
            complex_structure,
            "charge",
            charge_lambdas,
            complex_checkpoint,
            complex_reference,
            anchor_index,
        ),
        "complex_vdw": (
            complex_restrained_mdp,
            complex_topology,
            complex_structure,
            "vdw",
            vdw_lambdas,
            complex_checkpoint,
            complex_reference,
            anchor_index,
        ),
        "complex_restraint": (
            complex_release_mdp,
            complex_topology,
            complex_structure,
            "transform",
            restraint_lambdas,
            complex_checkpoint,
            complex_reference,
            anchor_index,
        ),
    }
    legs = {}
    for name, (
        mdp,
        topology,
        structure,
        mode,
        lambdas,
        checkpoint,
        reference,
        index,
    ) in leg_definitions.items():
        setup_fep.run(
            _setup_args(
                mdp=mdp,
                topology=topology,
                structure=structure,
                outdir=outdir / name,
                mode=mode,
                moltype=None if mode == "transform" else args.moltype,
                windows=len(lambdas),
                lambdas=lambdas,
                nstdhdl=args.nstdhdl,
                checkpoint=checkpoint,
                reference=reference,
                index=index,
                args=args,
            )
        )
        legs[name] = name

    manifest = {
        "schema_version": ABFE_SCHEMA_VERSION,
        "workflow": "abfe",
        "temperature": args.temperature,
        "molecule_type": args.moltype,
        "anchors": anchors,
        "geometry": geometry,
        "springs": {
            "distance": args.distance_spring,
            "angle": args.angle_spring,
            "dihedral": args.dihedral_spring,
        },
        "long_range_method": "LJ-PME",
        "lj_pme_comb_rule": args.lj_pme_comb_rule,
        "gmx_executable": args.gmx,
        "legs": legs,
    }
    (outdir / ABFE_MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Generated ABFE cycle in %s", outdir)
