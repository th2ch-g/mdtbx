import argparse
import json
from pathlib import Path
import sys
from ..config import SYSTEM_NAME
from ..logger import generate_logger
from ..utils.tleap import run_tleap
from .calc_ion_conc import calc_ion_conc_from_volume
from .packmol_water import repack_water

LOGGER = generate_logger(__name__)


def _center_pdb_in_box(input_path, output_path, boxsize):
    lines = Path(input_path).read_text().splitlines()
    coordinates = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            coordinates.append(
                tuple(float(line[start : start + 8]) for start in (30, 38, 46))
            )
    if not coordinates:
        raise ValueError(f"No PDB coordinates found in {input_path}")

    minima = [min(point[axis] for point in coordinates) for axis in range(3)]
    maxima = [max(point[axis] for point in coordinates) for axis in range(3)]
    spans = [maximum - minimum for minimum, maximum in zip(minima, maxima)]
    if any(span >= edge for span, edge in zip(spans, boxsize)):
        raise ValueError("Solute dimensions must be smaller than --boxsize")
    shifts = [
        edge / 2 - (minimum + maximum) / 2
        for minimum, maximum, edge in zip(minima, maxima, boxsize)
    ]

    output = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            xyz = [float(line[start : start + 8]) for start in (30, 38, 46)]
            moved = [value + shift for value, shift in zip(xyz, shifts)]
            line = (
                f"{line[:30]}{moved[0]:8.3f}{moved[1]:8.3f}{moved[2]:8.3f}{line[54:]}"
            )
        output.append(line)
    Path(output_path).write_text("\n".join(output) + "\n")


def add_subcmd(subparsers):
    """
    mdtbx build_solution -i input_structure.pdb -o ./outdir --ion_conc 0.15 --cation Na+ --anion Cl- --ligparam FRCMOD:LIB --boxsize 100 100 100
    """
    parser = subparsers.add_parser(
        "build_solution",
        help="Build solution system",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Input Structure file (.pdb)",
    )

    parser.add_argument(
        "-o", "--outdir", default="./", type=str, help="Output file path"
    )

    parser.add_argument(
        "--water",
        default="tip3p",
        type=str,
        help="Water model",
        choices=["tip3p", "opc"],
    )

    parser.add_argument(
        "--ion_conc", default=0.15, type=float, help="Ion concentration [M]"
    )

    parser.add_argument("--cation", default="Na+", type=str, help="Cation name")

    parser.add_argument("--anion", default="Cl-", type=str, help="Anion name")

    parser.add_argument("--noions", action="store_true", help="No ions")

    parser.add_argument(
        "--ligparam", type=str, help="Ligand parameter. e.g. --ligparam FRCMOD:LIB"
    )

    parser.add_argument(
        "--boxsize",
        nargs=3,
        type=float,
        default=[100, 100, 100],
        help="Box size [angstrom, angstrom, angstrom]",
    )

    parser.add_argument(
        "--addprecmd",
        type=str,
        help="Additional pre command before load structure in tleap (e.g. bond SYS.1.SG SYS.2.SG)",
    )

    parser.add_argument(
        "--addpostcmd",
        type=str,
        help="Additional command after load structure in tleap (e.g. bond SYS.1.SG SYS.2.SG)",
    )

    parser.add_argument(
        "--template-tleap",
        default=str(Path(__file__).parent.parent / "utils" / "template_tleap.in"),
        type=str,
        help="Template file for tleap",
    )

    parser.add_argument(
        "--keepfiles", action="store_true", help="Keep intermediate files"
    )

    parser.add_argument(
        "--water-seed",
        default=2026,
        type=int,
        help="Packmol random seed for water placement",
    )

    parser.add_argument(
        "--packmol-tolerance",
        default=2.0,
        type=float,
        help="Packmol intermolecular distance tolerance [angstrom]",
    )

    parser.set_defaults(func=run)


def run(args):
    if args.packmol_tolerance <= 0:
        raise ValueError("--packmol-tolerance must be positive")
    if args.ion_conc < 0:
        raise ValueError("--ion_conc must be non-negative")
    if args.water_seed < 0:
        raise ValueError("--water-seed must be non-negative")
    if any(edge <= args.packmol_tolerance for edge in args.boxsize):
        raise ValueError("Each --boxsize edge must exceed --packmol-tolerance")

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    input_path = (
        Path(args.input).expanduser().resolve() if args.input is not None else None
    )
    tleap_input_path = input_path
    centered_input = None
    if input_path is not None:
        centered_input = outdir / "tleap_centered.pdb"
        if centered_input.exists():
            raise FileExistsError(f"Refusing to overwrite {centered_input}")
        _center_pdb_in_box(input_path, centered_input, args.boxsize)
        tleap_input_path = centered_input

    # tleap
    lines = []
    with open(args.template_tleap) as f:
        for line in f:
            line = line.rstrip()
            if "LOADPDB" in line:
                if tleap_input_path is not None:
                    line = line.replace(
                        "LOADPDB",
                        f"{SYSTEM_NAME} = loadpdb {tleap_input_path}",
                    )
                else:
                    LOGGER.warning("No input structure")
                    LOGGER.warning("System will be water system")
                    line = line.replace(
                        "LOADPDB",
                        f"{SYSTEM_NAME} = createunit '{SYSTEM_NAME}'",
                    )
            if "SYSTEM_NAME" in line:
                line = line.replace("SYSTEM_NAME", SYSTEM_NAME)
            if "OUT_DIR" in line:
                line = line.replace("OUT_DIR", str(outdir))
            if "BOX_SIZE" in line:
                line = line.replace("BOX_SIZE", " ".join(map(str, args.boxsize)))
            if "LIGAND_PARAMS" in line:
                if args.ligparam is not None:
                    parts = args.ligparam.split(":")
                    if len(parts) != 2:
                        LOGGER.error("--ligparam must be in FRCMOD:LIB format")
                        sys.exit(1)
                    frcmod, lib = (
                        str(Path(value).expanduser().resolve()) for value in parts
                    )
                    cmd = f"""
loadamberparams {frcmod}
loadoff {lib}
                    """
                    line = line.replace("LIGAND_PARAMS", cmd)
                else:
                    line = ""
            if "WATER_MODEL" in line:
                water_name = args.water.upper() if "solvatebox" in line else args.water
                line = line.replace("WATER_MODEL", water_name)
            if "ADDION" in line:
                if args.noions:
                    LOGGER.info("Ions will not be added")
                    line = ""
                else:
                    ion_num = calc_ion_conc_from_volume(
                        args.boxsize[0] * args.boxsize[1] * args.boxsize[2],
                        args.ion_conc,
                    )  # cubic method
                    cmd = f"""
addionsrand {SYSTEM_NAME} {args.cation} {ion_num}
addionsrand {SYSTEM_NAME} {args.anion} 0
                    """
                    line = line.replace("ADDION", cmd)
            if "ADDPRECMD" in line:
                if args.addprecmd is not None:
                    line = line.replace("ADDPRECMD", args.addprecmd)
                else:
                    line = ""
            if "ADDPOSTCMD" in line:
                if args.addpostcmd is not None:
                    line = line.replace("ADDPOSTCMD", args.addpostcmd)
                else:
                    line = ""
            lines.append(line)

    cmd_tleap = "\n".join(lines)
    run_tleap(
        cmd_tleap,
        keepfiles=args.keepfiles,
        extra_cleanup=(centered_input,) if centered_input is not None else (),
        cwd=outdir,
    )

    packmol_report = repack_water(
        outdir / "leap.parm7",
        outdir / "leap.rst7",
        outdir / "leap.pdb",
        seed=args.water_seed,
        tolerance=args.packmol_tolerance,
    )
    if packmol_report is not None:
        (outdir / "packmol_transfer.json").write_text(
            json.dumps(packmol_report, indent=2) + "\n"
        )

    LOGGER.info(
        f"{args.outdir}/leap.parm7 {args.outdir}/leap.rst7 {args.outdir}/leap.pdb generated"
    )
