import argparse
from pathlib import Path
import subprocess
from ..config import *  # NOQA
from ..logger import generate_logger
from .calc_ion_conc import calc_ion_conc_from_volume

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx build_solution -f input_structure.pdb -o output_structure.pdb --ion_conc 0.15 --cation Na+ --anion Cl- --ligparam FRCMOD:LIB --boxsize 100 100 100
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
        default=Path(__file__).parent / "template_tleap.in",
        type=str,
        help="Template file for tleap",
    )

    parser.add_argument(
        "--keepfiles", action="store_true", help="Keep intermediate files"
    )


def run(args):
    # tleap
    lines = []
    with open(args.template_tleap) as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if "LOADPDB" in line:
                if args.input is not None:
                    line = line.replace(
                        "LOADPDB",
                        f"{SYSTEM_NAME} = loadpdb {args.input}",  # NOQA
                    )
                else:
                    LOGGER.warn("No input structure")
                    LOGGER.warn("System will be water system")
                    line = line.replace(
                        "LOADPDB",
                        f"{SYSTEM_NAME} = createunit '{SYSTEM_NAME}'",  # NOQA
                    )
            if "SYSTEM_NAME" in line:
                line = line.replace("SYSTEM_NAME", SYSTEM_NAME)  # NOQA
            if "OUT_DIR" in line:
                line = line.replace("OUT_DIR", args.outdir)
            if "BOX_SIZE" in line:
                line = line.replace("BOX_SIZE", " ".join(map(str, args.boxsize)))
            if "LIGAND_PARAMS" in line:
                if args.ligparam is not None:
                    frcmod = args.ligparam.split(":")[0]
                    lib = args.ligparam.split(":")[1]
                    cmd = f"""
loadamberparams {frcmod}
loadoff {lib}
                    """
                    line = line.replace("LIGAND_PARAMS", cmd)
                else:
                    line = ""
            if "WATER_MODEL" in line:
                if "solvatebox" in line:
                    water_name = args.water.upper()
                else:
                    water_name = args.water
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
                    """  # NOQA
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
    with open("tleap.in", "w") as f:
        f.write(cmd_tleap)
    cmd = "tleap -f tleap.in"
    subprocess.run(cmd, shell=True, check=True)

    LOGGER.info(
        f"{args.outdir}/leap.parm7 {args.outdir}/leap.rst7 {args.outdir}/leap.pdb generated"
    )

    if not args.keepfiles:
        cmd = "rm -f leap.log tleap.in"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info("leap.log tleap.in removed")
