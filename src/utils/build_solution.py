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
        "-f",
        "--input",
        required=True,
        type=str,
        help="Input Structure file (.pdb)",
    )

    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Output file path"
    )

    parser.add_argument(
        "--ion_conc", default=0.15, type=float, help="Ion concentration [M]"
    )

    parser.add_argument("--cation", default="Na+", type=str, help="Cation name")

    parser.add_argument("--anion", default="Cl-", type=str, help="Anion name")

    parser.add_argument(
        "--ligparam", type=str, help="Ligand parameter. e.g. --ligparam FRCMOD:LIB"
    )

    parser.add_argument(
        "--boxsize", nargs=3, type=float, help="Box size [angstrom, angstrom, angstrom]"
    )
    parser.add_argument(
        "--addcmd",
        type=str,
        help="Additional command in tleap (e.g. bond SYS.1.SG SYS.2.SG)",
    )


def run(args):
    # tleap
    # SYSTEM_NAME, ION_NUM, INPUT_PDB, LIGAND_PARAMS, SSBONDS
    lines = []
    with open(Path(__file__).parent / "template_tleap.in") as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if "INPUT_PDB" in line:
                line = line.replace("INPUT_PDB", args.input)
            if "SYSTEM_NAME" in line:
                line = line.replace("SYSTEM_NAME", SYSTEM_NAME)  # NOQA
            if "OUT_DIR" in line:
                line = line.replace("OUT_DIR", args.output)
            if "CATION" in line:
                line = line.replace("CATION", args.cation)
            if "ANION" in line:
                line = line.replace("ANION", args.anion)
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
            if "ION_NUM" in line:
                ionnum = calc_ion_conc_from_volume(
                    args.boxsize[0] * args.boxsize[1] * args.boxsize[2], args.ion_conc
                )  # cubic method
                line = line.replace("ION_NUM", str(ionnum))
            if "ADDCMD" in line:
                if args.addcmd is not None:
                    line = line.replace("ADDCMD", args.addcmd)
                else:
                    line = ""
            lines.append(line)

    cmd_tleap = "\n".join(lines)
    with open("tleap.in", "w") as f:
        f.write(cmd_tleap)
    cmd = "tleap -f tleap.in"
    subprocess.run(cmd, shell=True, check=True)

    LOGGER.info(
        f"{args.output}/leap.parm7 {args.output}/leap.rst7 {args.output}/leap.pdb generated"
    )

    cmd = "rm -f leap.log tleap.in"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("leap.log tleap.in removed")
