import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_resp -s structure.mol -m multiplicity -c charge
    """
    parser = subparsers.add_parser(
        "gen_resp",
        help="Generate RESP charges",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Compound(Structure) file (.mol or .mol2 or .pdb?)",
    )

    parser.add_argument(
        "-r",
        "--resname",
        required=True,
        type=str,
        help="Residue name",
    )

    parser.add_argument(
        "-m", "--multiplicity", required=True, type=int, help="Multiplicity"
    )

    parser.add_argument("-c", "--charge", required=True, type=int, help="Charge")

    parser.add_argument(
        "--memory", default="60GB", type=str, help="Memory for Gaussian"
    )

    parser.add_argument(
        "--threads", default=16, type=int, help="Number of threads for Gaussian"
    )


def run(args):
    # structure optimization
    cmd = f"obabel -i {filetype} {args.structure} -o gjf > structure_optimization.gjf"
    subprocess.run(cmd, shell=True, check=True)

    with open("structure_optimization.gjf") as ref:
        lines = ref.readlines()

    lines[0] = f"%chk=structure_optimization.chk"
    lines[1] = f"%mem={args.memory}GB"
    lines.insert(2, f"%nprocshared={args.threads}")
    lines.insert(3, STRUCTURE_OPTIMIZATION)

    with open("structure_optimization.gjf", "w") as f:
        f.writelines(lines)

    cmd = f"g16 < structure_optimization.gjf > structure_optimization.log"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("structure_optimization.log generated")

    # single point
    cmd = f"obabel -i g16 structure_optimization.log -o gjf > single_point_calculation.gjf"
    subprocess.run(cmd, shell=True, check=True)

    with open(f"single_point_calculation.gjf") as ref:
        lines = ref.readlines()

    lines[0] = f"%chk=single_point_calculation.chk"
    lines[1] = f"%mem={args.memory}GB"
    lines.insert(2, f"%nprocshared={args.threads}")
    lines.insert(3, single_point_calculation)

    with open("single_point_calculation.gjf", "w") as f:
        f.writelines(lines)

    cmd = f"g16 < single_point_calculation.gjf > single_point_calculation.log"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("single_point_calculation.log generated")

    # RESP
    cmd = f"antechamber -i single_point_calculation.log -fi gout -o {args.resname}.mol2 -fo mol2 -c resp -at gaff2 -nc {args.charge} -m {args.multiplicity} -rn {args.resname} -pf y"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.mol2 generated")

    cmd = "rm -f QOUT esout punch qout"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("QOUT esout punch qout removed")
