import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_am1bcc -s structure.mol -m multiplicity -c charge
    """
    parser = subparsers.add_parser(
        "gen_am1bcc",
        help="Generate AM1BCC charges",
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
        default="UNK",
        type=str,
        help="Residue name",
    )

    parser.add_argument(
        "-m", "--multiplicity", default=1, type=int, help="Multiplicity"
    )

    parser.add_argument("-c", "--charge", default=0, type=int, help="Charge")


def run(args):
    filetype = Path(args.structure).suffix[1:]
    cmd = f"antechamber -i {args.structure} -fi {filetype} -o {args.resname}.mol2 -fo mol2 -c bcc -s 2 -nc {args.charge} -m {args.multiplicity} -rn {args.resname} -pf y"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.mol2 generated")

    cmd = f"parmchk2 -i {args.output_prefix}.mol2 -f mol2 -o {args.resname}.frcmod -s gaff2"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.frcmod generated")

    cmd_tleap = """
source leaprc.gaff2
loadamberparams {args.resname}.frcmod
{args.resname} = loadmol2 {args.resname}.mol2
saveoff {args.resname} {args.resname}.lib
quit
    """
    cmd = f"echo '{cmd_tleap}' | tleap -f -"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.lib generated")

    cmd = "rm -f leap.log"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("leap.log removed")
