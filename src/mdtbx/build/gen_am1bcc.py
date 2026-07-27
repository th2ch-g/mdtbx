import argparse
from pathlib import Path

from ..logger import generate_logger
from ..utils.proc import run_cmd
from ..utils.tleap import run_tleap

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_am1bcc -s structure.mol -m multiplicity -c charge
    """
    parser = subparsers.add_parser(
        "gen_am1bcc",
        help="Generate AM1BCC charges for Ligand",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Compound(Structure) file (.mol or .mol2)",
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

    parser.set_defaults(func=run)


def run(args):
    filetype = Path(args.structure).suffix[1:]
    # antechamber uses "mdl" as the format flag for MDL .mol files
    if filetype == "mol":
        filetype = "mdl"
    cmd = [
        "antechamber",
        "-i",
        args.structure,
        "-fi",
        filetype,
        "-o",
        f"{args.resname}.mol2",
        "-fo",
        "mol2",
        "-c",
        "bcc",
        "-s",
        "2",
        "-nc",
        str(args.charge),
        "-m",
        str(args.multiplicity),
        "-rn",
        args.resname,
        "-pf",
        "y",
    ]
    run_cmd(cmd, log=f"{args.resname}.mol2 generated")

    cmd = [
        "parmchk2",
        "-i",
        f"{args.resname}.mol2",
        "-f",
        "mol2",
        "-o",
        f"{args.resname}.frcmod",
        "-s",
        "gaff2",
    ]
    run_cmd(cmd, log=f"{args.resname}.frcmod generated")

    cmd_tleap = f"""
source leaprc.gaff2
loadamberparams {args.resname}.frcmod
{args.resname} = loadmol2 {args.resname}.mol2
saveoff {args.resname} {args.resname}.lib
quit
    """
    run_tleap(cmd_tleap)
    LOGGER.info(f"{args.resname}.lib generated")
