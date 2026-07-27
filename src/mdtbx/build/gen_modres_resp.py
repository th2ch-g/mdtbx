import argparse
from pathlib import Path

from ..config import GAUSSIAN_CMD, SINGLE_POINT_CALCULATION
from ..logger import generate_logger
from ..utils.proc import run_cmd
from .gaussian import (
    configure_gaussian_input,
    convert_to_gaussian_input,
    run_gaussian,
)

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_modres_resp
    """
    parser = subparsers.add_parser(
        "gen_modres_resp",
        help="Generate modified residue parameters with Gaussian",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Compound(Structure) file (.mol or .mol2). ",
    )

    parser.add_argument(
        "-r",
        "--resname",
        default="UNK",
        type=str,
        help="Residue name",
    )

    parser.add_argument(
        "--sepbond1",
        type=str,
        nargs=2,
        # default=["N", "C"],
        help="Bond1",
    )

    parser.add_argument(
        "--sepbond2",
        type=str,
        nargs=2,
        # default=["C", "N"],
        help="Bond2",
    )

    parser.add_argument(
        "--atomcharge",
        type=str,
        nargs="*",
        default=[],
        # default=["N:-0.4175", "HA:0.2719", "C:0.5973", "O:-0.5679"],
        help="Fixed atom charge",
    )

    parser.add_argument(
        "-m", "--multiplicity", default=1, type=int, help="Multiplicity"
    )

    parser.add_argument("-c", "--charge", default=0, type=int, help="Charge")

    parser.add_argument(
        "--memory", default="60", type=int, help="Memory(GB) for Gaussian"
    )

    parser.add_argument(
        "--threads", default=16, type=int, help="Number of threads for Gaussian"
    )

    parser.set_defaults(func=run)


def run(args):
    # ref: https://qiita.com/tacoma/items/02474d9aaa99b903e4ee
    # hint: Check Structure as PDB before running. You can check by antechamber or obabel
    # hint: Need cap atom for specifying SEP_BOND

    filetype = Path(args.structure).suffix[1:]
    convert_to_gaussian_input(
        args.structure,
        filetype,
        "single_point_calculation.gjf",
    )
    configure_gaussian_input(
        "single_point_calculation.gjf",
        checkpoint="single_point_calculation.chk",
        memory_gb=args.memory,
        threads=args.threads,
        route=SINGLE_POINT_CALCULATION,
        charge=args.charge,
        multiplicity=args.multiplicity,
    )
    run_gaussian(
        GAUSSIAN_CMD,
        "single_point_calculation.gjf",
        "single_point_calculation.log",
    )
    LOGGER.info("single_point_calculation.log generated")

    cmd = [
        "antechamber",
        "-fi",
        "gout",
        "-i",
        "single_point_calculation.log",
        "-fo",
        "ac",
        "-o",
        f"{args.resname}.ac",
        "-pf",
        "y",
        "-rn",
        args.resname,
        "-at",
        "amber",
        "-s",
        "2",
        "-nc",
        str(args.charge),
        "-m",
        str(args.multiplicity),
    ]
    run_cmd(cmd, log=f"{args.resname}.ac generated")

    cmd = [
        "espgen",
        "-i",
        "single_point_calculation.log",
        "-o",
        f"{args.resname}.esp",
    ]
    run_cmd(cmd, log=f"{args.resname}.esp generated")

    atom_charges = []
    for atom_charge in args.atomcharge:
        atom_name, separator, charge = atom_charge.partition(":")
        if not separator or not atom_name or not charge:
            raise ValueError("--atomcharge values must use ATOM:CHARGE format")
        atom_charges.append(f"{atom_name} {charge}")
    atom_charges = "\n".join(atom_charges)
    if args.sepbond1 is not None:
        sep_bond1 = f"SEP_BOND {args.sepbond1[0]} {args.sepbond1[1]}"
    else:
        sep_bond1 = ""

    if args.sepbond2 is not None:
        sep_bond2 = f"SEP_BOND {args.sepbond2[0]} {args.sepbond2[1]}"
    else:
        sep_bond2 = ""

    content = f"""
INPUT_FILE {args.resname}.ac
CONF_NUM 1
ESP_FILE {args.resname}.esp
{sep_bond1}
{sep_bond2}
{atom_charges}
NET_CHARGE {args.charge}
PREP_FILE: {args.resname}.prep
RESIDUE_FILE_NAME: {args.resname}.res
RESIDUE_SYMBOL: {args.resname}
    """

    with open(f"{args.resname}.in", "w") as f:
        f.write(content)

    cmd = ["residuegen", f"{args.resname}.in"]
    run_cmd(cmd, log=f"{args.resname}.prep generated")

    cmd = [
        "parmchk2",
        "-i",
        f"{args.resname}.prep",
        "-f",
        "prepi",
        "-o",
        f"{args.resname}_parm10.frcmod",
        "-s",
        "parm10",
    ]
    run_cmd(cmd, log=f"{args.resname}_parm10.frcmod generated")

    lines = []
    with open(f"{args.resname}_parm10.frcmod") as f:
        for line in f:
            line = line.rstrip()
            if "ATTN" not in line:
                lines.append(line)
    content = "\n".join(lines)
    with open(f"{args.resname}_parm10.frcmod", "w") as f:
        f.write(content + "\n")
    LOGGER.info(f"{args.resname}_parm10.frcmod updated")

    cmd = [
        "parmchk2",
        "-i",
        f"{args.resname}.prep",
        "-f",
        "prepi",
        "-o",
        f"{args.resname}_gaff2.frcmod",
        "-s",
        "gaff2",
    ]
    run_cmd(cmd, log=f"{args.resname}_gaff2.frcmod generated")

    LOGGER.warning(f"You need to modify {args.resname}.prep manually.")
    LOGGER.warning("Like ' 18  C8    C     E' => ' 18  C8    C     M'")

    """
    mol2 file should be as follows
    NH(mainchain)
    CA, HA, CB, HB1, HB2,...(sidechain)
    CO(mainchain)

    --addprecmd '
    loadAmberPrep {args.resname}.prep
    loadamberparams {args.resname}_gaff2.frcmod
    loadamberparams {args.resname}_parm10.frcmod'
    """
