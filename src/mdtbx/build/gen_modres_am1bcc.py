import argparse
from pathlib import Path

from ..logger import generate_logger
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_modres_am1bcc
    """
    parser = subparsers.add_parser(
        "gen_modres_am1bcc",
        help="Generate modified residue parameters with AM1BCC (WIP)",
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

    parser.add_argument("--headname", required=True, type=str, help="Atom name of Head")

    parser.add_argument("--tailname", required=True, type=str, help="Atom name of Tail")

    parser.add_argument(
        "--mainchain", nargs="*", type=str, help="Atom name of Mainchain"
    )

    parser.add_argument("--omitname", nargs="*", type=str, help="Atom name of Omit")

    parser.add_argument(
        "--preheadtype", default="C", type=str, help="Atom name of Prehead"
    )

    parser.add_argument(
        "--posttailtype", default="N", type=str, help="Atom name of Posttail"
    )

    parser.set_defaults(func=run)


def run(args):
    # ref: https://ambermd.org/tutorials/basic/tutorial5/index.php
    filetype = Path(args.structure).suffix[1:]
    # antechamber uses "mdl" as the format flag for MDL .mol files
    if filetype == "mol":
        filetype = "mdl"
    cmd = [
        "antechamber",
        "-fi",
        filetype,
        "-i",
        args.structure,
        "-bk",
        args.resname,
        "-fo",
        "ac",
        "-o",
        f"{args.resname}.ac",
        "-c",
        "bcc",
        "-at",
        "amber",
        "-pf",
        "y",
        "-s",
        "2",
        "-nc",
        str(args.charge),
        "-m",
        str(args.multiplicity),
    ]
    run_cmd(cmd, log=f"{args.resname}.ac generated")

    if args.mainchain is not None:
        main_chain = "\n".join([f"MAIN_CHAIN {mc}" for mc in args.mainchain])
    else:
        main_chain = ""

    if args.omitname is not None:
        omit_name = "\n".join([f"OMIT_NAME {om}" for om in args.omitname])
    else:
        omit_name = ""

    mc = f"""
HEAD_NAME {args.headname}
TAIL_NAME {args.tailname}
{main_chain}
{omit_name}
PRE_HEAD_TYPE {args.preheadtype}
POST_TAIL_TYPE {args.posttailtype}
CHARGE {args.charge}
    """

    with open(f"{args.resname}.mc", "w") as f:
        f.write(mc)

    cmd = [
        "prepgen",
        "-i",
        f"{args.resname}.ac",
        "-o",
        f"{args.resname}.prepin",
        "-m",
        f"{args.resname}.mc",
        "-rn",
        args.resname,
    ]
    run_cmd(cmd, log=f"{args.resname}.prepin generated")

    parm10_path = (
        Path(__file__).parent.parent.parent
        / ".pixi/envs/default/dat/leap/parm/parm10.dat"
    )
    cmd = [
        "parmchk2",
        "-i",
        f"{args.resname}.prepin",
        "-f",
        "prepi",
        "-o",
        f"{args.resname}_1.frcmod",
        "-a",
        "Y",
        "-p",
        str(parm10_path),
    ]
    run_cmd(cmd, log=f"{args.resname}_1.frcmod generated")

    lines = []
    with open(f"{args.resname}_1.frcmod") as f:
        for line in f:
            line = line.rstrip()
            if "ATTN" not in line:
                lines.append(line)
    content = "\n".join(lines)
    with open(f"{args.resname}_1.frcmod", "w") as f:
        f.write(content + "\n")
    LOGGER.info(f"{args.resname}_1.frcmod updated")

    cmd = [
        "parmchk2",
        "-i",
        f"{args.resname}.prepin",
        "-f",
        "prepi",
        "-o",
        f"{args.resname}_2.frcmod",
        "-s",
        "gaff2",
    ]
    run_cmd(cmd, log=f"{args.resname}_2.frcmod generated")
