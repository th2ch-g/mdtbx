import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_modres_am1bcc
    """
    parser = subparsers.add_parser(
        "gen_modres_am1bcc",
        help="Centering modified residue parameters with AM1BCC",
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


def run(args):
    # ref: https://ambermd.org/tutorials/basic/tutorial5/index.php
    filetype = Path(args.structure).suffix[1:]
    cmd = f"antechamber -fi {filetype} -i {args.structure} -bk {args.resname} -fo ac -o {args.resname}.ac -c bcc -at amber -pf y -s 2 -nc {args.charge} -m {args.multiplicity}"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.ac generated")

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

    cmd = f"prepgen -i {args.resname}.ac -o {args.resname}.prepin -m {args.resname}.mc -rn {args.resname}"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.prepin generated")

    cmd = f"parmchk2 -i {args.resname}.prepin -f prepi -o {args.resname}_1.frcmod -a Y -p {Path(__file__).parent.parent.parent}/.pixi/envs/default/dat/leap/parm/parm10.dat"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}_1.frcmod generated")

    lines = []
    with open(f"{args.resname}_1.frcmod") as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if "ATTN" not in line:
                lines.append(line)
    content = "\n".join(lines)
    with open(f"{args.resname}_1.frcmod", "w") as f:
        f.writelines(content)
    LOGGER.info(f"{args.resname}_1.frcmod updated")

    cmd = f"parmchk2 -i {args.resname}.prepin -f prepi -o {args.resname}_2.frcmod"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}_2.frcmod generated")
