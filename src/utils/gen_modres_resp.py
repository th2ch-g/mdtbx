import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx gen_modres_resp
    """
    parser = subparsers.add_parser(
        "gen_modres_resp",
        help="Centering modified residue parameters with Gaussian",
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
        "--postatomname",
        type=str,
        help="Post atom name for connection",
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


def run(args):
    # ref: https://qiita.com/tacoma/items/02474d9aaa99b903e4ee
    # hint: Check Structure as PDB before running. You can check by antechamber or obabel
    # hint: Need cap atom for specifying SEP_BOND

    filetype = Path(args.structure).suffix[1:]
    cmd = f"obabel -i {filetype} {args.structure} -o gjf > structure_optimization.gjf"
    subprocess.run(cmd, shell=True, check=True)

    with open("structure_optimization.gjf") as ref:
        lines = ref.readlines()

    lines[0] = "%chk=structure_optimization.chk\n"
    lines[1] = f"%mem={args.memory}GB\n"
    lines.insert(2, f"%nprocshared={args.threads}\n")
    lines.insert(3, f"{STRUCTURE_OPTIMIZATION}\n")  # NOQA

    for idx, line in enumerate(lines):
        line = line.strip()
        if len(line.split()) == 2:
            try:
                _charge = int(line.split()[0])
                _multiplicity = int(line.split()[1])
                target_idx = idx
                lines[target_idx] = f"{args.charge} {args.multiplicity}\n"
                break
            except Exception:
                continue

    with open("structure_optimization.gjf", "w") as f:
        f.writelines(lines)

    cmd = f"{GAUSSIAN_CMD} < structure_optimization.gjf > structure_optimization.log"  # NOQA
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("structure_optimization.log generated")

    # single point
    cmd = f"obabel -i {GAUSSIAN_CMD} structure_optimization.log -o gjf > single_point_calculation.gjf"  # NOQA
    subprocess.run(cmd, shell=True, check=True)

    # # single point calculation
    # filetype = Path(args.structure).suffix[1:]
    # cmd = f"obabel -i {filetype} {args.structure} -o gjf > single_point_calculation.gjf"
    # subprocess.run(cmd, shell=True, check=True)

    with open("single_point_calculation.gjf") as ref:
        lines = ref.readlines()

    lines[0] = "%chk=single_point_calculation.chk\n"
    lines[1] = f"%mem={args.memory}GB\n"
    lines.insert(2, f"%nprocshared={args.threads}\n")
    lines.insert(3, f"{SINGLE_POINT_CALCULATION}\n")  # NOQA

    for idx, line in enumerate(lines):
        line = line.strip()
        if len(line.split()) == 2:
            try:
                _charge = int(line.split()[0])
                _multiplicity = int(line.split()[1])
                target_idx = idx
                lines[target_idx] = f"{args.charge} {args.multiplicity}\n"
                break
            except Exception:
                continue

    with open("single_point_calculation.gjf", "w") as f:
        f.writelines(lines)

    cmd = (
        f"{GAUSSIAN_CMD} < single_point_calculation.gjf > single_point_calculation.log"  # NOQA
    )
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("single_point_calculation.log generated")

    cmd = f"antechamber -fi gout -i single_point_calculation.log -fo ac -o {args.resname}.ac -pf y -rn {args.resname} -at amber -s 2 -nc {args.charge} -m {args.multiplicity}"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.ac generated")

    cmd = f"espgen -i single_point_calculation.log -o {args.resname}.esp"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.esp generated")

    atom_charges = "\n".join(
        [
            f"{atom_charge.split(':')[0]} {atom_charge.split(':')[1]}"
            for atom_charge in args.atomcharge
        ]
    )
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

    cmd = f"residuegen {args.resname}.in"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}.prep generated")

    cmd = (
        f"parmchk2 -i {args.resname}.prep -f prepi -o {args.resname}_1.frcmod -s parm10"
    )
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

    cmd = (
        f"parmchk2 -i {args.resname}.prep -f prepi -o {args.resname}_2.frcmod -s gaff2"
    )
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.resname}_2.frcmod generated")

    # if args.postatomname is not None:
    #     section = None
    #     with open(f"{args.resname}.prep") as f:
    #         for idx, line in enumerate(f):
    #             line = line.rstrip()
    #             if line.startwith("CORRECT"):
    #                 section = "CORRECT"
    #                 continue
    #             if args.resname in line:

    LOGGER.warning(f"You need to modify {args.resname}.prep manually.")
    LOGGER.warning("Like ' 18  C8    C     E' => ' 18  C8    C     M'")

    print(
        f"use --addprecmd 'loadAmberPrep {args.resname}.prep\nloadamberparams {args.resname}_1.frcmod\nloadamberparams {args.resname}_2.frcmod'"
    )
