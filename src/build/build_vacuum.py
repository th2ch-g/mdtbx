import argparse
import sys
from pathlib import Path

from ..config import SYSTEM_NAME
from ..logger import generate_logger
from ..utils.tleap import run_tleap

LOGGER = generate_logger(__name__)

TLEAP_TEMPLATE = """\
source leaprc.protein.ff14SB
source leaprc.water.tip3p
{addprecmd}
{ligand_params}
{SYSTEM_NAME} = loadpdb {input}
{addpostcmd}
saveamberparm {SYSTEM_NAME} {outdir}/leap.parm7 {outdir}/leap.rst7
savepdb {SYSTEM_NAME} {outdir}/leap.pdb
quit
"""


def add_subcmd(subparsers):
    """
    mdtbx build_vacuum -i input.pdb -o ./
    """
    parser = subparsers.add_parser(
        "build_vacuum",
        help="Build vacuum system (no water, no ions, no box)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Input PDB file",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="./",
        type=str,
        help="Output directory",
    )

    parser.add_argument(
        "--ligparam",
        type=str,
        help="Ligand parameter in FRCMOD:LIB format",
    )

    parser.add_argument(
        "--addprecmd",
        type=str,
        help="Additional tleap commands before loadpdb (e.g. source leaprc.GLYCAM_06j-1)",
    )

    parser.add_argument(
        "--addpostcmd",
        type=str,
        help="Additional tleap commands after loadpdb (e.g. SS-bond settings)",
    )

    parser.add_argument(
        "--keepfiles",
        action="store_true",
        help="Keep intermediate files (tleap.in, leap.log)",
    )

    parser.set_defaults(func=run)


def run(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    addprecmd = args.addprecmd if args.addprecmd else ""
    addpostcmd = args.addpostcmd if args.addpostcmd else ""

    ligand_params = ""
    if args.ligparam:
        parts = args.ligparam.split(":")
        if len(parts) != 2:
            LOGGER.error("--ligparam must be in FRCMOD:LIB format")
            sys.exit(1)
        frcmod, lib = parts
        ligand_params = f"loadamberparams {frcmod}\nloadoff {lib}"

    tleap_input = TLEAP_TEMPLATE.format(
        addprecmd=addprecmd,
        ligand_params=ligand_params,
        SYSTEM_NAME=SYSTEM_NAME,
        input=args.input,
        addpostcmd=addpostcmd,
        outdir=str(outdir),
    )

    run_tleap(tleap_input, keepfiles=args.keepfiles)

    LOGGER.info(
        f"{args.outdir}/leap.parm7 {args.outdir}/leap.rst7 {args.outdir}/leap.pdb generated"
    )
