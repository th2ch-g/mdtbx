import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx convert -p leap.parm7 -x leap.rst7
    """
    parser = subparsers.add_parser(
        "convert",
        help="Convert files from Amber to Gromacs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p", "--parm", required=True, type=str, help="Parameter file (.parm7, .prmtop)"
    )

    parser.add_argument(
        "-x", "--rst", required=True, type=str, help="RST file (.rst7, .inpcrd)"
    )

    parser.add_argument(
        "--type",
        default="parmed",
        type=str,
        help="Type of conversion",
        choices=["parmed", "acpype"],
    )


def run(args):
    if args.type == "parmed":
        cmd = f"amb2gro_top_gro.py -p {args.parm} -c {args.rst} -t gmx.top -g gmx.gro -b gmx.pdb"
        subprocess.run(cmd, shell=True)
        LOGGER.info("gmx.gro generated")
        LOGGER.info("gmx.top generated")
        LOGGER.info("gmx.pdb generated")
    elif args.type == "acpype":
        cmd = f"acpype -p {args.parm} -x {args.rst}"
        subprocess.run(cmd, shell=True)
        stem = Path(args.parm).stem
        cmd = f"cp {stem}.amb2gmx/{stem}_GMX.gro gmx.gro"
        subprocess.run(cmd, shell=True)
        cmd = f"cp {stem}.amb2gmx/{stem}_GMX.top gmx.top"
        subprocess.run(cmd, shell=True)
        LOGGER.info("gmx.gro generated")
        LOGGER.info("gmx.top generated")
        cmd = f"rm -rf {stem}.amb2gmx/"
        subprocess.run(cmd, shell=True)
        LOGGER.info(f"{stem}.amb2gmx/ removed")
