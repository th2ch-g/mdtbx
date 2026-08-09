import argparse
import shutil
import tempfile
from pathlib import Path

from ..utils.gmx import gmx_tempfile
from ..utils.proc import run_cmd
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx amb2gro -p leap.parm7 -x leap.rst7
    """
    parser = subparsers.add_parser(
        "amb2gro",
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

    parser.add_argument(
        "--no-editconf", action="store_true", help="Do not run gmx editconf"
    )

    parser.set_defaults(func=run)


def run(args):
    if args.type == "parmed":
        cmd = [
            "amb2gro_top_gro.py",
            "-p",
            args.parm,
            "-c",
            args.rst,
            "-t",
            "gmx.top",
            "-g",
            "gmx.gro",
            "-b",
            "gmx.pdb",
        ]
        run_cmd(cmd, log="gmx.gro generated")
        LOGGER.info("gmx.top generated")
        LOGGER.info("gmx.pdb generated")
    # acpype may create wrong gro file
    # because of overflow of residue numbers or atom numbers
    elif args.type == "acpype":
        stem = Path(args.parm).stem
        with tempfile.TemporaryDirectory(prefix="mdtbx-acpype-") as tempdir:
            cmd = [
                "acpype",
                "-p",
                str(Path(args.parm).resolve()),
                "-x",
                str(Path(args.rst).resolve()),
            ]
            run_cmd(cmd, cwd=tempdir)
            generated_dir = Path(tempdir) / f"{stem}.amb2gmx"
            shutil.copy2(generated_dir / f"{stem}_GMX.gro", "gmx.gro")
            shutil.copy2(generated_dir / f"{stem}_GMX.top", "gmx.top")
        LOGGER.info("gmx.gro generated")
        LOGGER.info("gmx.top generated")

    if not args.no_editconf:
        with gmx_tempfile(".gro") as edited_path:
            cmd = [
                "gmx",
                "editconf",
                "-f",
                "gmx.gro",
                "-o",
                edited_path,
                "-resnr",
                "1",
            ]
            run_cmd(cmd, log="gmx editconf residue renumbering completed")
            Path(edited_path).replace("gmx.gro")
