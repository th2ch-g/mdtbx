import argparse
import json
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

    parser.add_argument(
        "-o", "--outdir", default=".", type=str, help="Output directory"
    )

    parser.set_defaults(func=run)


def run(args):
    parm = Path(args.parm).expanduser().resolve()
    rst = Path(args.rst).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    gro = outdir / "gmx.gro"
    topology = outdir / "gmx.top"
    pdb = outdir / "gmx.pdb"
    itp = None
    if args.type == "parmed":
        cmd = [
            "amb2gro_top_gro.py",
            "-p",
            str(parm),
            "-c",
            str(rst),
            "-t",
            str(topology),
            "-g",
            str(gro),
            "-b",
            str(pdb),
        ]
        run_cmd(cmd, log="gmx.gro generated")
        LOGGER.info("gmx.top generated")
        LOGGER.info("gmx.pdb generated")
    # acpype may create wrong gro file
    # because of overflow of residue numbers or atom numbers
    elif args.type == "acpype":
        with tempfile.TemporaryDirectory(prefix="mdtbx-acpype-") as tempdir:
            cmd = [
                "acpype",
                "-p",
                str(parm),
                "-x",
                str(rst),
            ]
            run_cmd(cmd, cwd=tempdir)
            generated_dirs = sorted(Path(tempdir).glob("*.amb2gmx"))
            if len(generated_dirs) != 1:
                raise RuntimeError(
                    "ACPYPE must generate exactly one *.amb2gmx directory"
                )
            generated_dir = generated_dirs[0]
            generated_stem = generated_dir.name.removesuffix(".amb2gmx")
            shutil.copy2(generated_dir / f"{generated_stem}_GMX.gro", gro)
            source_topology = generated_dir / f"{generated_stem}_GMX.top"
            source_itp = generated_dir / f"{generated_stem}_GMX.itp"
            if source_itp.is_file():
                itp = outdir / "gmx.itp"
                shutil.copy2(source_itp, itp)
                topology.write_text(
                    source_topology.read_text().replace(source_itp.name, itp.name)
                )
            else:
                shutil.copy2(source_topology, topology)
        LOGGER.info("gmx.gro generated")
        LOGGER.info("gmx.top generated")

    if not args.no_editconf:
        # Keep the temp file next to the final gro so the replace below never
        # crosses filesystems.
        with gmx_tempfile(".gro", directory=outdir) as edited_path:
            cmd = [
                "gmx",
                "editconf",
                "-f",
                str(gro),
                "-o",
                edited_path,
                "-resnr",
                "1",
            ]
            run_cmd(cmd, log="gmx editconf residue renumbering completed")
            Path(edited_path).replace(gro)
        run_cmd(
            ["gmx", "editconf", "-f", str(gro), "-o", str(pdb)],
            log="gmx.pdb generated",
        )

    manifest = {
        "schema_version": 1,
        "workflow": "amber-to-gromacs",
        "conversion": args.type,
        "parm": str(parm),
        "rst": str(rst),
        "topology": str(topology),
        "structure": str(gro),
        "pdb": str(pdb) if pdb.is_file() else None,
        "itp": str(itp) if itp is not None else None,
    }
    (outdir / "conversion_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
