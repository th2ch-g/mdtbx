import argparse
import json
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

    parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        type=str,
        help="Output directory",
    )

    parser.set_defaults(func=run)


def run(args):
    structure = Path(args.structure).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    mol2 = outdir / f"{args.resname}.mol2"
    frcmod = outdir / f"{args.resname}.frcmod"
    library = outdir / f"{args.resname}.lib"
    pdb = outdir / f"{args.resname}.pdb"
    filetype = structure.suffix[1:]
    # antechamber uses "mdl" as the format flag for MDL .mol files
    if filetype == "mol":
        filetype = "mdl"
    cmd = [
        "antechamber",
        "-i",
        str(structure),
        "-fi",
        filetype,
        "-o",
        str(mol2),
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
    run_cmd(cmd, log=f"{mol2} generated", cwd=outdir)

    cmd = [
        "parmchk2",
        "-i",
        str(mol2),
        "-f",
        "mol2",
        "-o",
        str(frcmod),
        "-s",
        "gaff2",
    ]
    run_cmd(cmd, log=f"{frcmod} generated", cwd=outdir)

    cmd = [
        "antechamber",
        "-i",
        str(mol2),
        "-fi",
        "mol2",
        "-o",
        str(pdb),
        "-fo",
        "pdb",
        "-pf",
        "y",
    ]
    run_cmd(cmd, log=f"{pdb} generated", cwd=outdir)

    cmd_tleap = f"""
source leaprc.gaff2
loadamberparams {frcmod}
{args.resname} = loadmol2 {mol2}
saveoff {args.resname} {library}
quit
    """
    run_tleap(cmd_tleap, cwd=outdir)
    manifest = {
        "schema_version": 1,
        "workflow": "gaff2-am1bcc",
        "structure": str(structure),
        "resname": args.resname,
        "charge": args.charge,
        "multiplicity": args.multiplicity,
        "mol2": str(mol2),
        "frcmod": str(frcmod),
        "library": str(library),
        "pdb": str(pdb),
    }
    (outdir / "parameterization_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    LOGGER.info(f"{library} and {pdb} generated")
