import argparse
import subprocess
from pathlib import Path

from ..logger import generate_logger

LOGGER = generate_logger(__name__)

ONE_TO_THREE = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "E": "GLU",
    "Q": "GLN",
    "G": "GLY",
    "H": "HIE",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

TLEAP_TEMPLATE = """\
source leaprc.protein.{ff}
mol = sequence {{ {residues} }}
savepdb mol {output}
quit
"""


def _to_three_letter(sequence: str) -> list[str]:
    """Convert one-letter or space-separated three-letter sequence to list of three-letter codes."""
    if " " in sequence:
        return sequence.upper().split()
    return [ONE_TO_THREE[aa] for aa in sequence.upper()]


def add_subcmd(subparsers):
    """
    mdtbx gen_loop_aa ACGST -o loop.pdb
    mdtbx gen_loop_aa "ALA CYS GLY SER THR" -o loop.pdb
    """
    parser = subparsers.add_parser(
        "gen_loop_aa",
        help="Generate linear (loop) peptide structure from amino acid sequence using tleap",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "sequence",
        type=str,
        help="Amino acid sequence: one-letter codes (e.g. ACGST) or space-separated three-letter codes (e.g. 'ALA CYS GLY')",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="loop.pdb",
        type=str,
        help="Output PDB file",
    )

    parser.add_argument(
        "--ff",
        default="ff14SB",
        choices=["ff14SB", "ff19SB"],
        help="Force field for tleap",
    )

    parser.add_argument(
        "--cap",
        action="store_true",
        help="Add ACE/NME capping groups to N/C termini",
    )

    parser.add_argument(
        "--keepfiles",
        action="store_true",
        help="Keep intermediate files (tleap.in, leap.log)",
    )

    parser.set_defaults(func=run)


def run(args):
    residues = _to_three_letter(args.sequence)

    if args.cap:
        residues = ["ACE"] + residues + ["NME"]

    tleap_input = TLEAP_TEMPLATE.format(
        ff=args.ff,
        residues=" ".join(residues),
        output=args.output,
    )

    Path("tleap.in").write_text(tleap_input)
    subprocess.run("tleap -f tleap.in", shell=True, check=True)
    LOGGER.info(f"{args.output} generated")

    if not args.keepfiles:
        subprocess.run("rm -f leap.log tleap.in", shell=True, check=True)
        LOGGER.info("leap.log tleap.in removed")
