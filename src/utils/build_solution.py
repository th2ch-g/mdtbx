import argparse
from pathlib import Path
from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx build_solution -f input_structure.pdb -o output_structure.pdb --ion_conc 0.15 --cation Na+ --anion Cl- --ligparam FRCMOD:LIB --boxsize 100 100 100
    """
    parser = subparsers.add_parser(
        "build_solution",
        help="Build solution system",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-f",
        "--input",
        required=True,
        type=str,
        help="Input Structure file (.pdb)",
    )

    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Output file path"
    )

    parser.add_argument(
        "--ion_conc", default=0.15, type=float, help="Ion concentration [M]"
    )

    parser.add_argument(
        "--cation", default="Na+", type=str, help="Cation name"
    )

    parser.add_argument(
        "--anion", default="Cl-", type=str, help="Anion name"
    )

    parser.add_argument(
        "--ligparam", type=str, help="Ligand parameter"
    )

    parser.add_argument(
        "--boxsize", nargs=3, type=float, help="Box size [angstrom, angstrom, angstrom]"
    )


def run(args):
    # tleap
    with open(Path(__file__).parent / "template_tleap.in") as ref:
        lines = ref.readlines()

    # ION_CONC, INPUT_PDB, LIGAND_PARAMS, SSBONDS


    pass
