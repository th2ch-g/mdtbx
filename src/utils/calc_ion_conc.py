import argparse
import re

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx calc_ion_conc -p structure.pdb
    """
    parser = subparsers.add_parser(
        "calc_ion_conc",
        help="Calculate ion concentration",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--pdb",
        type=str,
        help="PDB file (.pdb)",
    )

    parser.add_argument(
        "-v",
        "--volume",
        type=float,
        help="Volume [A^3]",
    )

    parser.add_argument(
        "-c",
        "--concentration",
        default=0.15,
        type=float,
        help="Ion concentration [M]",
    )


def get_boxsize_from_pdb(pdb_file) -> tuple[float, float, float]:  # angstrom^3
    with open(pdb_file) as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if "CRYST" in line:
                parsed = re.findall(r"\S+", line[0])
                return float(parsed[1]), float(parsed[2]), float(parsed[3])
    raise Exception("CRYST line not found")


def run(args):
    # How to calculate ion concentration?
    # ref1: https://ambermd.org/tutorials/basic/tutorial8/index.php
    # ref2: https://github.com/YoshitakaMo/preparemd/blob/c6f21dd49e6a2aa0e3d12bd6ed60db37316baaa7/preparemd/amber/top/leapin.py#L33

    if args.volume is None and args.pdb is not None:
        boxsize = get_boxsize_from_pdb(args.pdb)
        volume = boxsize[0] * boxsize[1] * boxsize[2]
    if args.volume is not None and args.pdb is None:
        volume = args.volume

    if args.pdb is None and args.volume is None:
        raise Exception("pdb or volume is required")

    # concentration [M] = concentration [mol/L]
    # volume [A^3] = volume x 10^-3 [nm^3] = volume x 10^-3 x 10^-24 [L]
    # AVOGADRO = 6.022 * 10^23 [mol^-1]
    # -> concentration [M] * volume [A^3] * 1/10000
    ionnum = volume * args.concentration * AVOGADRO_CONST // 10000 # NOQA
    ionnum = int(ionnum)
    LOGGER.info(f"Number of ions that should be added: # {ionnum}")
    print(f"ionnum: {ionnum}")
