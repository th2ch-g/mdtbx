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
        help="Volume [A^3] (if this is provided, pdb is not required, directly calculate ion concentration)",
    )

    parser.add_argument(
        "-c",
        "--concentration",
        default=0.15,
        type=float,
        help="Ion concentration [M]",
    )

    parser.add_argument(
        "--water_name",
        default="WAT",
        type=str,
        help="Water name (e.g. WAT, HOH, TIP3)",
    )

    parser.add_argument(
        "--method",
        default="cubic",
        type=str,
        help="Method for calculating ion concentration (if pdb is provided)",
        choices=["cubic", "water", "optimize"],
    )

    # cubic: calculate from system volume(assume cubic system)
    # water: calculate from water volume(recomended if lipid system)
    # optimize: consider charge of system


def get_boxsize_from_pdb(args) -> tuple[float, float, float]:  # angstrom^3
    with open(args.pdb) as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if "CRYST" in line:
                parsed = re.findall(r"\S+", line)
                print(parsed)
                return float(parsed[1]), float(parsed[2]), float(parsed[3])
    raise Exception("CRYST line not found")


def get_water_number_from_pdb(args) -> int:
    count = 0
    with open(args.pdb) as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if "O" in line and args.water_name in line:
                count += 1
    return count


def calc_ion_conc_from_volume(volume: float, concentration: float) -> int:
    # concentration [M] = concentration [mol/L]
    # volume [A^3] = volume x 10^-3 [nm^3] = volume x 10^-3 x 10^-24 [L]
    # AVOGADRO = 6.022 * 10^23 [mol^-1]
    # -> concentration [M] * volume [A^3] * 1/10000
    ionnum = volume * concentration * AVOGADRO_CONST // 10000  # NOQA
    ionnum = int(ionnum)
    return ionnum


def run(args):
    if args.pdb is None and args.volume is None:
        raise Exception("pdb or volume is required")

    if args.volume is not None:
        ionnum = calc_ion_conc_from_volume(args.volume, args.concentration)
    else:
        if args.method == "cubic":
            # How to calculate ion concentration?
            # ref1: https://ambermd.org/tutorials/basic/tutorial8/index.php
            # ref2: https://github.com/YoshitakaMo/preparemd/blob/c6f21dd49e6a2aa0e3d12bd6ed60db37316baaa7/preparemd/amber/top/leapin.py#L33
            boxsize = get_boxsize_from_pdb(args)
            volume = boxsize[0] * boxsize[1] * boxsize[2]
            LOGGER.info(f"Volume of system: {volume}")
            ionnum = calc_ion_conc_from_volume(volume, args.concentration)
        elif args.method == "water":
            num_water = get_water_number_from_pdb(args)
            LOGGER.info(f"Number of water molecules: # {num_water}")
            volume = num_water * WATER_VOLUME * 1000  # NOQA # convert from nm^3 to A^3
            # volume = num_water * TIP3P_VOLUME * 1000
            LOGGER.info(f"Volume of water molecules: {volume}")
            ionnum = calc_ion_conc_from_volume(volume, args.concentration)
        else:
            raise NotImplementedError
    LOGGER.info(f"Number of ions that should be added: # {ionnum}")
    print(f"ionnum: {ionnum}")
