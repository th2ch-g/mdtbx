import argparse
import re
from pathlib import Path
from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx mv_crds_mol2 -r reference.mol2 -c coordinates.mol2
    """
    parser = subparsers.add_parser(
        "mv_crds_mol2",
        help="Move coordinates to mol2 file (coordinates -> reference)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        type=str,
        help="Reference mol2 file",
    )

    parser.add_argument(
        "-c",
        "--coordinates",
        required=True,
        type=str,
        help="Coordinates mol2 file",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="output.mol2",
        type=str,
        help="Output mol2 file",
    )

def run(args):
    atomname2crds = {}

    section = None
    with open(args.coordinates) as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if line.startwith("@TRIPOS>"):
                section = line
                continue
            if section == "@TRIPOS>ATOM":
                parsed = re.findall(r"\S+", line)
                atomname = parsed[1]
                crds = parsed[2:5]
                atomname2crds[atomname] = crds

    section = None
    dummy = ["XXX", "YYY", "ZZZ"]
    new_lines = []
    with open(args.reference) as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if line.startwith("@TRIPOS>"):
                section = line
                new_lines.append(line)
                continue
            if section == "@TRIPOS>ATOM":
                parsed = re.findall(r"\S+", line)
                atomname = parsed[1]
                old_crds = parsed[2:5]
                new_crds = atomname2crds[atomname]
                for i in range(len(old_crds)):
                    line = line.replace(old_crds[i], dummy[i], 1)
                for i in range(len(old_crds)):
                    line = line.replace(dummy[i], new_crds[i], 1)
                new_lines.append(line)

    with open(args.output, "w") as f:
        f.writelines(new_lines)

    LOGGER.info(f"{args.output} updated")
