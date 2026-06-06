import argparse
import re
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

    parser.set_defaults(func=run)


def _replace_atom_coords(line, new_crds):
    """Replace only the x/y/z coordinate fields of a mol2 ATOM line.

    Rewrites by token span from right to left so earlier spans stay valid,
    preserving the surrounding column layout and avoiding the fragile
    substring substitution (replace via dummy tokens) used previously.
    """
    matches = list(re.finditer(r"\S+", line))
    for ci in (4, 3, 2):
        m = matches[ci]
        line = line[: m.start()] + new_crds[ci - 2] + line[m.end() :]
    return line


def run(args):
    # First pass: map atom name -> xyz tokens from the coordinate mol2.
    atomname2crds = {}
    section = None
    with open(args.coordinates) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("@<TRIPOS>"):
                section = line
                continue
            if section == "@<TRIPOS>ATOM":
                parsed = re.findall(r"\S+", line)
                atomname2crds[parsed[1]] = parsed[2:5]

    # Second pass: copy the reference verbatim, swapping only ATOM coordinates.
    section = None
    new_lines = []
    with open(args.reference) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("@<TRIPOS>"):
                section = line
                new_lines.append(line)
                continue
            if section == "@<TRIPOS>ATOM":
                parsed = re.findall(r"\S+", line)
                new_lines.append(_replace_atom_coords(line, atomname2crds[parsed[1]]))
            else:
                new_lines.append(line)

    with open(args.output, "w") as f:
        f.write("\n".join(new_lines) + "\n")

    LOGGER.info(f"{args.output} updated")
