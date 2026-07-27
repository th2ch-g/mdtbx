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
    parser.add_argument(
        "--match-by",
        choices=["name", "index"],
        default="name",
        help="Atom field used to match the coordinate and reference files",
    )

    parser.set_defaults(func=run)


def _replace_atom_coords(line, new_crds):
    """Replace only the x/y/z coordinate fields of a mol2 ATOM line.

    Rewrites by token span from right to left so earlier spans stay valid,
    preserving the surrounding column layout and avoiding the fragile
    substring substitution (replace via dummy tokens) used previously.
    """
    matches = list(re.finditer(r"\S+", line))
    if len(matches) < 5:
        raise ValueError(f"Malformed MOL2 atom line: {line!r}")
    for ci in (4, 3, 2):
        m = matches[ci]
        line = line[: m.start()] + new_crds[ci - 2] + line[m.end() :]
    return line


def run(args):
    match_by = getattr(args, "match_by", "name")
    key_index = 1 if match_by == "name" else 0

    # First pass: map the requested atom field to xyz tokens.
    atomkey2crds = {}
    section = None
    with open(args.coordinates) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("@<TRIPOS>"):
                section = line
                continue
            if section == "@<TRIPOS>ATOM":
                parsed = re.findall(r"\S+", line)
                if not parsed:
                    continue
                if len(parsed) < 5:
                    raise ValueError(f"Malformed MOL2 atom line: {line!r}")
                atom_key = parsed[key_index]
                if atom_key in atomkey2crds:
                    raise ValueError(
                        f"Duplicate atom {match_by} in coordinate MOL2: {atom_key}"
                    )
                atomkey2crds[atom_key] = parsed[2:5]

    # Second pass: copy the reference verbatim, swapping only ATOM coordinates.
    section = None
    new_lines = []
    reference_atom_keys = set()
    with open(args.reference) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("@<TRIPOS>"):
                section = line
                new_lines.append(line)
                continue
            if section == "@<TRIPOS>ATOM":
                parsed = re.findall(r"\S+", line)
                if not parsed:
                    new_lines.append(line)
                    continue
                if len(parsed) < 5:
                    raise ValueError(f"Malformed MOL2 atom line: {line!r}")
                atom_key = parsed[key_index]
                if atom_key in reference_atom_keys:
                    raise ValueError(
                        f"Duplicate atom {match_by} in reference MOL2: {atom_key}"
                    )
                if atom_key not in atomkey2crds:
                    raise ValueError(
                        f"Atom {match_by} {atom_key} is missing from coordinate MOL2"
                    )
                reference_atom_keys.add(atom_key)
                new_lines.append(_replace_atom_coords(line, atomkey2crds[atom_key]))
            else:
                new_lines.append(line)

    extra_atom_keys = set(atomkey2crds) - reference_atom_keys
    if extra_atom_keys:
        extras = ", ".join(sorted(extra_atom_keys))
        raise ValueError(
            f"Coordinate MOL2 has atom {match_by}s missing from reference MOL2: {extras}"
        )

    with open(args.output, "w") as f:
        f.write("\n".join(new_lines) + "\n")

    LOGGER.info(f"{args.output} updated")
