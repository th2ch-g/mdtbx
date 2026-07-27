import argparse
import math
import os
import re
from pathlib import Path

import mdtraj as md

from ..utils.common_args import add_topology_arg
from ..logger import generate_logger

LOGGER = generate_logger(__name__)

CONST_HATENA = 1
CONST_FUNCT = 1


def _ensure_list(values):
    if isinstance(values, (list, tuple)):
        return list(values)
    return [values]


def _expand_bounds(values, n_selection, label):
    values = _ensure_list(values)
    if len(values) == 1:
        return [values[0] for _ in range(n_selection)]
    if len(values) != n_selection:
        raise ValueError(f"number of selection and {label} should be same")
    return values


def _validate_bounds(lower_bounds, upper_bounds1, upper_bounds2):
    """Validate finite, ordered distance-restraint bounds."""
    for index, (lower, upper1, upper2) in enumerate(
        zip(lower_bounds, upper_bounds1, upper_bounds2, strict=True), start=1
    ):
        if not all(math.isfinite(value) for value in (lower, upper1, upper2)):
            raise ValueError(f"bounds for selection {index} must be finite")
        if lower < 0 or not lower <= upper1 <= upper2:
            raise ValueError(
                f"bounds for selection {index} must satisfy "
                "0 <= lower_bound <= upper_bound1 <= upper_bound2"
            )


def _macro_name(output_path):
    """Create a valid preprocessor symbol from an output filename."""
    name = re.sub(r"[^A-Za-z0-9_]", "_", output_path.stem).upper()
    if not name:
        raise ValueError("--output-prefix must contain a filename")
    if name[0].isdigit():
        name = f"_{name}"
    return name


def _append_include_once(topology_path, output_path):
    """Append a topology-relative include unless it already exists."""
    relative_path = os.path.relpath(
        output_path.resolve(), topology_path.parent.resolve()
    )
    include_line = f'#include "{Path(relative_path).as_posix()}"'
    topology_text = topology_path.read_text()
    if any(line.strip() == include_line for line in topology_text.splitlines()):
        return
    separator = "" if not topology_text or topology_text.endswith("\n") else "\n"
    with topology_path.open("a") as topology:
        topology.write(f"{separator}\n{include_line}\n")


def add_subcmd(subparsers):
    """
    mdtbx gen_distres -g structure.gro -p topology.top
    """
    parser = subparsers.add_parser(
        "gen_distres",
        help="Generate DISTRES",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-g", "--gro", required=True, type=str, help="GRO file (.gro)")

    add_topology_arg(parser, help="Topology file (.top)")

    parser.add_argument(
        "-s1",
        "--selection1",
        required=True,
        type=str,
        help="Selection for distance restraints. (MDtraj atom selection language). Use comma ',' if multiple selection",
    )

    parser.add_argument(
        "-s2",
        "--selection2",
        required=True,
        type=str,
        help="Selection for distance restraints. (MDtraj atom selection language). Use comma ',' if multiple selection",
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        "--output_prefix",
        dest="output_prefix",
        default="distres",
        type=str,
        help="Output file prefix (This also will be constant name)",
    )

    parser.add_argument(
        "-lo",
        "--lower-bound",
        "--lower_bound",
        dest="lower_bound",
        default=[0.0],
        type=float,
        nargs="*",
        help="Lower bound [nm]. Use multiple value if multiple selection. Single value will be applied to all selection if single value is given",
    )

    parser.add_argument(
        "-up1",
        "--upper-bound1",
        "--upper_bound1",
        dest="upper_bound1",
        default=[0.3],
        type=float,
        nargs="*",
        help="Upper bound1 [nm] Use multiple value if multiple selection. Single value will be applied to all selection if single value is given",
    )

    parser.add_argument(
        "-up2",
        "--upper-bound2",
        "--upper_bound2",
        dest="upper_bound2",
        default=[0.4],
        type=float,
        nargs="*",
        help="Upper bound2 [nm] Use multiple value if multiple selection. Single value will be applied to all selection if single value is given",
    )

    parser.set_defaults(func=run)


def run(args):
    # Parse paired atom selections.
    atom_selector1 = [s.strip() for s in args.selection1.split(",") if s.strip()]
    atom_selector2 = [s.strip() for s in args.selection2.split(",") if s.strip()]
    if not atom_selector1 or not atom_selector2:
        raise ValueError("selection1 and selection2 must not be empty")
    if len(atom_selector1) != len(atom_selector2):
        raise ValueError(
            "number of selection should be same for both selection1 and selection2"
        )

    lower_bound = _expand_bounds(args.lower_bound, len(atom_selector1), "lower bound")
    upper_bound1 = _expand_bounds(
        args.upper_bound1, len(atom_selector1), "upper bound1"
    )
    upper_bound2 = _expand_bounds(
        args.upper_bound2, len(atom_selector1), "upper bound2"
    )
    _validate_bounds(lower_bound, upper_bound1, upper_bound2)

    gro_path = Path(args.gro).expanduser()
    if not gro_path.is_file():
        raise FileNotFoundError(f"Structure file not found: {gro_path}")
    topology_path = Path(args.topology).expanduser()
    if not topology_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_path}")

    gro = md.load(str(gro_path))
    target_atom_indices1 = []
    target_atom_indices2 = []

    for selector1, selector2 in zip(atom_selector1, atom_selector2, strict=True):
        sele1 = gro.top.select(selector1)
        if len(sele1) != 1:
            raise ValueError(f"selection {selector1} should be single atom: {sele1}")
        sele2 = gro.top.select(selector2)
        if len(sele2) != 1:
            raise ValueError(f"selection {selector2} should be single atom: {sele2}")

        # convert 0-based mdtraj index to 1-based GROMACS atom id
        target_atom_indices1.append(sele1[0] + 1)
        target_atom_indices2.append(sele2[0] + 1)

    output_path = Path(f"{args.output_prefix}.itp").expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    const = _macro_name(output_path)
    force_const = const + "_FC"
    with output_path.open("w") as f:
        f.write(f"#ifdef {const}\n")
        f.write("[ intermolecular_interactions ]\n")
        f.write("[ distance_restraints ]\n")
        f.write(
            ";   i     j ? label      funct         lo        up1        up2     weight\n"
        )
        restraints = zip(
            target_atom_indices1,
            target_atom_indices2,
            lower_bound,
            upper_bound1,
            upper_bound2,
            atom_selector1,
            atom_selector2,
            strict=True,
        )
        for label, restraint in enumerate(restraints):
            atom1, atom2, lower, upper1, upper2, selector1, selector2 = restraint
            f.write(
                f"{atom1} {atom2} {CONST_HATENA} {label} {CONST_FUNCT} "
                f"{lower} {upper1} {upper2} {force_const} ; "
                f"'{selector1}'-'{selector2}'\n"
            )
        f.write("#endif\n")

    _append_include_once(topology_path, output_path)
