import argparse
import mdtraj as md

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)

CONST_HATENA = 1
CONST_FUNCT = 1


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

    parser.add_argument(
        "-p", "--topology", required=True, type=str, help="Topology file (.top)"
    )

    parser.add_argument(
        "-s1",
        "--selection1",
        required=True,
        type=str,
        help="Selection for distance restraints. (MDtraj atom selection language). Use camma ',' if multiple selection",
    )

    parser.add_argument(
        "-s2",
        "--selection2",
        required=True,
        type=str,
        help="Selection for distance restraints. (MDtraj atom selection language). Use camma ',' if multiple selection",
    )

    parser.add_argument(
        "-o",
        "--output_prefix",
        default="distres",
        type=str,
        help="Output file prefix (This also will be constant name)",
    )

    parser.add_argument(
        "-lo", "--lower_bound", default=0.0, type=float, help="Lower bound [nm]"
    )

    parser.add_argument(
        "-up1", "--upper_bound1", default=0.3, type=float, help="Upper bound1 [nm]"
    )

    parser.add_argument(
        "-up2", "--upper_bound2", default=0.4, type=float, help="Upper bound2 [nm]"
    )


def run(args):
    # generate posres.itp
    atom_selector1 = args.selection1.split(",")
    atom_selector2 = args.selection2.split(",")
    assert len(atom_selector1) == len(atom_selector2), (
        "number of selection should be same for both selection1 and selection2"
    )

    gro = md.load(args.gro)
    target_atom_indices1 = []
    target_atom_indices2 = []

    for i in range(len(atom_selector1)):
        sele1 = gro.top.select(atom_selector1[i])
        assert len(sele1) == 1, (
            f"selection {atom_selector1[i]} should be single atom: {sele1}"
        )
        sele2 = gro.top.select(atom_selector2[i])
        assert len(sele2) == 1, (
            f"selection {atom_selector2[i]} should be single atom: {sele2}"
        )

        # index start from 0
        target_atom_indices1.append(sele1[0] + 1)
        target_atom_indices2.append(sele2[0] + 1)

    const = args.output_prefix.upper()
    force_const = const + "_FC"
    with open(args.output_prefix + ".itp", "w") as f:
        f.write(f"#ifdef {const}\n")
        f.write("[ intermolecular_interactions ]\n")
        f.write("[ distance_restraints ]\n")
        f.write(
            ";   i     j ? label      funct         lo        up1        up2     weight\n"
        )
        for i in range(len(target_atom_indices1)):
            f.write(
                f"{target_atom_indices1[i]} {target_atom_indices2[i]} {CONST_HATENA} {i} {CONST_FUNCT} {args.lower_bound} {args.upper_bound1} {args.upper_bound2} {force_const}\n"
            )
        f.write("#endif\n")

    # insert distres.itp into topology.top
    # system section treats as global
    with open(args.topology, "a") as f:
        f.write(f'\n#include "{args.output_prefix}.itp"\n')
