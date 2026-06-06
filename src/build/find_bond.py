import argparse

from pymol import cmd

from ..config import SYSTEM_NAME
from ..logger import generate_logger
from ..utils.pymol_session import pymol_session

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx find_bond -s structure -s1 selection -s2 selection -c cutoff -o output
    """
    parser = subparsers.add_parser(
        "find_bond",
        help="Find bond (e.g. CYS-CYS)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Structure file (can be parsed in PyMOL)",
    )

    parser.add_argument(
        "-s1",
        "--selection1",
        type=str,
        default="resn CYS and name SG",
        help="Selection 1 (can be parsed in PyMOL)",
    )

    parser.add_argument(
        "-s2",
        "--selection2",
        default="resn CYS and name SG",
        type=str,
        help="Selection 2 (can be parsed in PyMOL)",
    )

    parser.add_argument(
        "-c",
        "--cutoff",
        default=3.0,
        type=float,
        help="Cutoff distance [angstrom]",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name",
    )

    parser.add_argument(
        "-op",
        "--output-pdb",
        type=str,
        help="Output PDB file with CYS renamed to CYM for SS-bonded residues",
    )

    parser.set_defaults(func=run)


def run(args):
    with pymol_session(cmd, args.structure):
        pairs = cmd.find_pairs(args.selection1, args.selection2, cutoff=args.cutoff)
        bonds_set = set()
        # Keep the bonded SG atom indices together with their resi so the CYM rename
        # can target each residue uniquely (atom index), avoiding cross-chain CYS
        # with the same residue number being renamed by accident.
        bonds = []
        for idx1, idx2 in pairs:
            if idx1[1] in bonds_set or idx2[1] in bonds_set:
                continue
            bonds_set.add(idx1[1])
            bonds_set.add(idx2[1])
            tmp = []
            cmd.iterate_state(1, f"index {idx1[1]}", "tmp.append(resi)", space=locals())
            cmd.iterate_state(1, f"index {idx2[1]}", "tmp.append(resi)", space=locals())
            bonds.append(
                {
                    "index1": idx1[1],
                    "index2": idx2[1],
                    "resi1": tmp[0],
                    "resi2": tmp[1],
                }
            )
        bonds_str_lines = []
        for bond in bonds:
            bonds_str_lines.append(
                f"bond {SYSTEM_NAME}.{bond['resi1']}.SG {SYSTEM_NAME}.{bond['resi2']}.SG"  # NOQA
            )
        bonds_str = "\n".join(bonds_str_lines)

        if args.output is not None:
            with open(args.output, "w") as f:
                f.write(bonds_str)
            LOGGER.info(f"{args.output} generated")
        else:
            if len(bonds) == 0:
                LOGGER.info("No bond found")
            else:
                LOGGER.info(f"{len(bonds)} bonds found")
                print(bonds_str)

        if args.output_pdb is not None:
            if len(bonds) == 0:
                LOGGER.info("No SS-bond found; saving original structure to output PDB")
                cmd.save(args.output_pdb, "target")
            else:
                # Rename only the residues that actually own a bonded SG atom.
                # Target by the SG atom index (byres) so a CYS in another chain that
                # happens to share the residue number is left untouched.
                bonded_indices = set()
                for bond in bonds:
                    bonded_indices.add(bond["index1"])
                    bonded_indices.add(bond["index2"])
                for idx in sorted(bonded_indices):
                    cmd.alter(
                        f"target and resn CYS and byres (index {idx})", "resn='CYM'"
                    )
                    LOGGER.info(f"CYS SG (index {idx}) residue renamed to CYM")
                cmd.save(args.output_pdb, "target")
                LOGGER.info(f"{args.output_pdb} generated")
