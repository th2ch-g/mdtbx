import argparse

from pymol import cmd, editor

from ..logger import generate_logger
from ..utils.pymol_session import pymol_session

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "addnme",
        help="Add NME to protein",
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
        "-o", "--output_prefix", default="out_nme", type=str, help="Output file prefix"
    )

    parser.set_defaults(func=run)


def run(args):
    output = f"{args.output_prefix}.pdb"
    with pymol_session(cmd, args.structure):
        cmd.select("oxts", "name OXT")
        cmd.remove("oxts")
        LOGGER.info("OXTs removed")
        for chain in cmd.get_chains("target and polymer.protein"):
            if chain:
                selection = f"last (chain {chain}) and name C"
            else:
                selection = "last polymer.protein and name C"
            editor.attach_amino_acid(selection, "nme")
            LOGGER.info(f"NME added to chain '{chain}'")
        cmd.set("retain_order", 0)
        cmd.sort()
        cmd.save(output)

    # Post-process the PDB once in memory: read -> edit the list -> write once.
    # Re-enumerating the file on disk while insert/pop-ing the in-memory list
    # desynchronizes indices, so `lines` is the single mutable source of truth.
    with open(output) as f:
        lines = f.readlines()

    # 1) Rename NME cap atom names (HH3x -> Hx, CH3 -> C), preserving columns.
    for i, line in enumerate(lines):
        if "NME" in line:
            for k in range(1, 4):
                lines[i] = lines[i].replace(f"HH3{k}", f" H{k} ", 1)
            lines[i] = lines[i].replace("CH3", "C  ", 1)

    # 2) Insert a TER after each NME residue's last hydrogen (H3). Walk in
    #    reverse so earlier insertions do not shift later target indices.
    for i in range(len(lines) - 1, -1, -1):
        if (
            lines[i].startswith(("ATOM", "HETATM"))
            and "NME" in lines[i]
            and lines[i][12:16].strip() == "H3"
        ):
            lines.insert(i + 1, "TER\n")

    # 3) Drop a stray TER immediately before an NME backbone nitrogen so the cap
    #    stays bonded to the preceding residue. Match the atom-name column: a bare
    #    `"N" in line` is always true because the residue name NME contains 'N'.
    for i in range(len(lines) - 1, 0, -1):
        if (
            lines[i].startswith(("ATOM", "HETATM"))
            and "NME" in lines[i]
            and lines[i][12:16].strip() == "N"
            and lines[i - 1].startswith("TER")
        ):
            lines.pop(i - 1)

    with open(output, "w") as f:
        f.writelines(lines)
