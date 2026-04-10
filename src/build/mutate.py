import argparse
from pathlib import Path

from pymol import cmd as pymol_cmd

from ..logger import generate_logger

LOGGER = generate_logger(__name__)

SUPPORTED_MUTANTS = {
    "ALA",
    "ARG",
    "ARGN",
    "ASN",
    "ASP",
    "ASPH",
    "CYS",
    "GLN",
    "GLU",
    "GLUH",
    "GLY",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYSN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def add_subcmd(subparsers):
    """
    mdtbx mutate -s input.pdb --selection "chain A and resi 42" --mutant VAL -o mutated.pdb
    """
    parser = subparsers.add_parser(
        "mutate",
        help="Mutate one residue with the PyMOL mutagenesis wizard",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Input structure file (PDB etc.)",
    )
    parser.add_argument(
        "--selection",
        required=True,
        type=str,
        help="PyMOL atom selection for the target residue",
    )
    parser.add_argument(
        "-m",
        "--mutant",
        required=True,
        type=str,
        help="Target residue name (e.g. VAL, HIE, HIP)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="mutant.pdb",
        type=str,
        help="Output structure file",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        help="Deprecated alias for output prefix. Writes <prefix>.pdb",
    )
    parser.add_argument(
        "--rotamer",
        default=1,
        type=int,
        help="Rotamer index (1-based) shown by the mutagenesis wizard",
    )
    parser.add_argument(
        "--hydrogens",
        default="auto",
        choices=["auto", "keep", "none"],
        help="Hydrogen handling mode in the mutagenesis wizard",
    )

    parser.set_defaults(func=run)


def _count_selected_residues(selection_name: str) -> int:
    space = {"residues": set()}
    pymol_cmd.iterate(
        f"byres ({selection_name})",
        "residues.add((model, segi, chain, resi, resn))",
        space=space,
    )
    return len(space["residues"])


def run(args):
    mutant = args.mutant.upper()
    if mutant not in SUPPORTED_MUTANTS:
        raise ValueError(f"Unsupported mutant residue: {mutant}")
    if args.rotamer < 1:
        raise ValueError("--rotamer must be 1 or larger")

    structure_path = Path(args.structure)
    if args.output_prefix is not None:
        output_pdb = Path(f"{args.output_prefix}.pdb")
    else:
        output_pdb = Path(args.output)
    object_name = "target"
    selection_name = "_mutate_input"

    pymol_cmd.reinitialize()
    pymol_cmd.load(str(structure_path), object_name)
    pymol_cmd.select(selection_name, f"({object_name}) and ({args.selection})")

    residue_count = _count_selected_residues(selection_name)
    if residue_count != 1:
        raise ValueError(
            f"Selection must match exactly one residue, but got {residue_count}: {args.selection}"
        )

    pymol_cmd.wizard("mutagenesis")
    wizard = pymol_cmd.get_wizard()
    wizard.set_mode(mutant)
    wizard.set_hyd(args.hydrogens)
    wizard.do_select(selection_name)
    if args.rotamer > 1:
        pymol_cmd.frame(args.rotamer)
    wizard.apply()
    pymol_cmd.set_wizard()
    pymol_cmd.save(str(output_pdb), object_name)

    LOGGER.info(f"{output_pdb} generated")
