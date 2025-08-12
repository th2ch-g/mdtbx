import argparse
import subprocess
import mdtraj as md

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx add_make_ndx
    """
    parser = subparsers.add_parser(
        "add_make_ndx",
        help="Add new index group",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-g", "--gro", required=True, type=str, help="Gro file (.gro)")

    parser.add_argument(
        "-s",
        "--selection",
        type=str,
        help="Selection for new index group (MDtraj atom selection language)",
    )

    parser.add_argument(
        "--name", default="NEWGROUP", type=str, help="Name for new index group"
    )

    parser.add_argument("-g", "--gmx", default="gmx", type=str, help="gmx command")

    parser.add_argument("-n", "--index", type=str, help="Input Index file")

    parser.add_argument(
        "-o", "--output", default="index.ndx", type=str, help="Output index file"
    )


def make_default_index(args):
    cmd = f"echo q | {args.gmx} make_ndx -f {args.structure} -o {args.index}"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.index} newly generated")


def add_index(args):
    gro = md.load(args.gro)
    target_atom_indices = gro.topology.select(args.selection)
    target_atom_indices = [i + 1 for i in target_atom_indices]  # index start from 0
    target_atom_indices = " ".join([str(i) for i in target_atom_indices])
    count = count_index_group(args)
    cmd = f"""
echo '
a {target_atom_indices}
name {count + 1} {args.name}
q
' | {args.gmx} make_ndx -f {args.structure} -n {args.index} -o {args.output}"""
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.output} updated")


def count_index_group(args) -> int:
    count = 0
    with open(args.index) as f:
        for idx, line in enumerate(f):
            line = line.rstrip()
            if line.startswith("["):
                count += 1
    return count


def run(args):
    if args.selection is None and args.index is not None:
        LOGGER.warn("use --selection")
    if args.selection is None and args.index is None:
        make_default_index(args)
        args.index = args.output
    if args.selection is not None and args.index is None:
        make_default_index(args)
        args.index = args.output
        add_index(args)
    if args.selection is not None and args.index is not None:
        add_index(args)
