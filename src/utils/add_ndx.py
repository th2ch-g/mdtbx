import argparse
import subprocess
import mdtraj as md
from itertools import groupby

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx add_ndx
    """
    parser = subparsers.add_parser(
        "add_ndx",
        help="Add new index group with gmx make_ndx",
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

    parser.add_argument("--gmx", default="gmx", type=str, help="gmx command")

    parser.add_argument("-n", "--index", type=str, help="Input Index file")

    parser.add_argument(
        "-o", "--output", default="index.ndx", type=str, help="Output index file"
    )


def make_default_index(args):
    cmd = f"echo q | {args.gmx} make_ndx -f {args.gro} -o {args.output}"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.output} newly generated")


def add_index(args):
    gro = md.load(args.gro)
    target_atom_indices = gro.topology.select(args.selection)
    target_atom_indices = [i + 1 for i in target_atom_indices]  # index start from 0
    target_atom_indices = sorted(target_atom_indices)
    grouped_numbers = groupby(enumerate(target_atom_indices), lambda x: x[0] - x[1])
    result_parts = []
    for _, group in grouped_numbers:
        group_list = [item[1] for item in group]
        if len(group_list) > 1:
            result_parts.append(f"{group_list[0]}-{group_list[-1]}")
        else:
            result_parts.append(str(group_list[0]))
    target_atom_indices = " ".join(result_parts)
    print(target_atom_indices)
    count = count_index_group(args)
    cmd = f"""
echo '
a {target_atom_indices}
name {count + 1} {args.name}
q
' | {args.gmx} make_ndx -f {args.gro} -n {args.index} -o {args.output}"""
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
