import argparse
import sys
import mdtraj as md
from itertools import groupby

from ..utils.proc import run_cmd
from ..utils.gmx import to_gmx_index
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

    parser.add_argument("-n", "--index", type=str, help="Input Index file")

    parser.add_argument(
        "-o", "--output", default="index.ndx", type=str, help="Output index file"
    )

    parser.set_defaults(func=run)


def make_default_index(args):
    cmd = f"echo q | gmx make_ndx -f {args.gro} -o {args.output}"
    run_cmd(cmd, log=f"{args.output} newly generated")


def count_index_group(args) -> int:
    count = 0
    with open(args.index) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("["):
                count += 1
    return count


def add_index(args):
    gro = md.load(args.gro)
    target_atom_indices = gro.topology.select(args.selection)
    target_atom_indices = to_gmx_index(target_atom_indices)  # index start from 0
    grouped_numbers = groupby(enumerate(target_atom_indices), lambda x: x[0] - x[1])
    result_parts = []
    for _, group in grouped_numbers:
        group_list = [item[1] for item in group]
        if len(group_list) > 1:
            result_parts.append(f"a {group_list[0]}-{group_list[-1]}")
        elif len(group_list) == 1:
            result_parts.append(f"a {group_list[0]}")
    if not result_parts:
        LOGGER.error(f"selection matched no atoms: {args.selection}")
        sys.exit(1)
    target_atom_indices = " | ".join(result_parts)
    LOGGER.info(f"Target selection in make_ndx: {target_atom_indices=}")
    count = count_index_group(args)
    LOGGER.info(f"Number of Current index group: {count=}")
    cmd = f"""
echo '
{target_atom_indices}
name {count} {args.name}
q
' | gmx make_ndx -f {args.gro} -n {args.index} -o {args.output}"""
    run_cmd(cmd, log=f"{args.output} updated")


def run(args):
    if args.selection is None and args.index is not None:
        LOGGER.error("--selection is required when --index is given")
        sys.exit(1)
    if args.index is None:
        make_default_index(args)
        args.index = args.output
    if args.selection is not None:
        add_index(args)
