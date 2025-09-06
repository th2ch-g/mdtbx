import argparse
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx mod_mdp
    """
    parser = subparsers.add_parser(
        "mod_mdp",
        help="Modify mdp files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--path", type=str, help="Path to the directory", default=".")

    parser.add_argument(
        "-t",
        "--target_variable",
        type=str,
        required=True,
        help="Target variable. If variable is not found, it will be added.",
    )

    parser.add_argument("-v", "--new_value", type=str, required=True, help="New value")

    parser.add_argument("--exclude", type=str, nargs="+", help="Exclude name of files")

    parser.add_argument(
        "-lj", "--ljust", type=int, default=23, help="Ljust for new variable line"
    )


def run(args):
    for mdp in Path(args.path).glob("*.mdp"):
        key_skip = False
        if args.exclude is not None:
            for exclude in args.exclude:
                if exclude in mdp.name:
                    LOGGER.info(f"{mdp} excluded")
                    key_skip = True
                    continue
        if key_skip:
            continue
        mod_mdp(args.target_variable, args.new_value, mdp, args.ljust)


def mod_mdp(target_variable, new_value, mdp, ljust):
    new_lines = []
    added_key = False
    with open(mdp) as f:
        for idx, line in enumerate(f):
            line_ = line.strip()
            if line_.startswith(";"):
                new_lines.append(line)
                continue
            if line_.startswith(target_variable):
                current_value = line_.split("=")[1].strip()
                new_line = line.replace(current_value, new_value)
                new_lines.append(new_line)
                added_key = True
            else:
                new_lines.append(line)

    if not added_key:
        new_lines.append(f"{target_variable.ljust(ljust)} = {new_value}\n")
        LOGGER.info(f"new variable {target_variable} added to {mdp}")
    else:
        LOGGER.info(f"{mdp} modified")

    with open(mdp, "w") as f:
        f.writelines(new_lines)
