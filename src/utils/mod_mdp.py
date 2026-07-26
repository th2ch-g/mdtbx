import argparse
from pathlib import Path

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

    parser.set_defaults(func=run)


def run(args):
    excluded_names = args.exclude or []
    for mdp in sorted(Path(args.path).glob("*.mdp")):
        if any(exclude in mdp.name for exclude in excluded_names):
            LOGGER.info("%s excluded", mdp)
            continue
        mod_mdp(args.target_variable, args.new_value, mdp, args.ljust)


def mod_mdp(
    target_variable: str,
    new_value: str,
    mdp: str | Path,
    ljust: int,
) -> None:
    if (
        not target_variable
        or target_variable != target_variable.strip()
        or any(char in target_variable for char in "=\r\n;")
    ):
        raise ValueError("Target variable must be a non-empty MDP key")
    if "\n" in new_value or "\r" in new_value:
        raise ValueError("New value must be a single line")
    if ljust < 0:
        raise ValueError("ljust must be non-negative")

    mdp_path = Path(mdp)
    new_lines: list[str] = []
    added_key = False
    with mdp_path.open() as f:
        for line in f:
            setting, comment_separator, comment = line.partition(";")
            stripped_setting = setting.strip()
            if not stripped_setting:
                new_lines.append(line)
                continue
            key, sep, _rest = setting.partition("=")
            if sep and key.strip() == target_variable:
                new_line = f"{key.rstrip()} = {new_value}\n"
                if comment_separator:
                    new_line = f"{key.rstrip()} = {new_value} ;{comment.rstrip()}\n"
                new_lines.append(new_line)
                added_key = True
            else:
                new_lines.append(line)

    if not added_key:
        if new_lines and not new_lines[-1].endswith("\n"):
            new_lines[-1] += "\n"
        new_lines.append(f"{target_variable.ljust(ljust)} = {new_value}\n")
        LOGGER.info("New variable %s added to %s", target_variable, mdp_path)
    else:
        LOGGER.info("%s modified", mdp_path)

    with mdp_path.open("w") as f:
        f.writelines(new_lines)
