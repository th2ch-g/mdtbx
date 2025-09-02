import argparse
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx shell_hook
    """
    _parser = subparsers.add_parser(
        "shell_hook",
        help="Generate shell hook",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )


def run(args):
    hook = f"""
----- BEGIN OF MDTBX SHELL HOOK -----
# for alias
shopt -s expand_aliases
alias mdtbx="pixi run --manifest-path {Path(__file__).parent.parent.parent} mdtbx"
alias pymol="pixi run --manifest-path {Path(__file__).parent.parent.parent} pymol"

# for function
mdtbx() {{
    pixi run --manifest-path {Path(__file__).parent.parent.parent} mdtbx "$@"
}}

pymol() {{
    pixi run --manifest-path {Path(__file__).parent.parent.parent} pymol "$@"
}}

# pymol template
pymol -c -p <<EOF
import pymol_plugins
from pymol_plugins import *
from pymol import cmd

EOF
----- END OF MDTBX SHELL HOOK -----
    """
    print(hook)
