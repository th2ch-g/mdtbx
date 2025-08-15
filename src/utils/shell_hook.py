import argparse
import subprocess
from pathlib import Path

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx shell_hook
    """
    parser = subparsers.add_parser(
        "shell_hook",
        help="Generate shell hook",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

def run(args):
    hook = f"""
# begin of mdtbx shell hook
alias mdtbx = "pixi run --manifest-path {Path(__file__).parent.parent.parent} mdtbx"
alias pymol = "pixi run --manifest-path {Path(__file__).parent.parent.parent} pymol"

pymol -c -p <<EOF
import myplugins
from pymol import cmd
EOF

# end of mdtbx shell hook
    """
    print(hook)
